#Scatter plot for VCF
#!/usr/bin/env Rscript

# ===============================
# SNV Annotation & Visualization
# ===============================

suppressPackageStartupMessages({
  library(data.table)
  library(VariantAnnotation)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(IRanges)
  library(ggplot2)
  library(viridis)
  library(txdbmaker)  # <- makeTxDbFromGFF moved here
})

# --- Inputs (edit as needed) ---
snv_dir  <- "/cluster/home/t922316uhn/output/GATK"  # directory with VCFs
gtf_file <- "/cluster/tools/data/genomes/human/hg38/iGenomes/Annotation/Genes/genes.gtf" # GTF path

# --------------------------------
# Helpers
# --------------------------------
# make chromosome style match "chr1..chr22,chrX,chrY,chrM"
normalize_seqlevels_chr <- function(x) {
  x <- as.character(x)
  has_chr <- any(grepl("^chr", x))
  if (!has_chr) {
    x <- paste0("chr", x)
  }
  # Common MT naming fixes
  x <- sub("^chrMT$", "chrM", x, perl = TRUE)
  x <- sub("^chrM[T]?$", "chrM", x, perl = TRUE)
  x
}

# Extract attribute value from GTF attribute column; returns NA if missing
attr_get <- function(attr, key) {
  # pattern like key "VALUE";
  m <- regexpr(paste0("(^|;\\s*)", key, "\\s+\"([^\"]+)\""), attr, perl = TRUE)
  out <- rep(NA_character_, length(attr))
  hits <- m > 0
  if (any(hits)) {
    ext <- regmatches(attr, m)
    out[hits] <- sub(paste0(".*", key, '\\s+"([^"]+)".*'), "\\1", ext)
  }
  out
}

# --------------------------------
# Parse VCFs → data.table (sample, chr, pos)
# --------------------------------
parse_snv_for_scatter <- function(dir){
  files <- list.files(dir, pattern="\\.vcf(\\.gz)?$", full.names=TRUE)
  if(length(files) == 0) stop("No VCF files found in directory: ", dir)

  out <- vector("list", length(files))
  names(out) <- basename(files)

  for(i in seq_along(files)){
    f <- files[i]
    sample <- sub("\\.vcf(\\.gz)?$", "", basename(f))

    message("Parsing SNV VCF: ", basename(f))
    vcf <- tryCatch(
      readVcf(f, "hg38"),
      error = function(e) { message("  - Skipping (read error): ", e$message); NULL }
    )
    if (is.null(vcf) || length(vcf) == 0) next

    vr <- rowRanges(vcf)
    dt <- data.table(
      sample = sample,
      chr    = as.character(seqnames(vr)),
      pos    = as.integer(start(vr))
    )
    out[[i]] <- dt
  }

  res <- rbindlist(out, use.names=TRUE, fill=TRUE)
  if (nrow(res) == 0) stop("No variants parsed from VCFs in: ", dir)
  # normalize chromosome names to 'chr*'
  res[, chr := normalize_seqlevels_chr(chr)]
  res
}

# --------------------------------
# Robust GTF reader → gene GRanges with gene_name fallback
# Works with GTFs that only have exons (aggregates to gene span)
# --------------------------------
gtf_to_gene_granges <- function(gtf_file) {
  # Read GTF without using comment.char (not supported by fread on some builds)
  gtf_dt <- fread(
    gtf_file,
    header = FALSE, sep = "\t", quote = "",
    col.names = c("seqname","source","feature","start","end","score","strand","frame","attribute"),
    data.table = TRUE, showProgress = FALSE
  )

  # Drop comment lines (start with '#')
  gtf_dt <- gtf_dt[!startsWith(seqname, "#")]

  # Extract IDs/names from attributes
  gtf_dt[, gene_id   := attr_get(attribute, "gene_id")]
  gtf_dt[, gene_name := attr_get(attribute, "gene_name")]

  # If we have explicit gene rows, use them; else aggregate exons
  if ("gene" %in% gtf_dt$feature) {
    genes_dt <- gtf_dt[feature == "gene", .(
      seqname, start = as.integer(start), end = as.integer(end),
      strand, gene_id, gene_name
    )]
  } else {
    # Need exons to aggregate
    if (!("exon" %in% gtf_dt$feature)) {
      stop("GTF does not contain 'gene' or 'exon' features; cannot derive gene ranges.")
    }
    exons <- gtf_dt[feature == "exon" & !is.na(gene_id)]
    if (nrow(exons) == 0) stop("No exon rows with gene_id found in GTF.")
    # Aggregate exon mins/max to gene span
    genes_dt <- exons[, .(
      start = as.integer(min(start)),
      end   = as.integer(max(end)),
      strand = data.table::first(na.omit(strand))
    ), by = .(seqname, gene_id)]
    # add gene_name: first non-NA per gene_id from anywhere in the GTF
    name_map <- gtf_dt[!is.na(gene_id) & !is.na(gene_name),
                       .(gene_name = data.table::first(gene_name)), by = gene_id]
    genes_dt <- merge(genes_dt, name_map, by = "gene_id", all.x = TRUE)
  }

  # Fallback: if gene_name missing, use gene_id
  genes_dt[is.na(gene_name) | gene_name == "", gene_name := gene_id]

  # Normalize seqnames to chr*
  genes_dt[, seqname := normalize_seqlevels_chr(seqname)]

  # Build GRanges
  GRanges(
    seqnames = genes_dt$seqname,
    ranges   = IRanges(start = genes_dt$start, end = genes_dt$end),
    strand   = genes_dt$strand,
    gene_id  = genes_dt$gene_id,
    gene_name = genes_dt$gene_name
  )
}

# --------------------------------
# Main
# --------------------------------
# 1) Parse SNVs
snv_dt <- parse_snv_for_scatter(snv_dir)

# Save parsed SNV data
save(snv_dt, file = "snv_dt.RData")
fwrite(snv_dt, "snv_dt.tsv", sep = "\t")

# Order chromosomes (factor) for plotting
chr_levels <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
snv_dt[, chr := factor(chr, levels = chr_levels)]

# 2) Build TxDb (not strictly needed for overlap but kept for completeness)
txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")  # from txdbmaker package

# 3) Create gene GRanges (robust to exon-only GTFs)
genes_named_gr <- gtf_to_gene_granges(gtf_file)

# 4) Annotate SNVs with gene_name (overlap)
snv_gr <- GRanges(
  seqnames = as.character(snv_dt$chr),
  ranges   = IRanges(start = snv_dt$pos, end = snv_dt$pos)
)

# Keep only seqlevels in common to avoid warnings
common_seqs <- intersect(seqlevels(snv_gr), seqlevels(genes_named_gr))
snv_gr <- keepSeqlevels(snv_gr, common_seqs, pruning.mode = "coarse")
genes_named_gr <- keepSeqlevels(genes_named_gr, common_seqs, pruning.mode = "coarse")

hits <- findOverlaps(snv_gr, genes_named_gr, ignore.strand = TRUE)

snv_dt[, gene_name := NA_character_]
if (length(hits) > 0) {
  snv_dt$gene_name[queryHits(hits)] <- mcols(genes_named_gr)$gene_name[subjectHits(hits)]
}

# 5) Filter for MM driver genes
mm_driver_genes <- c("MYC","KRAS","NRAS","CDKN2C","FAM46C","TRAF3","RB1",
                     "SP140","EMSY","MFF","PRR14L","PRSS2","SDCCAG8","TP53")

snv_mm <- snv_dt[!is.na(gene_name) & gene_name %in% mm_driver_genes]

# Handle case with no hits to avoid empty-plot errors
if (nrow(snv_mm) == 0) {
  message("Note: No SNVs overlapped MM driver genes in this dataset.")
  # create a dummy row to avoid ggplot errors, then immediately filter legends off
  # but we'll instead skip plot 1 gracefully
  make_plot1 <- FALSE
} else {
  make_plot1 <- TRUE
}

# 6) Plots
samples <- unique(snv_dt$sample)
n       <- length(samples)
cols    <- viridis(n, option="D")
shapes  <- rep(c(15, 17, 8), length.out = n)

# Plot 1: MM Driver Gene Scatter Plot (only if we have data)
if (make_plot1) {
  p1 <- ggplot(snv_mm, aes(x = pos, y = gene_name, color = sample, shape = sample)) +
    geom_point(alpha = 0.8, size = 3) +
    scale_color_manual(values = cols) +
    scale_shape_manual(values = shapes) +
    theme_bw() +
    labs(title = "SNV Scatter Plot - MM Driver Genes",
         x = "Genomic Position", y = "Gene") +
    theme(axis.text.y = element_text(size=10),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.title = element_text(size=10),
          legend.text  = element_text(size=8))
  ggsave("plot_mm_driver_genes.png", p1, width=12, height=6, dpi=150)
  ggsave("plot_mm_driver_genes.pdf", p1, width=12, height=6)
}

# Plot 2: Genome-wide SNV Density (faceted by sample)
p2 <- ggplot(snv_dt, aes(x = pos, y = chr)) +
  stat_bin2d(bins = 200) +
  scale_fill_viridis_c() +
  theme_bw() +
  labs(title = "Genome-wide SNV density per sample",
       x = "Genomic Position", y = "Chromosome") +
  facet_wrap(~sample, ncol = 2)
ggsave("plot_snv_density_per_sample.png", p2, width=12, height=8, dpi=150)
ggsave("plot_snv_density_per_sample.pdf", p2, width=12, height=8)

# Plot 3: Genome-wide SNV Scatter by Chromosome
p3 <- ggplot(snv_dt, aes(x = pos, y = chr, color = sample, shape = sample)) +
  geom_point(alpha = 0.5, size = 1.5) +
  scale_color_manual(values = cols) +
  scale_shape_manual(values = shapes) +
  theme_bw() +
  labs(title="Genome-wide SNV Scatter - by Chromosome",
       x="Genomic Position", y="Chromosome") +
  theme(axis.text.y = element_text(size=10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size=10),
        legend.text = element_text(size=8))
ggsave("plot_snv_scatter_by_chr.png", p3, width=12, height=8, dpi=150)
ggsave("plot_snv_scatter_by_chr.pdf", p3, width=12, height=8)

message("✅ All plots saved as PNG and PDF. Annotation complete.")
