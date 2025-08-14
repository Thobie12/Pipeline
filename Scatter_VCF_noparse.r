#Scatter plot for VCF with no parsing
#!/usr/bin/env Rscript

# ===============================
# SNV Annotation & Visualization (from saved table)
# ===============================

library(data.table)
library(GenomicRanges)
library(txdbmaker)   # replaces makeTxDbFromGFF
library(ggplot2)
library(viridis)

# --- Input ---
snv_tsv  <- "snv_dt.tsv"  # pre-parsed SNV table
gtf_file <- "/cluster/tools/data/genomes/human/hg38/iGenomes/Annotation/Genes/genes.gtf"

# ============================================================
# Load SNV table
# ============================================================
snv_dt <- fread(snv_tsv)
chr_levels <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
snv_dt[, chr := factor(chr, levels = chr_levels)]

# ============================================================
# Load gene annotations from GTF (handle exon-only GTF)
# ============================================================
#gtf_dt <- fread(gtf_file, header = FALSE, sep = "\t", skip = "#")
gtf_dt <- fread(gtf_file, header = FALSE, sep = "\t")
gtf_dt <- gtf_dt[!grepl("^#", V1)]
setnames(gtf_dt, c("seqname","source","feature","start","end","score","strand","frame","attribute"))

# Extract gene_id and gene_name
extract_attr <- function(attr_string, key){
  sub(paste0('.*', key, ' "([^"]+)".*'), "\\1", attr_string)
}
gtf_dt[, gene_id   := extract_attr(attribute, "gene_id")]
gtf_dt[, gene_name := extract_attr(attribute, "gene_name")]

# Collapse to gene-level coordinates if only exon entries
genes_named <- gtf_dt[feature == "gene"]
if (nrow(genes_named) == 0) {
  message("No 'gene' entries found in GTF — collapsing exons to gene-level coords.")
  genes_named <- gtf_dt[feature == "exon", .(
    start = min(start),
    end   = max(end),
    gene_name = unique(na.omit(gene_name))[1]
  ), by = .(seqname, strand, gene_id)]
}

# Make GRanges
genes_named_gr <- GRanges(
  seqnames = genes_named$seqname,
  ranges   = IRanges(start = genes_named$start, end = genes_named$end),
  gene_name = genes_named$gene_name
)

# ============================================================
# Annotate SNVs with gene names
# ============================================================
snv_gr <- GRanges(
  seqnames = snv_dt$chr,
  ranges   = IRanges(start = snv_dt$pos, end = snv_dt$pos)
)

hits <- findOverlaps(snv_gr, genes_named_gr)

snv_dt$gene_name <- NA_character_
snv_dt$gene_name[queryHits(hits)] <- mcols(genes_named_gr)$gene_name[subjectHits(hits)]

# ============================================================
# Filter for MM driver genes
# ============================================================
mm_driver_genes <- c("MYC","KRAS","NRAS","CDKN2C","FAM46C","TRAF3","RB1",
                     "SP140","EMSY","MFF","PRR14L","PRSS2","SDCCAG8","TP53")

snv_mm <- snv_dt[gene_name %in% mm_driver_genes]

# ============================================================
# Plot 1: MM Driver Gene Scatter Plot
# ============================================================
samples <- unique(snv_mm$sample)
n       <- length(samples)
cols    <- viridis(n, option="D")
shapes  <- rep(c(15, 17, 8), length.out = n)

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

# ============================================================
# Plot 2: Genome-wide SNV Density (faceted by sample)
# ============================================================
p2 <- ggplot(snv_dt, aes(x = pos, y = chr)) +
  stat_bin2d(bins = 200) +
  scale_fill_viridis_c() +
  theme_bw() +
  labs(title = "Genome-wide SNV density per sample") +
  facet_wrap(~sample, ncol = 2)

ggsave("plot_snv_density_per_sample.png", p2, width=12, height=8, dpi=150)
ggsave("plot_snv_density_per_sample.pdf", p2, width=12, height=8)

# ============================================================
# Plot 3: Genome-wide SNV Scatter by Chromosome
# ============================================================
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
