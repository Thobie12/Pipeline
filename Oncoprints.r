#Oncoprint plots for CNA + SNV + SV
# --- Libraries ---
library(data.table)
library(dplyr)
library(tidyr)
library(vcfR)
library(VariantAnnotation)
library(GenomicRanges)
library(GenomicFeatures)
library(ComplexHeatmap)
library(circlize)
library(IRanges)
library(S4Vectors)
#library(org.Hs.eg.db)
#library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# --- User Inputs ---
#snv_dir <- "~/Desktop/UofT/PhD/Pugh Lab/Projects/MM_Patient_Samples/UltimaWGS/GATK"
snv_dir <- "/cluster/home/t922316uhn/output/GATK"
#sv_dir  <- "~/Desktop/UofT/PhD/Pugh Lab/Projects/MM_Patient_Samples/UltimaWGS/Delly"
sv_dir <- "/cluster/home/t922316uhn/output/Delly"
#sv_dir2 <- "~/Desktop/UofT/PhD/Pugh Lab/Projects/MM_Patient_Samples/UltimaWGS/GRIDSS"
sv_dir2 <- "/cluster/home/t922316uhn/output/GRIDSS"

#---- create ichorCNA segment file ----
# run this only once to merge segment files
#seg_dir <- "/cluster/home/t922316uhn/output/Test_R"
#seg_files <- list.files(seg_dir, pattern = "\\.seg$", full.names = TRUE)

#all_seg <- rbindlist(lapply(seg_files, fread), fill = TRUE)

#Standardize copy number column name for parse_ichor()
#setnames(all_seg, old = "copy.number", new = "copyNumber")

# Save merged segment file
#merged_file <- file.path(seg_dir, "merged_segments.tsv")
#fwrite(all_seg, merged_file, sep = "\t")

#message("Merged segment file saved to: ", merged_file)
#ichor_segment_file <- "~/Desktop/UofT/PhD/Pugh Lab/Projects/MM_Patient_Samples/UltimaWGS/ichorCNA/Test_R/merged_segments.tsv"
ichor_segment_file <- "/cluster/home/t922316uhn/output/Test_R/merged_segments.tsv"
ichor_gene_file <- NULL
output_prefix <- "oncoprint_result"
genome_build <- "hg38"
cn_amp_threshold <- 2.5
cn_del_threshold <- 1.5

# ===========================
# --- Helper Functions ---
# ===========================

# --- Load gene ranges for MAC---
#get_gene_ranges <- function(){
#  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#  genes_gr <- genes(txdb)
#  mcols(genes_gr)$gene_id <- as.character(mcols(genes_gr)$gene_id)
#  return(genes_gr)
# }

# --- Load gene ranges from GTF file for cluster usage ---
library(GenomicRanges)
library(rtracklayer)
library(dplyr)

gtf_path <- "/cluster/tools/data/genomes/human/hg38/iGenomes/Annotation/Genes/genes.gtf"

get_gene_ranges <- function(gtf_file = gtf_path) {
  # Import GTF
  gtf <- import(gtf_file)
  
  # Ensure type is character
  gtf$type <- as.character(gtf$type)
  
  # Keep only exons (CDS could be used too, but exons are safer)
  exons <- gtf[gtf$type == "exon"]
  
  # Ensure gene_id exists
  if (!"gene_id" %in% colnames(mcols(exons))) {
    stop("gene_id column not found in GTF!")
  }
  
  # Aggregate exons to get gene ranges
  exons_df <- as.data.frame(exons)
  library(dplyr)
  gene_ranges_df <- exons_df %>%
    group_by(seqnames, strand, gene_id, gene_name) %>%
    summarize(start = min(start),
              end = max(end),
              .groups = "drop")
  
  # Convert back to GRanges
  genes_gr <- GRanges(
    seqnames = gene_ranges_df$seqnames,
    ranges = IRanges(start = gene_ranges_df$start, end = gene_ranges_df$end),
    strand = gene_ranges_df$strand,
    gene_id = gene_ranges_df$gene_id,
    gene_name = gene_ranges_df$gene_name
  )
  
  return(genes_gr)
}

# Usage
gene_ranges <- get_gene_ranges()
gene_ranges

# --- Parse SNV VCFs ---
parse_snv_dir <- function(dir, genes_gr){
  files <- list.files(dir, pattern="\\.vcf(\\.gz)?$", full.names=TRUE)
  if(length(files) == 0) return(data.table())

  out <- lapply(files, function(f){
    sample <- sub("\\.vcf(\\.gz)?$", "", basename(f))
    message("Parsing SNV VCF: ", sample)
    vcf <- tryCatch(readVcf(f, genome_build), error=function(e) return(NULL))
    if(is.null(vcf) || length(vcf) == 0) return(NULL)
    
    vr <- rowRanges(vcf)
    gr <- GRanges(seqnames=seqnames(vr), ranges=ranges(vr), strand="*")
    hits <- findOverlaps(gr, genes_gr)
    if(length(hits) == 0) return(NULL)
    
    df <- data.table(
      sample = sample,
      gene = as.character(mcols(genes_gr)$gene_id[subjectHits(hits)])
    )
    unique(df[, .(sample, gene)])[ , event := "SNV"]
  })
  rbindlist(out, use.names=TRUE, fill=TRUE)
}

# --- Parse SV VCFs ---
parse_sv_dir <- function(dir, genes_gr){
  files <- list.files(dir, pattern="\\.vcf(\\.gz)?$", full.names=TRUE)
  if(length(files) == 0) return(data.table())

  out <- lapply(files, function(f){
    sample <- sub("\\.vcf(\\.gz)?$", "", basename(f))
    message("Parsing SV VCF: ", sample)
    vcf <- tryCatch(readVcf(f, genome_build), error=function(e) return(NULL))
    if(is.null(vcf) || length(vcf) == 0) return(NULL)
    
    vr <- rowRanges(vcf)
    info_df <- info(vcf)
    sv_end <- ifelse(!is.na(info_df$END), as.integer(info_df$END), end(ranges(vr)))
    gr <- GRanges(seqnames=seqnames(vr), ranges=IRanges(start=start(ranges(vr)), end=sv_end), strand="*")
    svtypes <- as.character(info_df$SVTYPE)
    
    hits <- findOverlaps(gr, genes_gr)
    if(length(hits) == 0) return(NULL)
    
    df <- data.table(
      sample = sample,
      gene = as.character(mcols(genes_gr)$gene_id[subjectHits(hits)]),
      svtype = svtypes[queryHits(hits)]
    )
    df <- unique(df[, .(sample, gene, svtype)])
    df[, event := paste0("SV:", svtype)]
  })
  rbindlist(out, use.names=TRUE, fill=TRUE)
}

# --- Parse ichorCNA ---
parse_ichor <- function(segment_file=NULL, gene_file=NULL, genes_gr=NULL){
  if(!is.null(gene_file) && file.exists(gene_file)){
    dt <- fread(gene_file)
    cn_col <- intersect(c("cn","total_cn","copyNumber","copy_number","CN"), names(dt))[1]
    if(is.na(cn_col)) stop("Cannot find CN column.")
    dt[, event := ifelse(get(cn_col) >= cn_amp_threshold, "CNA:AMP",
                         ifelse(get(cn_col) <= cn_del_threshold, "CNA:DEL", NA))]
    return(unique(dt[!is.na(event), .(sample,gene,event)]))
  }
  else if(!is.null(segment_file) && file.exists(segment_file) && !is.null(genes_gr)){
    seg <- fread(segment_file)
    chr_col <- intersect(c("chrom","chr","Chromosome"), names(seg))[1]
    start_col <- intersect(c("start","seg_start","startpos"), names(seg))[1]
    end_col <- intersect(c("end","seg_end","endpos"), names(seg))[1]
    cn_col <- intersect(c("copyNumber","cn","segmean","tcn"), names(seg))[1]
    if(any(is.na(c(chr_col,start_col,end_col,cn_col)))) stop("Segment file missing columns")
    if(!"sample" %in% names(seg)) stop("Segment file must include 'sample'")
    
    seg_gr <- GRanges(seqnames = seg[[chr_col]], ranges=IRanges(seg[[start_col]], seg[[end_col]]))
    hits <- findOverlaps(genes_gr, seg_gr)
    if(length(hits) == 0) return(data.table())
    
    mapped <- data.table(
      gene = as.character(mcols(genes_gr)$gene_id[queryHits(hits)]),
      sample = seg$sample[subjectHits(hits)],
      cn = as.numeric(seg[[cn_col]][subjectHits(hits)])
    )
    mapped[, event := ifelse(cn >= cn_amp_threshold, "CNA:AMP",
                             ifelse(cn <= cn_del_threshold, "CNA:DEL", NA))]
    unique(mapped[!is.na(event), .(sample,gene,event)])
  } else {
    return(data.table())
  }
}

# --- Build alteration matrix ---
build_alteration_matrix <- function(...){
  combined <- rbindlist(list(...), use.names=TRUE, fill=TRUE)
  if(nrow(combined) == 0) return(matrix(character(0),0,0))
  mat <- dcast(combined, gene ~ sample, value.var="event", fun.aggregate=function(x) paste(unique(x), collapse=";"))
  mat <- as.matrix(mat[,-1, drop=FALSE])
  rownames(mat) <- dcast(combined, gene ~ sample, value.var="event")[,1]
  return(mat)
}

# ===========================
# --- Run Parsing ---
# ===========================

genes_gr <- get_gene_ranges()
snv_dt <- parse_snv_dir(snv_dir, genes_gr)
sv_dt  <- parse_sv_dir(sv_dir, genes_gr)
cna_dt <- parse_ichor(ichor_segment_file, ichor_gene_file, genes_gr)

# Save combined events
combined_dt <- rbindlist(list(
 snv_dt[, .(sample,gene,event)],
 sv_dt[, .(sample,gene,event)],
 cna_dt[, .(sample,gene,event)]
), use.names=TRUE, fill=TRUE)
if(nrow(combined_dt) > 0){
  fwrite(combined_dt, paste0(output_prefix,"_combined_events.tsv"), sep="\t")
  message("Saved combined event table.")
}

#--------------------------------------------------------------------------------------------------
library(data.table)
library(ComplexHeatmap)
library(circlize)
#library(org.Hs.eg.db)

output_prefix <- "oncoprint_result"

# --- Skip parsing, just load combined events --- do this if using parsed file
#combined_dt <- fread("~/Desktop/oncoprint_result_combined_events.tsv")

# --- Build alteration matrix ---
build_alteration_matrix <- function(combined_dt) {
  combined_unique <- combined_dt[, .(event = paste(unique(event), collapse = ";")),
                                 by = .(gene, sample)]
  mat <- dcast.data.table(combined_unique, gene ~ sample, value.var = "event", fill = "")
  mat <- as.data.frame(mat)
  rownames(mat) <- mat$gene
  mat$gene <- NULL
  return(mat)
}

mat <- build_alteration_matrix(combined_dt)

# --- Convert to onco_mat list ---
onco_mat <- apply(mat, c(1,2), function(cell){
  if(cell=="") return(NULL)
  tokens <- unique(strsplit(cell, ";")[[1]])
  lapply(tokens, function(t){
    if(grepl("^SV:", t)) "SV" else if(t=="SNV") "SNV" else t
  })
})

# Keep original gene names as rownames
rownames(onco_mat) <- rownames(onco_mat)

# --- Filter for MM driver genes ---
mm_driver_genes <- c("MYC","KRAS","NRAS","CDKN2C","FAM46C","TRAF3","RB1",
                     "SP140","EMSY","MFF","PRR14L","PRSS2","SDCCAG8","TP53")
onco_mat <- onco_mat[rownames(onco_mat) %in% mm_driver_genes, , drop=FALSE]

# --- Top Annotation: Tumor Fraction ---
# Read tumor fraction TSV
tumor_file <- "/cluster/home/t922316uhn/output/Test_R/tumorfraction.tsv"
tumor_dt <- fread(tumor_file)

# Ensure sample names match the column names of onco_mat
tumor_fraction_vec <- setNames(tumor_dt$tumor_fraction[match(colnames(onco_mat), tumor_dt$sample)],
                               colnames(onco_mat))

# Replace NAs with 0 (optional)
tumor_fraction_vec[is.na(tumor_fraction_vec)] <- 0

# Create the top annotation
top_anno <- HeatmapAnnotation(
  TumorFraction = anno_barplot(tumor_fraction_vec, gp=gpar(fill="#1f78b4"), border=FALSE)
)

# --- Right Annotation: % Altered ---
pct_vec <- apply(onco_mat, 1, function(x) sum(!sapply(x,is.null))/ncol(onco_mat)*100)
right_anno <- rowAnnotation(`%Altered` = anno_barplot(pct_vec, gp=gpar(fill="#FF69B4"), border=FALSE, ylim=c(0,100)))


#-------------------------------------------------------------------------------------------------
library(ComplexHeatmap)
library(grid)

# --- 1. Define all possible alteration types ---
alter_types <- c("SV", "CNA", "AMP", "DEL", "INS", "INV", "DUP", "BND", "SNV")

# --- 2. Convert your onco_mat into boolean matrix per type ---
# Ensure your matrix is character
onco_mat <- apply(as.matrix(onco_mat), c(1,2), as.character)

# Create an empty list to store boolean matrices
mat_list <- list()
for (atype in alter_types) {
  mat_list[[atype]] <- matrix(
    sapply(onco_mat, function(x) grepl(atype, x)),
    nrow = nrow(onco_mat),
    ncol = ncol(onco_mat),
    dimnames = dimnames(onco_mat)
  )
}

# --- 3. Define alter_fun for each type ---
alter_fun <- list(
  background = function(x, y, w, h) grid.rect(x, y, w, h, gp=gpar(fill="#CCCCCC", col=NA)),
  SV = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp=gpar(fill="#FFD700", col=NA)),
  CNA = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp=gpar(fill="#1E90FF", col=NA)),
  AMP = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp=gpar(fill="#32CD32", col=NA)),
  DEL = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp=gpar(fill="#FF4500", col=NA)),
  INS = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp=gpar(fill="#8A2BE2", col=NA)),
  INV = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp=gpar(fill="#FF69B4", col=NA)),
  DUP = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp=gpar(fill="#00CED1", col=NA)),
  BND = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp=gpar(fill="#FFA500", col=NA)),
  SNV = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp=gpar(fill="#A52A2A", col=NA))
)

# --- 4. Combine boolean matrices into a single "oncoPrint-ready" matrix ---
onco_print_mat <- matrix(nrow=nrow(onco_mat), ncol=ncol(onco_mat))
dimnames(onco_print_mat) <- dimnames(onco_mat)

for (i in seq_len(nrow(onco_mat))) {
  for (j in seq_len(ncol(onco_mat))) {
    types_present <- alter_types[sapply(mat_list, function(m) m[i,j])]
    if (length(types_present) == 0) {
      onco_print_mat[i,j] <- NA
    } else if (length(types_present) == 1) {
      onco_print_mat[i,j] <- types_present
    } else {
      onco_print_mat[i,j] <- paste(types_present, collapse=",")
    }
  }
}

# --- 5. Define colors ---
col <- c(
  SV = "#FFD700",
  CNA = "#1E90FF",
  AMP = "#32CD32",
  DEL = "#FF4500",
  INS = "#8A2BE2",
  INV = "#FF69B4",
  DUP = "#00CED1",
  BND = "#FFA500",
  SNV = "#A52A2A"
)

# --- 6. Draw oncoprint ---
ht <- oncoPrint(
  onco_print_mat,
  alter_fun = alter_fun,
  alter_fun_is_vectorized = TRUE,
  col = col,
  remove_empty_rows = TRUE,
  remove_empty_columns = FALSE,
  show_column_names = TRUE,
  show_row_names = TRUE,
  row_names_side = "left",
  pct_side = "right",
  column_title = "Oncoprint - SNV + CNA + SV",
  heatmap_legend_param = list(title="Alterations", at=names(col), labels=names(col))
)


# --- Save PDF & PNG ---
pdf(paste0(output_prefix, "_oncoprint.pdf"), width=12, height=8)
draw(ht, heatmap_legend_side = "right")
dev.off()

png(paste0(output_prefix, "_oncoprint.png"), width=1600, height=1200, res=150)
draw(ht, heatmap_legend_side = "right")
dev.off()
