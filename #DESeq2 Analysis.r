#DESeq2 Analysis.R
#includes plotMA, Volcano Plot, PCA, Heatmap, and GSEA.

# install libraries if not already installed
#BiocManager::install("apeglm")
#BiocManager::install("EnhancedVolcano")
#BiocManager::install("DESeq2")
#if (!requireNamespace("clusterProfiler"))
#install.packages("BiocManager"); BiocManager::install("clusterProfiler")
#BiocManager::install(c("org.Hs.eg.db", "enrichplot"))
#BiocManager::install(c("ReactomePA", "DOSE"))
#install.packages("circlize")

# Load libraries
setwd("/cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/KMS26/RSEM") # Set working directory here
library(ggplot2)
library(pheatmap)
library(apeglm)
library(EnhancedVolcano)
library(DESeq2)
library(clusterProfiler)
library(enrichplot)
library(circlize)
library(dplyr)
library(readr)
library(stringr)
library(org.Hs.eg.db)
library(ReactomePA)
library(DOSE)

counts_FR4all <- read.csv("2025-04-21_FR4_KMS26_Apr2025_rsem_expected_counts.tsv", sep = "\t", row.names = 1, header = TRUE, stringsAsFactors = TRUE) #link for RSEM expected count data
coldata_FR4all <- read.csv("metadata.tsv", sep = "\t", header = TRUE, row.names = 1 ) #link for metadata if it exists; make as a tsv file using nanos and tab to separate entry
nrow(coldata_FR4all) # Check number of rows in coldata
colnames(coldata_FR4all) # Check column names in coldata
geneinfo_FR4all <- counts_FR4all[, 1:3] # Extract gene information
CleanedCounts_FR4all <- counts_FR4all[, -(1:2)] # Remove gene information columns
CleanedCounts_FR4all <- round(CleanedCounts_FR4all) # Round counts to integers
all(colnames(CleanedCounts_FR4all) == rownames(coldata_FR4all)) # Check if column names of counts match row names of coldata; if FALSE, check the metadata file and fix
coldata_FR4all$Condition <- as.factor(coldata_FR4all$Condition) # Convert condition to factor
dds <- DESeqDataSetFromMatrix(countData = CleanedCounts_FR4all,
                              colData = coldata_FR4all,
                              design = ~ Condition)
#Clustering before dds

#Clustering independently: without any grouping
# 1. Calculate sample-to-sample distances
vsd <- vst(dds)
sampleDists <- dist(t(assay(vsd)))  # transpose so samples are compared, not genes

# 2. Turn distance into a matrix
sampleDistMatrix <- as.matrix(sampleDists)

# 3. Label rows and columns with sample names
rownames(sampleDistMatrix) <- colnames(vsd)
colnames(sampleDistMatrix) <- colnames(vsd)

pdf("Preplots.pdf") #test plots for clustering
# 4. Cluster and plot heatmap
pheatmap(sampleDistMatrix, 
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Sample-to-sample distance heatmap")
# Plot PCA
#Perform vst transformation
vsd <- vst(dds)

# Perform PCA without any grouping variable (no intgroup argument)
# Perform PCA without grouping
pcaData <- prcomp(t(assay(vsd)))$x
pcaData <- as.data.frame(pcaData)
pcaData$Sample <- rownames(pcaData)  # Add sample names for labeling

# Create PCA plot using ggplot
library(ggplot2)
p <- ggplot(pcaData, aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = rownames(pcaData)), size = 3) +
  labs(title = "PCA of Samples (no grouping)", 
       x = paste0("PC1: ", round(100 * pcaData$percentVar[1]), "% variance"),
       y = paste0("PC2: ", round(100 * pcaData$percentVar[2]), "% variance")) +
  theme_minimal()
print(p)
# Perform PCA with grouping variable
#PCA Clustering with grouping
vsd <- vst(dds)

# PCA plot
plotPCA(vsd, intgroup = "Condition")

# Heatmap clustering (example)
library(pheatmap)
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 500)
pheatmap(assay(vsd)[topVarGenes, ])

dev.off()

#Run DESeq2
dds <- DESeq(dds) # Run DESeq2
#check results
res <- results(dds)
res
summary(res)
res_sig <- res[which(res$padj < 0.05), ] # Get significant genes using p < 0.05
head(res_sig)
summary(res_sig)
res_sig <- res[which(res$padj < 0.05), ]
head(res_sig)
sink("summary.txt") #save to file
summary(res)
summary(res_sig)
sink()

# bonferroni correction if needed/planned
#res$padj_bonferroni <- p.adjust(res$pvalue, method = "bonferroni")
#res_ordered <- res[order(res$pad_bonferroni), ]
#res$pad_bonferroni <- p.adjust(res$pvalue, method = "bonferroni")
#res_ordered <- res[order(res$pad_bonferroni), ]
#sig_genes <- subset(res_ordered, padj_bonferroni < 0.05)
#summary(sig_genes)

#Making Plots and Saving to PDF

#plotMA
plotMA(res, ylim = c(-5, 5))
png("MA_plot.png", width = 800, height = 600)
plotMA(res, ylim = c(-5, 5))

# Volcano Plot
EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'padj',
    pCutoff = 0.05,
    FCcutoff = 1,
    pointSize = 2.0,
    labSize = 3.0)
#PCA Plot
vsd <- vst(dds, blind = FALSE)  # or use rlog(dds) 
plotPCA(vsd, intgroup = "Condition")  # replace with your group


# Get top 20 significant genes and plot heatmap
top_genes <- order(res$padj, na.last = NA)[1:16]
mat <- assay(vsd)[top_genes, ]
mat <- mat - rowMeans(mat)
pheatmap(mat, annotation_col = coldata_FR4all)

#plot dispersion
plotDispEsts(dds)

# Begin saving to PDF
pdf("Plots.pdf")

## 1. MA Plot
plotMA(res, ylim = c(-5, 5), main = "MA Plot")

## 2. Volcano Plot
EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'padj',
    pCutoff = 0.05,
    FCcutoff = 1,
    pointSize = 2.0,
    labSize = 2.5,
    title = 'Volcano Plot')

## 3. PCA Plot
plotPCA(vsd, intgroup = "Condition")  # Change "condition" to your column name

## 4. Dispersion Plot
plotDispEsts(dds, main = "Dispersion Estimates")

## 5. Heatmap of Top Genes
top_genes <- head(order(res$padj), 20)
mat <- assay(vsd)[top_genes, ]
mat <- mat - rowMeans(mat)
pheatmap(mat, annotation_col = as.data.frame(colData(dds)[, "Condition", drop = FALSE]),
         main = "Top 20 Differentially Expressed Genes")


#save plots to PDF

## 6. GSEA Analysis
# Step 1: Remove version numbers (e.g., ".1" at the end)
sig_genes <- rownames(res[which(res$padj < 0.05), ])
sig_genes_clean <- gsub("\\..*$", "", sig_genes)
entrez_ids <- mapIds(org.Hs.eg.db, keys = sig_genes_clean, keytype = "ENSEMBL", column = "ENTREZID")

# Step 2: Convert to Entrez
entrez_ids <- mapIds(org.Hs.eg.db,
                     keys = sig_genes_clean,
                     keytype = "ENSEMBL",
                     column = "ENTREZID",
                     multiVals = "first")
go_results <- enrichGO(gene         = entrez_ids,
                       OrgDb        = org.Hs.eg.db,
                       keyType      = "ENTREZID",
                       ont          = "BP",       # BP = Biological Process, also "MF", "CC"
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       readable      = TRUE)
# Step 3: KEGG # needs internet
kegg_res <- enrichKEGG(gene         = entrez_ids,
                       organism     = "hsa",
                       pvalueCutoff = 0.05)
# Step 4: Reactome
reactome_res <- enrichPathway(gene = entrez_ids,
                              organism = "human",
                              pvalueCutoff = 0.05,
                              readable = TRUE)
# Step 5: GSEA
gene_list <- res$log2FoldChange
names(gene_list) <- ifelse(grepl("\\..*$", rownames(res)), 
                           gsub("\\..*$", "", rownames(res)), 
                           rownames(res))  # clean names only if version suffix exists
gene_list <- mapIds(org.Hs.eg.db,
                    keys = names(gene_list),
                    keytype = "ENSEMBL",
                    column = "ENTREZID",
                    multiVals = "first")
gene_list <- na.omit(gene_list)
names(res$log2FoldChange) <- gene_list
gene_list <- sort(res$log2FoldChange, decreasing = TRUE)
gsea_res <- gseGO(geneList = gene_list,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  keyType = "ENTREZID",
                  minGSSize = 10,
                  maxGSSize = 500,
                  pvalueCutoff = 0.05,
                  verbose = FALSE)

go_res <- enrichGO(gene         = entrez_ids,
                   OrgDb        = org.Hs.eg.db,
                   keyType      = "ENTREZID",
                   ont          = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   readable      = TRUE)


barplot(go_res, showCategory = 20, title = "GO Biological Process")
dotplot(go_res, showCategory = 20, title = "GO BP (dotplot)")
# KEGG
barplot(kegg_res, showCategory = 20, title = "KEGG Pathways")
dotplot(kegg_res, showCategory = 20, title = "KEGG (dotplot)")

# Reactome
barplot(reactome_res, showCategory = 20, title = "Reactome Pathways")
dotplot(reactome_res, showCategory = 20, title = "Reactome (dotplot)")

# GSEA
dotplot(gsea_res, showCategory = 20, title = "GSEA (GO BP)")
ridgeplot(gsea_res, showCategory = 20)

# STEP 4: Export results as CSV
write.csv(as.data.frame(go_res), "GO_results.csv")
write.csv(as.data.frame(kegg_res), "KEGG_results.csv")
write.csv(as.data.frame(reactome_res), "Reactome_results.csv")
write.csv(as.data.frame(gsea_res), "GSEA_results.csv")

## Fusion Analysis
fusion_data <- read_tsv("2025-04-17_FR4TP53_fusioncatcher_output_long.tsv") # path to fusioncatcher output long format
fusion_clean <- fusion_data %>%
  filter(!is.na(`Gene_1_symbol`) & !is.na(`Gene_2_symbol`)) %>%
  select(Gene_1 = `Gene_1_symbol`,
         Chr1   = `Gene_1_chromosome`,
         Pos1   = `Gene_1_position`,
         Gene_2 = `Gene_2_symbol`,
         Chr2   = `Gene_2_chromosome`,
         Pos2   = `Gene_2_position`)
circos_data <- fusion_clean %>%
  mutate(
    Pos1_start = Pos1 - 1000,
    Pos1_end   = Pos1 + 1000,
    Pos2_start = Pos2 - 1000,
    Pos2_end   = Pos2 + 1000
  ) %>%
  select(Chr1, Pos1_start, Pos1_end, Chr2, Pos2_start, Pos2_end)

colnames(fusion_data)

#input
nrow(circos_data)
head(fusion_data$`Fusion_point_for_gene_1.5end_fusion_partner.`, 10)
fusion_clean <- fusion_data %>%
  filter(!is.na(`Fusion_point_for_gene_1.5end_fusion_partner.`),
         !is.na(`Fusion_point_for_gene_2.3end_fusion_partner.`)) %>%
  mutate(
    Chr1 = paste0("chr", str_extract(`Fusion_point_for_gene_1.5end_fusion_partner.`, "^[^:]+")),
    Pos1 = as.numeric(str_extract(`Fusion_point_for_gene_1.5end_fusion_partner.`, "(?<=:)[0-9]+")),
    Chr2 = paste0("chr", str_extract(`Fusion_point_for_gene_2.3end_fusion_partner.`, "^[^:]+")),
    Pos2 = as.numeric(str_extract(`Fusion_point_for_gene_2.3end_fusion_partner.`, "(?<=:)[0-9]+")),
    Gene_1 = `Gene_1_symbol.5end_fusion_partner.`,
    Gene_2 = `Gene_2_symbol.3end_fusion_partner.`
  ) %>%
  filter(!is.na(Chr1), !is.na(Chr2), !is.na(Pos1), !is.na(Pos2))

circos_data <- fusion_clean %>%
  mutate(
    Pos1_start = Pos1 - 1000,
    Pos1_end   = Pos1 + 1000,
    Pos2_start = Pos2 - 1000,
    Pos2_end   = Pos2 + 1000
  ) %>%
  select(Chr1, Pos1_start, Pos1_end, Chr2, Pos2_start, Pos2_end)

cat("Number of fusion links: ", nrow(circos_data), "\n")

circos_data <- fusion_clean %>%
  mutate(
    Pos1_start = Pos1 - 1000,
    Pos1_end = Pos1 + 1000,
    Pos2_start = Pos2 - 1000,
    Pos2_end = Pos2 + 1000
  ) %>%
  select(Chr1, Pos1_start, Pos1_end, Chr2, Pos2_start, Pos2_end)

circos.clear()
circos.par("track.height" = 0.05)
circos.initializeWithIdeogram(species = "hg19")

circos.genomicLink(
  region1 = circos_data[, 1:3],
  region2 = circos_data[, 4:6],
  col = rand_color(nrow(circos_data), transparency = 0.4)
)
# Save the circos plot
circos.clear()
circos.par("track.height" = 0.05)
circos.initializeWithIdeogram(species = "hg19")

circos.genomicLink(
  region1 = circos_data[, 1:3],
  region2 = circos_data[, 4:6],
  col = rand_color(nrow(circos_data), transparency = 0.4)
)

dev.off()
#savehistory("DEgenecode.R") - code to save history
