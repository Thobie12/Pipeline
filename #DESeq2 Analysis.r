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

# Load libraries
setwd("/cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/R") # Set working directory here
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(apeglm)
library(EnhancedVolcano)
library(DESeq2)
library(clusterProfiler)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(enrichplot)
library(ReactomePA)
library(DOSE)

counts_FR4all <- read.csv("/cluster/home/t922316uhn/R/RSEMfr4.tsv", sep = "\t", row.names = 1, header = TRUE, stringsAsFactors = TRUE) #link for RSEM expected count data
coldata_FR4all <- read.csv("metadata2.tsv", sep = "\t", header = TRUE, row.names = 1 ) #link for metadata if it exists; make as a tsv file using nanos and tab to separate entry
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
plotPCA(vsd, intgroup = "condition")  # replace with your group


# Get top 20 significant genes and plot heatmap
top_genes <- head(order(res$padj), 20)
mat <- assay(vsd)[top_genes, ]
mat <- mat - rowMeans(mat)
pheatmap(mat, annotation_col = coldata_FR4all)

#plot dispersion
plotDispEsts(dds)

# Begin saving to PDF
pdf("Plots.pdf", width = 10, height = 8)
plotMA(res, ylim = c(-5, 5), main = "MA Plot")
EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'padj',
    pCutoff = 0.05,
    FCcutoff = 1,
    pointSize = 2.0,
    labSize = 2.5,
    title = 'Volcano Plot')
plotPCA(vsd, intgroup = "Condition")
plotDispEsts(dds, main = "Dispersion Estimates")
top_genes <- head(order(res$padj), 20)
mat <- assay(vsd)[top_genes, ]
mat <- mat - rowMeans(mat)
pheatmap(mat, annotation_col = as.data.frame(colData(dds)[, "Condition", drop = FALSE]),
         main = "Top 20 Differentially Expressed Genes")

dev.off()
# Begin saving to PDF
pdf("Plots.pdf", width = 10, height = 8)

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
dev.off()

## 6. GSEA Analysis
sig_genes <- rownames(res[which(res$padj < 0.05), ])
entrez_ids <- mapIds(org.Hs.eg.db, keys = sig_genes, keytype = "ENSEMBL", column = "ENTREZID")
# Step 1: Remove version numbers (e.g., ".1" at the end)
sig_genes_clean <- gsub("\\..*$", "", sig_genes)
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
# Step 3: KEGG
kegg_res <- enrichKEGG(gene         = entrez_ids,
                       organism     = "hsa",
                       pvalueCutoff = 0.05)
# Step 4: Reactome
reactome_res <- enrichPathway(gene = entrez_ids,
                              organism = "human",
                              pvalueCutoff = 0.05,
                              readable = TRUE)


BiocManager::install(c("ReactomePA", "DOSE"))
reactome_res <- enrichPathway(gene = entrez_ids,
                              organism = "human",
                              pvalueCutoff = 0.05,
                              readable = TRUE)
# Step 5: GSEA
gene_list <- res$log2FoldChange
names(gene_list) <- gsub("\\..*$", "", rownames(res))  # clean names
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

barplot(go_res, showCategory = 20, title = "GO Biological Process")
dotplot(go_res, showCategory = 20, title = "GO BP (dotplot)")

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
dev.off()
# STEP 4: Export results as CSV
write.csv(as.data.frame(go_res), "GO_results.csv")
write.csv(as.data.frame(kegg_res), "KEGG_results.csv")
write.csv(as.data.frame(reactome_res), "Reactome_results.csv")
write.csv(as.data.frame(gsea_res), "GSEA_results.csv")

#savehistory("DEgenecode.R") - code to save history
