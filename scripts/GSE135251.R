rm(list=ls())

library(readr)
library(tidyverse)
library(GEOquery)
library(DESeq2)
library(org.Hs.eg.db)
library(dplyr)
library(clusterProfiler)



# List all count files
files <- list.files("C:/Users/chukw/OneDrive/Desktop/capstone project/Rstudio capstone/capstone/counts/", pattern="\\.counts.txt$", full.names=TRUE)

# Read and merge all count files
count_list <- lapply(files, function(f) {
  read.delim(f, header=FALSE, col.names=c("Gene", basename(f)))
})

# Reduce to single data frame by merging on Gene column
count_matrix_2 <- purrr::reduce(count_list, dplyr::full_join, by="Gene")



# Save gene IDs separately as rownames instead of a column
rownames(count_matrix_2) <- count_matrix_2$Gene
count_matrix_2 <- count_matrix_2[,-1]

# Convert dataframe which can contain different types of data to a matrix that can only contain one type
count_matrix_2 <- as.matrix(count_matrix_2)

# Convert to integer matrix
mode(count_matrix_2) <- "integer"

# Quick check
head(count_matrix_2[,1:5])

# stripping uneccessary text
colnames(count_matrix_2) <- gsub("\\.counts\\.txt$", "", colnames(count_matrix_2))
colnames(count_matrix_2) <- gsub("(_.*)", "", colnames(count_matrix_2))  # optional: keep only GSM ID

gse2 <- getGEO("GSE135251", GSEMatrix = TRUE)
pheno2 <- pData(phenoData(gse2[[1]]))
head(pheno2)

# show the columns in the metadata
colnames(pheno2)

# check if the sample IDs in the column names fo the count matrix matches the GEO accession values in the pheno metadata
all(colnames(count_matrix_2) %in% pheno2$geo_accession)


# subset data
pheno2 <- pheno2[match(colnames(count_matrix_2), pheno2$geo_accession), ]

# change the disease:ch1 column to disease to remove the ":" because it interferes wiht parsing in DESeq2 code ran below
colnames(pheno2)[colnames(pheno2) == "disease:ch1"] <- "disease"
colnames(pheno2)[colnames(pheno2) == "fibrosis stage:ch1"] <- "fibrosis_stage"
colnames(pheno2)[colnames(pheno2) == "group in paper:ch1"] <- "group_in_paper"
colnames(pheno2)[colnames(pheno2) == "nas score:ch1"] <- "nas_score"
colnames(pheno2)[colnames(pheno2) == "Stage:ch1"] <- "Stage"

# get the different types of names int he disease.ch1 column of the pheno metadata
unique(pheno2$disease)
unique(pheno2$fibrosis_stage)
unique(pheno2$group_in_paper)
unique(pheno2$nas_score)
unique(pheno2$Stage)


# compare subgroups
# create column
pheno2$group_simple <- pheno2$group_in_paper

# Combine all NASH stages into one label
pheno2$group_simple[pheno2$group_in_paper %in% c("NASH_F0-F1", "NASH_F2", "NASH_F3", "NASH_F4")] <- "NASH"
pheno2$group_simple <- recode(pheno2$group_simple, "control" = "Control")

unique(pheno2$group_simple)



dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData   = pheno,
  design    = ~ group_simple   # your grouping variable
)

dds <- DESeq(dds)
res_NAFL <- results(dds, contrast = c("group_simple", "NAFL", "Control"))
res_NASH <- results(dds, contrast = c("group_simple", "NASH", "Control"))

res_NAFL["ENSG00000116044", ]
res_NASH["ENSG00000116044", ]

# HMOX1
res_NAFL["ENSG00000100292", ]
res_NASH["ENSG00000100292", ]

# SLC7A11
res_NAFL["ENSG00000151012", ]
res_NASH["ENSG00000151012", ]
# GPX4
res_NAFL["ENSG00000167468", ]
res_NASH["ENSG00000167468", ]

# FTH1
res_NAFL["ENSG00000167996", ]
res_NASH["ENSG00000167996", ]
# FTL
res_NAFL["ENSG00000087086", ]
res_NASH["ENSG00000087086", ]
# LOX
res_NAFL["ENSG00000113083", ]
res_NASH["ENSG00000113083", ]
# ACSL4
res_NAFL["ENSG00000068366", ]
res_NASH["ENSG00000068366", ]















# GSTA1
res_NAFL["ENSG00000243955", ]
res_NASH["ENSG00000243955", ]
# GSTA2
res_NAFL["ENSG00000244067", ]
res_NASH["ENSG00000244067", ]
# GSTM1
res_NAFL["ENSG00000134184", ]
res_NASH["ENSG00000134184", ]
# GSTM3
res_NAFL["ENSG00000134202", ]
res_NASH["ENSG00000134202", ]
# GSTM4
res_NAFL["ENSG00000168765", ]
res_NASH["ENSG00000168765", ]


j



dds2 <- DESeqDataSetFromMatrix(
  countData = count_matrix_2,
  colData   = pheno2,
  design    = ~ group_in_paper   # your grouping variable
)

dds2 <- DESeq(dds2)


res_NAFL <- results(dds2, contrast = c("group_in_paper", "NAFL", "control"))
res_F0F1 <- results(dds2, contrast = c("group_in_paper", "NASH_F0-F1", "control"))
res_F2   <- results(dds2, contrast = c("group_in_paper", "NASH_F2", "control"))
res_F3   <- results(dds2, contrast = c("group_in_paper", "NASH_F3", "control"))
res_F4   <- results(dds2, contrast = c("group_in_paper", "NASH_F4", "control"))


# SREBF target genes
srebp1_ensembl <- c(
  "ENSG00000169710",  # FASN
  "ENSG00000168477",  # SCD
  "ENSG00000196924",  # DGAT1
  "ENSG00000133136"   # DGAT2
)


res_NAFL[srebp1_ensembl, ]
res_F0F1[srebp1_ensembl, ]
res_F2[srebp1_ensembl, ]
res_F3[srebp1_ensembl, ]
res_F4[srebp1_ensembl, ]


# SREBP
res_NAFL["ENSG00000072310", ]
res_F0F1["ENSG00000072310", ]
res_F2["ENSG00000072310", ]
res_F3["ENSG00000072310", ]
res_F4["ENSG00000072310", ]

# ChREBP
res_NAFL["ENSG00000009950", ]
res_F0F1["ENSG00000009950", ]
res_F2["ENSG00000009950", ]
res_F3["ENSG00000009950", ]
res_F4["ENSG00000009950", ]


# check high activity NASH vs low activity NASH
# Convert NAS score to numeric
pheno$nas_score <- as.numeric(pheno$nas_score)






# Subset only NASH samples
pheno_nash <- pheno[pheno$group_simple == "NASH", ]

# Create high-activity (NAS >= 6) and low-activity (NAS <= 4) groups
pheno_nash$activity_group <- NA
pheno_nash$activity_group[pheno_nash$nas_score >= 6] <- "NASH_high_activity"
pheno_nash$activity_group[pheno_nash$nas_score <= 4] <- "NASH_low_activity"

# Remove samples that don't fall into either group (NAS = 5)
pheno_nash <- pheno_nash[!is.na(pheno_nash$activity_group), ]

table(pheno_nash$activity_group)

# Keep only columns (samples) that are in pheno_nash
count_nash <- count_matrix[, pheno_nash$geo_accession]

all(colnames(count_nash) == pheno_nash$geo_accession)

pheno_nash$activity_group <- factor(pheno_nash$activity_group,
                                    levels = c("NASH_low_activity", "NASH_high_activity"))

dds <- DESeqDataSetFromMatrix(
  countData = count_nash,
  colData = pheno_nash,
  design = ~ activity_group
)

dds <- DESeq(dds)
res_high_activity_vs_low_activity <- results(dds, contrast = c("activity_group", "NASH_high_activity", "NASH_low_activity"))
 

# compare stage
dds_stage <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData   = pheno,
  design    = ~ Stage
)


dds_stage <- DESeq(dds_stage)
res_early_vs_control <- results(dds_stage, contrast = c("Stage", "early", "control"))
res_moderate_vs_control <- results(dds_stage, contrast = c("Stage", "moderate", "control"))

res_early_vs_control["ENSG00000116044", ]
res_moderate_vs_control["ENSG00000116044", ]

# create nrf2 gene set
nrf2_genes <- c(
  "GCLM","SRXN1","GCLC","GSR","NQO1","SLC7A11",
  "AKR1C3","ME1","FTH1","FTL","OSGIN1","EPHX1",
  "PIR","ABHD4"
)

# ROS and Xenobiotic Detoxification
ros_detox_genes <- c(
  "NQO1",
  "PTGR1",
  "CYP2A6",
  "GSTA1", "GSTA2", "GSTA3", "GSTA5",
  "GSTP1", "GSTM1", "GSTM3", "MGST1",
  "UGT1A1", "UGT1A4", "UGT1A8", "UGT1A10",
  "CES1", "CBR1", "CBR3",
  "PRDX6", "PRDX1"
)


# GSH Production and Regeneration
GSH_genes <- c("GCLC", "GCLM", "GSR", "SLC7A11")

# Heme and Iron Metabolism
heme_genes <- c("FTL", "FTH1", "HMOX1")

# NADPH Regeneration
nadph_genes <- c("G6PD", "PGD", "TKT", "TALDO1", "ME1", "IDH1")

# Thioredoxin System
thioredoxin_genes <- c("TXN", "SRXN1", "TXNRD1")

# combine all the objects into a single vector
nrf2_genes <- c(
  ros_detox_genes,
  GSH_genes,
  heme_genes,
  nadph_genes,
  thioredoxin_genes
)

nrf2_sets <- data.frame(
  term = rep("NRF2_core", length(nrf2_genes)),
  gene = nrf2_genes
)


# create function
run_gsea_nrf2_custom <- function(res_df, nrf2_sets_df) {
  
  # convert rownames(ENSEMBL) → SYMBOL
  gene_symbols <- mapIds(
    org.Hs.eg.db,
    keys = rownames(res_df),
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  
  # attach SYMBOLs to res_df
  res_df$symbol <- gene_symbols
  
  # remove NA symbols
  res_df <- res_df[!is.na(res_df$symbol), ]
  
  # rank by log2FC using SYMBOL as names
  res_ranked <- res_df[order(res_df$log2FoldChange, decreasing = TRUE), ]
  gene_ranks <- res_ranked$log2FoldChange
  names(gene_ranks) <- res_ranked$symbol
  
  # remove duplicates (same symbol) keep max logFC
  gene_ranks <- tapply(gene_ranks, names(gene_ranks), max) |> sort(decreasing = TRUE)
  
  # print overlap
  overlap_genes <- intersect(names(gene_ranks), nrf2_sets_df$gene)
  cat("Number of overlapping genes:", length(overlap_genes), "\n")
  print(overlap_genes)
  
  # run GSEA
  gsea_result <- GSEA(
    geneList = gene_ranks,
    TERM2GENE = nrf2_sets_df,
    minGSSize = 5,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    verbose = FALSE
  )
  
  return(gsea_result)
}


# run gsea
gsea_NAFL     <- run_gsea_nrf2_custom(res_NAFL, nrf2_sets)
gsea_NASH     <- run_gsea_nrf2_custom(res_NASH, nrf2_sets)



# check individual genes
# check individual gene results
extract_nrf2_results <- function(res_df, nrf2_genes) {
  
  # map ENSEMBL rownames → SYMBOL temporarily
  gene_symbols <- mapIds(
    org.Hs.eg.db,
    keys = rownames(res_df),
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  
  res_df$symbol <- gene_symbols
  
  # now subset using SYMBOL not rownames
  res_sub <- res_df[res_df$symbol %in% nrf2_genes, ]
  
  # order by logFC
  res_sub <- res_sub[order(res_sub$log2FoldChange, decreasing = TRUE), ]
  
  # replace rownames with symbols (for display only)
  rownames(res_sub) <- res_sub$symbol
  
  return(res_sub)
}

nrf2_NAFL <- extract_nrf2_results(res_NAFL, nrf2_genes)
nrf2_NASH <- extract_nrf2_results(res_NASH, nrf2_genes)

nrf2_NAFL <- nrf2_NAFL[nrf2_genes, ]
nrf2_NASH <- nrf2_NASH[nrf2_genes, ]


rm(list=ls())
write.csv(nrf2_NAFL,
          file = "Sally_GSE135251_nrf2_NAFL_results.csv",
          row.names = TRUE)

write.csv(nrf2_NASH,
          file = "Sally_GSE135251_nrf2_NASH_results.csv",
          row.names = TRUE)

# NASH subgroup comparisons
nrf2_F0F1 <- extract_nrf2_results(res_F0F1, nrf2_genes)
nrf2_F2 <- extract_nrf2_results(res_F2, nrf2_genes)
nrf2_F3 <- extract_nrf2_results(res_F3, nrf2_genes)
nrf2_F4 <- extract_nrf2_results(res_F4, nrf2_genes)

nrf2_F0F1 <- nrf2_F0F1[nrf2_genes, ]
nrf2_F2 <- nrf2_F2[nrf2_genes, ]
nrf2_F3 <- nrf2_F3[nrf2_genes, ]
nrf2_F4 <- nrf2_F4[nrf2_genes, ]

nrf2_F4

write.csv(nrf2_F0F1,
          file = "Sally_GSE135251_nrf2_F0F1_results.csv",
          row.names = TRUE)

write.csv(nrf2_F2,
          file = "Sally_GSE135251_nrf2_F2_results.csv",
          row.names = TRUE)

write.csv(nrf2_F3,
          file = "Sally_GSE135251_nrf2_F3_results.csv",
          row.names = TRUE)

write.csv(nrf2_F4,
          file = "Sally_GSE135251_nrf2_F4_results.csv",
          row.names = TRUE)

nrf2_high_activity_vs_low_activity <- extract_nrf2_results(res_high_activity_vs_low_activity, nrf2_genes)
nrf2_high_activity_vs_low_activity <- nrf2_high_activity_vs_low_activity[nrf2_genes, ]

write.csv(nrf2_high_activity_vs_low_activity,
          file = "Sally_GSE135251_nrf2_high_activity_vs_low_activity_results.csv",
          row.names = TRUE)



# === correlate NRF2 expression with fibrosis + NAS ===

# normalize counts
vsd <- vst(dds, blind = TRUE)
expr_mat <- assay(vsd)

# NRF2 expression using ENSEMBL ID
nrf2_expr <- expr_mat["ENSG00000116044", ]


# make numeric (because they are character factors in this GEO)
pheno$fibrosis_stage <- as.numeric(pheno$fibrosis_stage)
pheno$nas_score      <- as.numeric(pheno$nas_score)

# remove samples with NA so cor.test doesn't break
valid_idx <- !is.na(pheno$fibrosis_stage) & !is.na(nrf2_expr)
valid_idx2 <- !is.na(pheno$nas_score) & !is.na(nrf2_expr)

# correlation spearman
cat("\nNRF2 vs fibrosis stage:\n")
print(cor.test(nrf2_expr[valid_idx], pheno$fibrosis_stage[valid_idx], method="spearman"))

cat("\nNRF2 vs NAS score:\n")
print(cor.test(nrf2_expr[valid_idx2], pheno$nas_score[valid_idx2], method="spearman"))




















# filter for genes that had a p value<0.05 and LFC > 1
sig_genes <- res %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%   # gene IDs in a column
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  arrange(padj, desc(abs(log2FoldChange)))  # most significant + strongest effect first

sig_genes_early_vs_control <- res_early_vs_control %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%   # gene IDs in a column
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  arrange(padj, desc(abs(log2FoldChange)))  # most significant + strongest effect first

nrow(sig_genes_early_vs_control)

sig_genes_moderate_vs_control <- res_moderate_vs_control %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%   # gene IDs in a column
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  arrange(padj, desc(abs(log2FoldChange)))  # most significant + strongest effect first

nrow(sig_genes_moderate_vs_control)

# how many significant DEGs
nrow(sig_genes)  


# export the DEGs
write.table(
  sig_genes$gene,                     # Data frame to export
  file = "sig_genesGSE135251.tsv",  # File name (you can change this)
  sep = "\t",                     # Tab-separated
  row.names = FALSE,              # Don’t include row numbers
  quote = FALSE                   # Don’t put quotes around text
)

# visualization
# volcano plot
library(EnhancedVolcano)

# genes ot be annotated on the volcano plot
topGenes <- rownames(resLFC)[resLFC$padj < 0.05 & abs(resLFC$log2FoldChange) > 1]

EnhancedVolcano(
  resLFC,
  lab = rownames(resLFC),
  x = "log2FoldChange",
  y = "padj",
  xlim = c(-6, 6),
  ylim = c(0, 25),
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 2.0,
  colAlpha = 0.8,
  title = "NAFLD vs Control Volcano Plot",
  subtitle = "DESeq2 results using apeglm shrinkage",
  caption = "Source: resLFC from lfcShrink()",
  selectLab = topGenes,
  labSize = 2      # <- smaller font size for labels (default is 4)
)


# MA plot
plotMA(res, ylim=c(-5,5))


# heat map
library(pheatmap)


top_genes <- rownames(sig_genes)[1:50]  # top 50 by padj
pheatmap(assay(dds)[top_genes,], scale="row", annotation_col=pheno)


# functional analysis
library("biomartr")
library(clusterProfiler)
library(org.Hs.eg.db)


# read in differentially expressed genes
diff_genes <- read_delim(
  file = "C:/Users/chukw/OneDrive/Desktop/capstone project/Rstudio capstone/capstone/sig_genesGSE135251.tsv",
  delim = "\t"
)

# Map ENSEMBL IDs to Entrez IDs
diff_genes_entrez <- bitr(diff_genes$x,   # <--- note the $x
                         fromType = "ENSEMBL",
                         toType = "ENTREZID",
                         OrgDb = org.Hs.eg.db)$ENTREZID



# create gene universe that contains all analysed genes
gene_universe <- rownames(res)

# export list fo genes
write.table(
  gene_universe,                     # Data frame to export
  file = "gene_universeGSE135251.tsv",  # File name (you can change this)
  sep = "\t",                     # Tab-separated
  row.names = FALSE,              # Don’t include row numbers
  quote = FALSE                   # Don’t put quotes around text
)

# import gene_universe
gene_universe <- read_delim(
  file = "C:\\Users\\chukw\\OneDrive\\Desktop\\capstone project\\Rstudio capstone\\capstone\\gene_universeGSE135251.tsv",
  delim = "\t"
)


# Map ENSEMBL IDs to Entrez IDs
gene_universe_entrez <- bitr(gene_universe$x, 
                             fromType = "ENSEMBL", 
                             toType = "ENTREZID", 
                             OrgDb = org.Hs.eg.db)$ENTREZID


# enrichment analysis
ora_analysis_bp <- enrichGO(gene          = diff_genes_entrez,
                universe      = gene_universe_entrez,
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP",         # Biological Process
                pAdjustMethod = "BH",
                qvalueCutoff  = 0.05,
                readable      = TRUE)


ora_analysis_mf <- enrichGO(
  gene          = diff_genes_entrez,
  universe      = gene_universe_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",  # Molecular Function
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

ora_analysis_cc <- enrichGO(
  gene          = diff_genes_entrez,
  universe      = gene_universe_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",   # Cellular Component
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)


# simplify proccesses that have similar names
ora_analysis_bp_simplified <- clusterProfiler::simplify(ora_analysis_bp) 
ora_analysis_mf_simplified <- clusterProfiler::simplify(ora_analysis_mf) 
ora_analysis_cc_simplified <- clusterProfiler::simplify(ora_analysis_cc) 

# export data
write_delim(
  x    = as.data.frame(ora_analysis_bp_simplified@result),
  file = "go_results_BP.tsv",
  delim = "\t"
)

write_delim(
  x    = as.data.frame(ora_analysis_mf_simplified@result),
  file = "go_results_MF.tsv",
  delim = "\t"
)

write_delim(
  x    = as.data.frame(ora_analysis_cc_simplified@result),
  file = "go_results_CC.tsv",
  delim = "\t"
)

dotplot(ora_analysis_bp_simplified)



# gene set enrichment analysis
# Load mapping libraries
library(org.Hs.eg.db)
library(dplyr)
library(clusterProfiler)

# sort DESeq2 results from most relevant to least relevant
res <- res[order(res$stat, decreasing = TRUE), ]
gene_ranks <- res$stat
names(gene_ranks) <- rownames(res)
gene_ranks <- gene_ranks[!is.na(gene_ranks)]

# Step 2: Convert Ensembl IDs to Gene Symbols
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = names(gene_ranks),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# keep only mapped genes
gene_ranks <- gene_ranks[!is.na(gene_symbols)]
gene_symbols <- gene_symbols[!is.na(gene_symbols)]

# Step 4: Replace Ensembl IDs with uppercase gene symbols
names(gene_ranks) <- toupper(gene_symbols)

# Remove duplicates by keeping the maximum stat for each gene
gene_ranks <- tapply(gene_ranks, names(gene_ranks), max)
gene_ranks <- sort(setNames(as.numeric(gene_ranks), names(gene_ranks)), decreasing = TRUE)
gene_ranks <- sort(gene_ranks, decreasing = TRUE)  # sort again



nrf2_sets <- read.gmt("C:/Users/chukw/OneDrive/Desktop/capstone project/NRF2_pathway_from_KEG.gmt")


gsea_result <- GSEA(
  geneList = gene_ranks,
  TERM2GENE = nrf2_sets,  
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

# making threshold less strict 0.05 to 0.3
gsea_result <- GSEA(
  geneList = gene_ranks,
  TERM2GENE = nrf2_sets,
  minGSSize = 5,
  maxGSSize = 500,
  pvalueCutoff = 0.3,   # less strict
  verbose = FALSE
)
# still there was no statistically significant upregulation or downregulation of nrf pathway

# checkign if the genes in the nrf2 target dataset are mostly in the ranked genes list
length(intersect(names(gene_ranks), nrf2_sets$gene))
intersect(names(gene_ranks), nrf2_sets$gene)
# they were mostly int he ranked genes list so that doesnt explain the lack of statistical significance

# view the ranking fo individual genes in the nrf2 pathway to try get an answer
nrf2_ranks <- gene_ranks[names(gene_ranks) %in% nrf2_sets$gene]
sort(nrf2_ranks, decreasing = TRUE)
# nrf2 pathway may nto have showed up as statistically significant due to there being genes in the pathway that were strongly upregulated vs genes in th epathway that were weakly upregulated so at the pathway level there wasnt upregulation or downregulation




