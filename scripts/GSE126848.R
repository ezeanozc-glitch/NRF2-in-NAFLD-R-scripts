# load libraries
library(tidyverse)
library(GEOquery)
library(org.Hs.eg.db)
library(DESeq2)

# load counts matrix
counts_matrix_1 <- read.delim("C:/Users/chukw/OneDrive/Desktop/capstone project/Rstudio capstone/capstone/GSE126848_Gene_counts_raw.txt/GSE126848_Gene_counts_raw.txt",
                     row.names = 1, check.names = FALSE)

# load metadata
gse1 <- getGEO("GSE126848", GSEMatrix = TRUE)
pheno1 <- pData(gse1[[1]])
head(pheno1)


# Remove leading zeros from counts column names
colnames(counts_matrix_1) <- sub("^0+", "", colnames(counts_matrix_1))

# match column names with the rownames
pheno1 <- pheno1[match(colnames(counts_matrix_1), pheno1$description), ]


# set the rownames to the description column
rownames(pheno1) <- pheno1$description

# check if the row names adn column names match
all(colnames(counts_matrix_1) == rownames(pheno1))  # should return TRUE


# remove the symbol from the column name
colnames(pheno1)[colnames(pheno1) == "disease:ch1"] <- "disease"
colnames(pheno1)[colnames(pheno1) == "gender:ch1"] <- "gender"

#  Combine healthy + obese into Control for iintegrative analysis
pheno1$condition <- pheno1$disease
pheno1$condition[pheno1$condition %in% c("healthy", "obese")] <- "Control"
pheno1$condition[pheno1$condition == "NAFLD"] <- "NAFL"
pheno1$condition <- factor(pheno1$condition, levels = c("Control", "NAFL", "NASH"))
unique(pheno1$condition)


#optional line of code
#-------- run this if you want to set the obeese as the reference
pheno1$disease <- factor(pheno1$disease, levels = c("obese", "healthy", "NAFLD", "NASH"))
pheno1$disease <- relevel(pheno1$disease, ref = "obese")

# create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix_1,
  colData = pheno1,
  design = ~ disease
)

# dds with covarriate
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = pheno,
  design = ~ gender + disease
)

# run differential expression analysis
dds <- DESeq(dds)
res_healthy_vs_obese <- results(dds, contrast = c("disease", "healthy", "obese"))
res_nafld_vs_healthy <- results(dds, contrast = c("disease", "NAFLD", "healthy"))
res_nash_vs_healthy <- results(dds, contrast = c("disease", "NASH", "healthy"))


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


# ------ create function to subset for nrf2 target genes
extract_nrf2_results <- function(res_df, nrf2_genes) {
  
  gene_symbols <- mapIds(
    org.Hs.eg.db,
    keys = rownames(res_df),
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  
  res_df$symbol <- gene_symbols
  
  res_sub <- res_df[res_df$symbol %in% nrf2_genes, ]
  
  res_sub <- res_sub[order(res_sub$log2FoldChange, decreasing = TRUE), ]
  
  rownames(res_sub) <- res_sub$symbol
  
  return(res_sub)
}

# run for each comparison
nrf2_obese  <- extract_nrf2_results(res_healthy_vs_obese, nrf2_genes)
nrf2_nafld  <- extract_nrf2_results(res_nafld_vs_healthy, nrf2_genes)
nrf2_nash   <- extract_nrf2_results(res_nash_vs_healthy,  nrf2_genes)


nrf2_obese <- nrf2_obese[nrf2_genes, ]
nrf2_nafld <- nrf2_nafld[nrf2_genes, ]
nrf2_nash <- nrf2_nash[nrf2_genes, ]

# export results
write.csv(nrf2_obese,
          file = "GSE126848_nrf2_obese",
          row.names = TRUE)

write.csv(nrf2_nafld,
          file = "GSE126848_nrf2_nafld",
          row.names = TRUE)

write.csv(nrf2_nash,
          file = "GSE126848_nrf2_nash",
          row.names = TRUE)



