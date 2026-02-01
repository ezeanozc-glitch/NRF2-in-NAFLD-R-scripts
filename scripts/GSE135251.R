# load libraries
library(tidyverse)
library(GEOquery)
library(DESeq2)
library(org.Hs.eg.db)
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

# Convert dataframe to a matrix
count_matrix_2 <- as.matrix(count_matrix_2)

# Convert to integer matrix
mode(count_matrix_2) <- "integer"

# stripping uneccessary text
colnames(count_matrix_2) <- gsub("\\.counts\\.txt$", "", colnames(count_matrix_2))
colnames(count_matrix_2) <- gsub("(_.*)", "", colnames(count_matrix_2))  # optional: keep only GSM ID

# load metadata
gse2 <- getGEO("GSE135251", GSEMatrix = TRUE)
pheno2 <- pData(phenoData(gse2[[1]]))
head(pheno2)

# subset data
pheno2 <- pheno2[match(colnames(count_matrix_2), pheno2$geo_accession), ]

# check if the sample IDs in the column names of the count matrix matches the GEO accession values in the pheno metadata
all(colnames(count_matrix_2) %in% pheno2$geo_accession)


# remove symbols
colnames(pheno2)[colnames(pheno2) == "disease:ch1"] <- "disease"
colnames(pheno2)[colnames(pheno2) == "fibrosis stage:ch1"] <- "fibrosis_stage"
colnames(pheno2)[colnames(pheno2) == "group in paper:ch1"] <- "group_in_paper"
colnames(pheno2)[colnames(pheno2) == "nas score:ch1"] <- "nas_score"
colnames(pheno2)[colnames(pheno2) == "Stage:ch1"] <- "Stage"

# view contents
unique(pheno2$disease)
unique(pheno2$fibrosis_stage)
unique(pheno2$group_in_paper)
unique(pheno2$nas_score)
unique(pheno2$Stage)


# create column
pheno2$group_simple <- pheno2$group_in_paper

# Combine all NASH stages into one label
pheno2$group_simple[pheno2$group_in_paper %in% c("NASH_F0-F1", "NASH_F2", "NASH_F3", "NASH_F4")] <- "NASH"
pheno2$group_simple <- recode(pheno2$group_simple, "control" = "Control")

unique(pheno2$group_simple)


# DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData   = pheno2,
  design    = ~ group_simple   # your grouping variable
)

# run differential expression analysis
dds <- DESeq(dds)
res_NAFL <- results(dds, contrast = c("group_simple", "NAFL", "Control"))
res_NASH <- results(dds, contrast = c("group_simple", "NASH", "Control"))



# create DESeq2
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






# compare stages
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



# create function to subset nrf2 target genes
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

# run
nrf2_NAFL <- extract_nrf2_results(res_NAFL, nrf2_genes)
nrf2_NASH <- extract_nrf2_results(res_NASH, nrf2_genes)

nrf2_NAFL <- nrf2_NAFL[nrf2_genes, ]
nrf2_NASH <- nrf2_NASH[nrf2_genes, ]


# export results
write.csv(nrf2_NAFL,
          file = "GSE135251_nrf2_NAFL.csv",
          row.names = TRUE)

write.csv(nrf2_NASH,
          file = "GSE135251_nrf2_NASH.csv",
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

# export results
write.csv(nrf2_F0F1,
          file = "GSE135251_nrf2_F0F1_results.csv",
          row.names = TRUE)

write.csv(nrf2_F2,
          file = "GSE135251_nrf2_F2_results.csv",
          row.names = TRUE)

write.csv(nrf2_F3,
          file = "GSE135251_nrf2_F3_results.csv",
          row.names = TRUE)

write.csv(nrf2_F4,
          file = "GSE135251_nrf2_F4_results.csv",
          row.names = TRUE)

