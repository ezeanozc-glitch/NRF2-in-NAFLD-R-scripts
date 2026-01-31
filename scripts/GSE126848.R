library(readr)
library(tidyverse)
library(GEOquery)
library(GEOquery)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(DESeq2)


counts_matrix_1 <- read.delim("C:/Users/chukw/OneDrive/Desktop/capstone project/Rstudio capstone/capstone/GSE126848_Gene_counts_raw.txt/GSE126848_Gene_counts_raw.txt",
                     row.names = 1, check.names = FALSE)


gse1 <- getGEO("GSE126848", GSEMatrix = TRUE)
pheno1<- pData(phenoData(gse1[[1]]))
head(pheno1)


# Remove leading zeros from counts column names
colnames(counts_matrix_1) <- sub("^0+", "", colnames(counts_matrix_1))

# match column names with the rownames
pheno1 <- pheno1[match(colnames(counts_matrix_1), pheno1$description), ]


# set the rownames to the description column
rownames(pheno1) <- pheno1$description

# check if the rownmaes adn column names match
all(colnames(counts_matrix_1) == rownames(pheno1))  # should return TRUE


# remove the symbol from the column name
colnames(pheno1)[colnames(pheno1) == "disease:ch1"] <- "disease"
colnames(pheno1)[colnames(pheno1) == "gender:ch1"] <- "gender"

#  Combine healthy + obese into Control
pheno1$condition <- pheno1$disease
pheno1$condition[pheno1$condition %in% c("healthy", "obese")] <- "Control"
pheno1$condition[pheno1$condition == "NAFLD"] <- "NAFL"
pheno1$condition <- factor(pheno1$condition, levels = c("Control", "NAFL", "NASH"))
unique(pheno1$condition)


#optional line of code
#-------- run this if you want to set the obeese as the reference
pheno1$disease <- factor(pheno1$disease, levels = c("obese", "healthy", "NAFLD", "NASH"))
pheno1$disease <- relevel(pheno1$disease, ref = "obese")


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

unique(pheno1$disease)

dds <- DESeq(dds)
res_healthy_vs_obese <- results(dds, contrast = c("disease", "healthy", "obese"))
res_nafld_vs_healthy <- results(dds, contrast = c("disease", "NAFLD", "healthy"))
res_nash_vs_healthy <- results(dds, contrast = c("disease", "NASH", "healthy"))

# optional obeese comparison
res_healthy_vs_obese <- results(dds, contrast = c("disease", "healthy", "obese"))
res_nafld_vs_obese <- results(dds, contrast = c("disease", "NAFLD", "obese"))
res_nash_vs_obese <- results(dds, contrast = c("disease", "NASH", "obese"))

# SREBP
res_healthy_vs_obese["ENSG00000072310", ]
res_nafld_vs_healthy["ENSG00000072310", ]
res_nash_vs_healthy["ENSG00000072310", ]

# ChREBP
res_healthy_vs_obese["ENSG00000009950", ]
res_nafld_vs_healthy["ENSG00000009950", ]
res_nash_vs_healthy["ENSG00000009950", ]





res_healthy_vs_obese["ENSG00000116044", ]
res_nafld_vs_healthy["ENSG00000116044", ]
res_nash_vs_healthy["ENSG00000116044", ]

# GSTM1
res_healthy_vs_obese["ENSG00000134184", ]
res_nafld_vs_healthy["ENSG00000134184", ]
res_nash_vs_healthy["ENSG00000134184", ]
#GSTM2
res_healthy_vs_obese["ENSG00000213366", ]
res_nafld_vs_healthy["ENSG00000213366", ]
res_nash_vs_healthy["ENSG00000213366", ]
#GSTM3
res_healthy_vs_obese["ENSG00000134202", ]
res_nafld_vs_healthy["ENSG00000134202", ]
res_nash_vs_healthy["ENSG00000134202", ]
# GSTM4
res_healthy_vs_obese["ENSG00000168765", ]
res_nafld_vs_healthy["ENSG00000168765", ]
res_nash_vs_healthy["ENSG00000168765", ]
# GSTM4
res_healthy_vs_obese["ENSG00000116044", ]
res_nafld_vs_healthy["ENSG00000116044", ]
res_nash_vs_healthy["ENSG00000116044", ]
# GSTM5
res_healthy_vs_obese["ENSG00000134201", ]
res_nafld_vs_healthy["ENSG00000134201", ]
res_nash_vs_healthy["ENSG00000134201", ]


# optional comparison
res_healthy_vs_obese["ENSG00000116044", ]
res_nafld_vs_obese["ENSG00000116044", ]
res_nash_vs_obese["ENSG00000116044", ]



### --- NRF2 gene set ---
nrf2_genes <- c(
  "GCLM","SRXN1","GCLC","GSR","NQO1","SLC7A11",
  "AKR1C3","ME1","FTH1","FTL","OSGIN1","EPHX1",
  "PIR","ABHD4"
)

nrf2_sets <- data.frame(
  term = rep("NRF2_core", length(nrf2_genes)),
  gene = nrf2_genes
)


### ---- GSEA function ----
run_gsea_nrf2_custom <- function(res_df, nrf2_sets_df) {
  
  # convert ENSEMBL → SYMBOL
  gene_symbols <- mapIds(
    org.Hs.eg.db,
    keys = rownames(res_df),
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  
  res_df$symbol <- gene_symbols
  res_df <- res_df[!is.na(res_df$symbol), ]
  
  # rank
  res_ranked <- res_df[order(res_df$log2FoldChange, decreasing = TRUE), ]
  gene_ranks <- res_ranked$log2FoldChange
  names(gene_ranks) <- res_ranked$symbol
  
  # remove duplicates
  gene_ranks <- tapply(gene_ranks, names(gene_ranks), max) |> sort(decreasing = TRUE)
  
  # check overlap
  overlap_genes <- intersect(names(gene_ranks), nrf2_sets_df$gene)
  cat("Number of overlapping genes:", length(overlap_genes), "\n")
  print(overlap_genes)
  
  # GSEA
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

### --- RUN GSEA for all 3 contrasts ---
gsea_obese_vs_healthy  <- run_gsea_nrf2_custom(res_healthy_vs_obese,  nrf2_sets)
gsea_nafld_vs_healthy  <- run_gsea_nrf2_custom(res_nafld_vs_healthy,  nrf2_sets)
gsea_nash_vs_healthy   <- run_gsea_nrf2_custom(res_nash_vs_healthy,   nrf2_sets)





# individual genes
# ------ get single gene results
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


# view tables
nrf2_obese
nrf2_nafld
nrf2_nash



# optional
# optional individual gene comparison for obeese
nrf2_healthy  <- extract_nrf2_results(res_healthy_vs_obese, nrf2_genes)
nrf2_nafld  <- extract_nrf2_results(res_nafld_vs_obese, nrf2_genes)
nrf2_nash   <- extract_nrf2_results(res_nash_vs_obese,  nrf2_genes)

nrf2_healthy
nrf2_nafld
nrf2_nash

nrf2_healthy <- nrf2_healthy[nrf2_genes, ]
nrf2_nafld <- nrf2_nafld[nrf2_genes, ]
nrf2_nash <- nrf2_nash[nrf2_genes, ]


write.csv(nrf2_healthy,
          file = "Sally2_GSE126848_nrf2_obese_redone.csv",
          row.names = TRUE)

write.csv(nrf2_nafld,
          file = "Sally2_GSE126848_nrf2_nafld_redone.csv",
          row.names = TRUE)

write.csv(nrf2_nash,
          file = "Sally2_GSE126848_nrf2_nash_redone.csv",
          row.names = TRUE)





library(org.Hs.eg.db)
library(dplyr)
library(clusterProfiler)



res_nafld_vs_healthy <- res_nafld_vs_healthy[order(res_nafld_vs_healthy$stat, decreasing = TRUE), ]
gene_ranks <- res_nafld_vs_healthy$log2FoldChange
names(gene_ranks) <- rownames(res_nafld_vs_healthy)
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


nrf2_ranks <- gene_ranks[names(gene_ranks) %in% nrf2_sets$gene]
sort(nrf2_ranks, decreasing = TRUE)
