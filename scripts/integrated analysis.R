library(biomaRt)
library(dplyr)
library(tidyr)
library(ggplot2)

# integrated analysis of dataset 1 and dataset 2
# find common genes between dataset 1 and 2
common_genes <- intersect(
  rownames(counts_matrix_1),
  rownames(count_matrix_2)
)

length(common_genes)

# subset count matrixes
counts_matrix_1 <- counts_matrix_1[common_genes, ]
count_matrix_2  <- count_matrix_2[common_genes, ]

# add condition and batch information
pheno1$condition <- factor(
  pheno1$condition,
  levels = c("Control", "NAFL", "NASH")
)
pheno1$batch <- "GSE126848"

pheno2$condition <- factor(
  pheno2$group_simple,
  levels = c("Control", "NAFL", "NASH")
)
pheno2$batch <- "GSE135251"

# merge count matrixes
counts_merged <- cbind(
  counts_matrix_1,
  count_matrix_2
)

# merge metadata
meta_merged <- bind_rows(
  pheno1[, c("condition", "batch")],
  pheno2[, c("condition", "batch")]
)
rownames(meta_merged) <- colnames(counts_merged)

dds <- DESeqDataSetFromMatrix(
  countData = counts_merged,
  colData   = meta_merged,
  design    = ~ batch + condition
)


dds <- DESeq(dds)

res_NAFL_vs_Control <- results(
  dds,
  contrast = c("condition", "NAFL", "Control")
)

res_NASH_vs_Control <- results(
  dds,
  contrast = c("condition", "NASH", "Control")
)

# NFE2L2
res_NAFL_vs_Control["ENSG00000116044", ]
res_NASH_vs_Control["ENSG00000116044", ]


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

# optional individual gene comparison for obeese
nrf2_nafld  <- extract_nrf2_results(res_NAFL_vs_Control, nrf2_genes)
nrf2_nash   <- extract_nrf2_results(res_NASH_vs_Control,  nrf2_genes)

nrf2_nafld <- nrf2_nafld[nrf2_genes, ]
nrf2_nash <- nrf2_nash[nrf2_genes, ]

write.csv(nrf2_nafld,
          file = "Dataset1_and_2_pooled_results_NAFL.csv",
          row.names = TRUE)

write.csv(nrf2_nash,
          file = "Dataset1_and_2_pooled_results_NASH.csv",
          row.names = TRUE)



#----- make graphs
genes_of_interest <- c(
  "GSTA1", "GSTA2", "GSTA3", "GSTA5",
  "GSTP1", "GSTM1", "GSTM3"
)

genes_of_interest <- c(
  "FTL", "FTH1", "HMOX1"
)

genes_of_interest <- c(
  "G6PD", "PGD", "TKT", "TALDO1",
  "ME1", "IDH1"
)

genes_of_interest <- c(
  "GCLC", "GCLM", "GSR", "SLC7A11"
)

# get normalized counts
norm_counts <- counts(dds, normalized = TRUE)

mart <- useMart(
  biomart = "ensembl",
  dataset = "hsapiens_gene_ensembl"
)

gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = rownames(norm_counts),
  mart = mart
)

norm_df <- as.data.frame(norm_counts)
norm_df$ensembl_gene_id <- rownames(norm_df)

norm_df <- norm_df %>%
  left_join(gene_map, by = "ensembl_gene_id") %>%
  filter(hgnc_symbol != "")


norm_subset <- norm_df %>%
  filter(hgnc_symbol %in% genes_of_interest)

unique(norm_subset$hgnc_symbol)


plot_df <- norm_subset %>%
  dplyr::select(-ensembl_gene_id) %>%
  pivot_longer(
    cols = -hgnc_symbol,
    names_to = "sample",
    values_to = "expression"
  ) %>%
  left_join(
    meta_merged %>% mutate(sample = rownames(meta_merged)),
    by = "sample"
  )


ggplot(plot_df, aes(x = condition, y = expression)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  facet_wrap(~ hgnc_symbol, scales = "free_y") +
  scale_y_log10() +
  labs(
    x = "Condition",
    y = "Normalized expression (log10)",
    title = "Expression of GST family genes across disease states"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "italic"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


# create graphs with different stages of fibrosis
pheno2$group_plot <- pheno2$group_in_paper

pheno2$group_plot <- factor(
  pheno2$group_plot,
  levels = c(
    "Control",
    "NAFL",
    "NASH F0–F1",
    "NASH F2",
    "NASH F3",
    "NASH F4"
  )
)

norm_counts <- counts(dds, normalized = TRUE)

mart <- useMart(
  biomart = "ensembl",
  dataset = "hsapiens_gene_ensembl"
)

gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = rownames(norm_counts),
  mart = mart
)

norm_df <- as.data.frame(norm_counts)
norm_df$ensembl_gene_id <- rownames(norm_df)

norm_df <- norm_df %>%
  left_join(gene_map, by = "ensembl_gene_id")


hmox_df <- norm_df %>%
  filter(hgnc_symbol == "HMOX1") %>%
  dplyr::select(-ensembl_gene_id) %>%
  pivot_longer(
    cols = -hgnc_symbol,
    names_to = "sample",
    values_to = "expression"
  ) %>%
  left_join(
    pheno2 %>% mutate(sample = rownames(pheno2)),
    by = "sample"
  )

ggplot(hmox_df, aes(x = group_plot, y = expression)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90") +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  scale_y_log10() +
  labs(
    x = "Disease stage",
    y = "Normalized HMOX1 expression (log10)",
    title = "HMOX1 expression across NAFLD and NASH fibrosis stages"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )
