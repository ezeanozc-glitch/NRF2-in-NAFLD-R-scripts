rm(list=ls())
# microarray with nrf2 KO in liver regeneration
library(affy)
library(limma)
library(GEOquery)
library(dplyr)
library(annotate)
library(mouse4302.db)
library(biomaRt)
library(ggplot2)
library(patchwork)


# fatty acid synthesis genes
nrf2_genes <- c(
  # Transcriptional regulators
  "SREBF1",    # SREBP-1c
  "MLXIPL",    # ChREBP
  "LXR",       # Liver X receptor
  
  # Key enzymatic genes
  "ACLY",      # ATP citrate lyase
  "ACSS2",     # Acetyl-CoA synthetase 2
  "ACACA",     # Acetyl-CoA carboxylase alpha
  "FASN",      # Fatty acid synthase
  "ELOVL6",    # Elongation of very long chain fatty acids protein 6
  "SCD1"      # Stearoyl-CoA desaturase 1
)

folder_path <- "C:/Users/chukw/OneDrive/Desktop/capstone project/Rstudio capstone/capstone/GSE8969_RAW"
list.files(folder_path)

# read data in
data <- ReadAffy(celfile.path = folder_path)

# get sample names
sampleNames(data)

# RMA preprocessing: background correction + normalization + summarization
eset <- rma(data)

# Check expression matrix
exprs_mat <- exprs(eset)
colnames(exprs_mat) <- sub("\\.CEL$", "", colnames(exprs_mat))


gse <- getGEO("GSE8969", GSEMatrix = TRUE)
pheno <- pData(phenoData(gse[[1]]))
head(pheno)

colnames(pheno)[colnames(pheno) == "genotype:ch1"] <- "genotype"
unique(pheno$genotype)

# check if sample names match
all(colnames(exprs_mat) == rownames(pheno))

# Extract genotype information
group <- factor(pheno$genotype)
levels(group) <- make.names(levels(group))


design <- model.matrix(~0 + group)
colnames(design) <- levels(group)


contrast_matrix <- makeContrasts(
  KO_vs_WT = Nrf2.ko - Nrf2.wt,
  levels = design
)


probe2gene <- getSYMBOL(rownames(exprs_mat), "mouse4302.db")
keep <- !is.na(probe2gene)

exprs_filt <- exprs_mat[keep, ]
genes_filt <- toupper(probe2gene[keep])

exprs_gene <- avereps(exprs_filt, ID = genes_filt)

fit <- lmFit(exprs_gene, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

results_gene <- topTable(
  fit2,
  coef = "KO_vs_WT",
  adjust.method = "BH",
  number = Inf
)

# Add GeneSymbol column back
results_gene$GeneSymbol <- rownames(results_gene)





# G6PD
G6PD_results <- results_gene[results_gene$GeneSymbol == "G6PDX", ]
G6PD_results


#NFE2L2
nfe2l2_results <- results_annot_noNA[results_annot_noNA$GeneSymbol == "NFE2L2", ]
nfe2l2_results




# check individual genes
extract_nrf2_results <- function(res_df, nrf2_genes) {
  # Subset rows matching nrf2_genes
  res_sub <- res_df[res_df$GeneSymbol %in% nrf2_genes, ]
  
  # Order by logFC
  res_sub <- res_sub[order(res_sub$logFC, decreasing = TRUE), ]
  
  return(res_sub)
}


nrf2_results <- extract_nrf2_results(results_gene, nrf2_genes)

nrf2_results_ordered <- nrf2_results[match(nrf2_genes, nrf2_results$GeneSymbol), ]

nrf2_results_ordered <- nrf2_results_ordered %>%
  dplyr::select(GeneSymbol, dplyr::everything())

head(nrf2_results_ordered)

write.csv(nrf2_results_ordered,
          file = "GSE8969_NRF2_KO_basal_averaged.csv",
          row.names = FALSE)

















# make box plots
genes_of_interest <- c("G6PDX","PGD","TKT","TALDO1","ME1","IDH1")

genes_of_interest <- c("GCLC", "GCLM", "GSR", "SLC7A11")

genes_of_interest <- c(
  "GSTA1", "GSTA2", "GSTA3", "GSTA5",
  "GSTM1", "GSTM3"
)

genes_of_interest <- c(
  "FTL1", "FTH1", "HMOX1"
)


# -------------------------------
# 3. Prepare expression dataframe
# -------------------------------
expr_df <- as.data.frame(exprs_mat)  # normalized expression
expr_df$ProbeID <- rownames(expr_df)

# Join gene symbols from annotation
expr_df <- expr_df %>%
  left_join(
    results_annot_noNA %>% dplyr::select(ProbeID, GeneSymbol),
    by = "ProbeID"
  ) %>%
  filter(!is.na(GeneSymbol) & GeneSymbol != "" & GeneSymbol %in% genes_of_interest)

# -------------------------------
# 4. Pivot to long format
# -------------------------------
expr_long <- expr_df %>%
  pivot_longer(
    cols = -c(ProbeID, GeneSymbol),
    names_to = "Sample",
    values_to = "Expression"
  )

# -------------------------------
# 5. Map sample to genotype
# -------------------------------
expr_long$Group <- factor(
  pheno$genotype[match(expr_long$Sample, rownames(pheno))],
  levels = c("Nrf2 wt", "Nrf2 ko")  # wt first, ko second
)

# -------------------------------
# 6. Aggregate probes per gene (optional)
# -------------------------------
expr_long <- expr_long %>%
  group_by(GeneSymbol, Sample, Group) %>%
  summarise(Expression = mean(Expression), .groups = "drop")


ggplot(expr_long, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(aes(color = Group), width = 0.1, size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "red") +
  facet_wrap(~GeneSymbol, scales = "free_y") +   # <- separate panel per gene
  labs(
    x = "Genotype",
    y = "Normalized expression (log2)",
    fill = "Genotype",
    color = "Genotype",
    title = "Expression of NADPH-regenerating genes in Nrf2 KO vs WT livers"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

















#------- over representation analsysis
# 1. Filter DEGs
DEGs_up <- results_gene[results_gene$logFC > 1 & results_gene$adj.P.Val < 0.05, ]
DEGs_down <- results_gene[results_gene$logFC < -1 & results_gene$adj.P.Val < 0.05, ]


library(clusterProfiler)
library(org.Mm.eg.db)





# Keep original gene symbols from annotation (no toupper)
genes_filt_DEG <- probe2gene[keep]

# Collapse probes to gene-level expression
exprs_gene_DEG <- avereps(exprs_filt, ID = genes_filt_DEG)

# Fit linear model and contrast
fit_DEG <- lmFit(exprs_gene_DEG, design)
fit2_DEG <- contrasts.fit(fit_DEG, contrast_matrix)
fit2_DEG <- eBayes(fit2_DEG)

# Get DEGs
results_gene_DEG <- topTable(
  fit2_DEG,
  coef = "KO_vs_WT",
  adjust.method = "BH",
  number = Inf
)

# Add GeneSymbol column back
results_gene_DEG$GeneSymbol <- rownames(results_gene_DEG)

# Filter for significant DEGs (logFC > 1 or < -1, adj.P.Val < 0.05)
DEGs_up_DEG <- results_gene_DEG[results_gene_DEG$logFC > 1 & results_gene_DEG$adj.P.Val < 0.05, ]
DEGs_down_DEG <- results_gene_DEG[results_gene_DEG$logFC < -1 & results_gene_DEG$adj.P.Val < 0.05, ]

# 1. Extract SYMBOL -> ENTREZID mapping from mouse4302.db
probe2entrez <- as.list(mouse4302ENTREZID)
probe2symbol <- as.list(mouse4302SYMBOL)



















# 2. Map upregulated DEGs
symbols_up <- DEGs_up_DEG$GeneSymbol
entrez_up_DEG <- sapply(symbols_up, function(s) {
  # get probes matching this symbol
  probes <- names(probe2symbol)[probe2symbol == s]
  # get Entrez IDs for those probes
  entrez_ids <- unlist(probe2entrez[probes])
  unique(entrez_ids)
})
entrez_up_DEG <- unique(unlist(entrez_up_DEG))  # final vector of Entrez IDs

# 3. Map downregulated DEGs
symbols_down <- DEGs_down_DEG$GeneSymbol
entrez_down_DEG <- sapply(symbols_down, function(s) {
  probes <- names(probe2symbol)[probe2symbol == s]
  entrez_ids <- unlist(probe2entrez[probes])
  unique(entrez_ids)
})
entrez_down_DEG <- unique(unlist(entrez_down_DEG))  # final vector of Entrez IDs

mapIds(org.Mm.eg.db, "Nfe2l2", "ENTREZID", "SYMBOL")

nfe2l2_entrez <- mapIds(org.Mm.eg.db, "Nfe2l2", "ENTREZID", "SYMBOL")

entrez_down_noNrf2 <- setdiff(entrez_down_DEG, nfe2l2_entrez)

# Check how many genes mapped
length(entrez_up_DEG)
length(entrez_down_noNrf2)


# Optional: see first few Entrez IDs
head(entrez_up_DEG)
head(entrez_down_noNrf2)











ego_up <- enrichGO(
  gene = entrez_up_DEG,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

ego_down_noNrf2 <- enrichGO(
  gene = entrez_down_noNrf2,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
















# Simplify upregulated GO results
ego_up_simpl <- simplify(
  ego_up,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)

ego_down_noNrf2_simpl <- simplify(
  ego_down_noNrf2,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)



# View top simplified results
head(ego_up_simpl)
head(ego_down_noNrf2_simpl)

nrow(ego_up_simpl)
nrow(ego_down_noNrf2_simpl)













# Upregulated
ego_up_simpl@result$geneSymbols <- sapply(
  ego_up_simpl@result$geneID,
  function(x) paste(
    mapIds(
      org.Mm.eg.db,
      keys = strsplit(x, "/")[[1]],
      column = "SYMBOL",
      keytype = "ENTREZID",
      multiVals = "first"
    ),
    collapse = "/"
  )
)

# Downregulated (NRF2 excluded)
ego_down_noNrf2_simpl@result$geneSymbols <- sapply(
  ego_down_noNrf2_simpl@result$geneID,
  function(x) paste(
    mapIds(
      org.Mm.eg.db,
      keys = strsplit(x, "/")[[1]],
      column = "SYMBOL",
      keytype = "ENTREZID",
      multiVals = "first"
    ),
    collapse = "/"
  )
)



df_ego_up <- as.data.frame(ego_up_simpl@result)
df_ego_down_noNrf2 <- as.data.frame(ego_down_noNrf2_simpl@result)

df_ego_up_clean <- df_ego_up[, c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "Count", "geneSymbols")]
df_ego_down_noNrf2_clean <- df_ego_down_noNrf2[, c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "Count", "geneSymbols")]

write.csv(df_ego_up_clean, "GSE8969_GO_BP_up_DEGs_simplified.csv", row.names = FALSE)
write.csv(df_ego_down_noNrf2_clean, "GSE8969_GO_BP_down_DEGs_noNrf2_simplified.csv", row.names = FALSE)






# Top 20 upregulated
top20_up <- df_ego_up_clean %>%
  arrange(p.adjust) %>%
  slice(1:20) %>%
  mutate(Description = factor(Description, levels = rev(Description)))


top20_up <- top20_up %>%
  mutate(GeneRatio = sapply(GeneRatio, function(x) {
    nums <- as.numeric(unlist(strsplit(x, "/")))
    nums[1] / nums[2]
  }))


ggplot(top20_up, aes(x = GeneRatio, y = Description)) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_gradient(low = "red", high = "blue") +
  labs(
    x = "Gene Ratio",
    y = "GO Biological Process",
    color = "Adjusted p-value",
    size = "Gene Count",
    title = "Top 20 Upregulated GO Pathways"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))


# Top 20 downregulated (NRF2 excluded)
top20_down_noNrf2 <- df_ego_down_noNrf2_clean %>%
  arrange(p.adjust) %>%
  slice(1:20) %>%
  mutate(Description = factor(Description, levels = rev(Description)))

top20_down_noNrf2 <- top20_down_noNrf2 %>%
  mutate(GeneRatio = sapply(GeneRatio, function(x) {
    nums <- as.numeric(unlist(strsplit(x, "/")))
    nums[1] / nums[2]
  }))

go_down_plot

go_down_plot <- ggplot(top20_down_noNrf2, aes(x = GeneRatio, y = Description)) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_gradient(low = "red", high = "blue") +
  labs(
    x = "Gene Ratio",
    y = "GO Biological Process",
    color = "Adjusted p-value",
    size = "Gene Count",
    title = "B Top 20 Downregulated GO Pathways (NRF2 Excluded)"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))

top20_up
top20_down_noNrf2











#------- KEG over representation analysis
ekegg_down_noNrf2 <- enrichKEGG(
  gene          = entrez_down_noNrf2,
  organism      = "mmu",
  keyType       = "ncbi-geneid",      # <-- FIXED
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)


# 1️⃣ Convert KEGG results to a dataframe
df_kegg_down_noNrf2 <- as.data.frame(ekegg_down_noNrf2@result)

# Check first few rows
head(df_kegg_down_noNrf2)

nrow(df_kegg_down_noNrf2)

# 2️⃣ Select columns of interest
df_kegg_down_noNrf2_clean <- df_kegg_down_noNrf2[, c(
  "ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "Count", "geneID"
)]

df_kegg_down_noNrf2_clean <- df_kegg_down_noNrf2_clean %>%
  filter(p.adjust <= 0.05)



# 3️⃣ Map gene IDs to symbols
df_kegg_down_noNrf2_clean$geneSymbols <- sapply(
  df_kegg_down_noNrf2_clean$geneID,
  function(x) paste(
    mapIds(
      org.Mm.eg.db,
      keys = strsplit(x, "/")[[1]],
      column = "SYMBOL",
      keytype = "ENTREZID",
      multiVals = "first"
    ),
    collapse = "/"
  )
)


# 4️⃣ Export results to CSV
write.csv(df_kegg_down_noNrf2_clean, "GSE8969_KEGG_down_DEGs_noNrf2.csv", row.names = FALSE)




















top20_kegg_down <- df_kegg_down_noNrf2_clean %>%
  arrange(p.adjust) %>%
  slice(1:20) %>%
  mutate(Description = factor(Description, levels = rev(Description)))

# Clean the Description column to remove " - Mus musculus (house mouse)"
top20_kegg_down <- top20_kegg_down %>%
  mutate(
    Description = gsub(" - Mus musculus \\(house mouse\\)", "", Description),
    Description = factor(Description, levels = rev(Description))  # keep the factor ordering
  )


# Convert GeneRatio from "a/b" to decimal
top20_kegg_down <- top20_kegg_down %>%
  mutate(
    GeneRatio = sapply(GeneRatio, function(x) {
      nums <- as.numeric(strsplit(x, "/")[[1]])
      nums[1] / nums[2]
    })
  )

kegg_plot

kegg_plot <- ggplot(top20_kegg_down, aes(x = GeneRatio, y = Description)) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_gradient(low = "red", high = "blue") +
  labs(
    x = "Gene Ratio",
    y = "KEGG Pathway",
    color = "Adjusted p-value",
    size = "Gene Count",
    title = "A Top 14 Downregulated KEGG Pathways (NRF2 Excluded)"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  )


ekegg_down_noNrf2@result

# Combine plots vertically
combined_plot <- kegg_plot / go_down_plot



combined_plot


# export
ggsave("GSE8969_combined_plot.pdf", combined_plot, width = 8, height = 12)

ggsave("GSE8969_GO_down_plot.pdf", go_down_plot, width = 8, height = 6)

ggsave("GSE8969_KEGG_plot.pdf", kegg_plot, width = 8, height = 6)







go_down_df <- read.csv("C:/Users/chukw/OneDrive/Desktop/capstone project/Rstudio capstone/capstone/GSE8969_GO_BP_down_DEGs_noNrf2_simplified.csv")
go_down_df <- read.csv("C:/Users/chukw/OneDrive/Desktop/capstone project/Rstudio capstone/capstone/GSE8969_GO_BP_down_DEGs_noNrf2_simplified.csv")

kegg_down_df <- read.csv("C:/Users/chukw/OneDrive/Desktop/capstone project/Rstudio capstone/capstone/GSE8969_KEGG_down_DEGs_noNrf2.csv")

kegg_down_df
go_down_df

nrow(kegg_down_df)
nrow(go_down_df)


