rm(list=ls())

# Load required libraries
library(limma)
library(data.table)
library(stringr)
library(GEOquery)
library(dplyr)
library(biomaRt)
library(clusterProfiler)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(fgsea)
library(readr)
library(biomartr)
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



# load data from folder
folder_path <- "C:/Users/chukw/OneDrive/Desktop/capstone project/Rstudio capstone/capstone/GSE33575_RAW"
files <- list.files(path = folder_path, pattern = "\\.txt$", full.names = TRUE)


agilent_data <- read.maimages(files, source = "agilent", green.only = TRUE)


# background correction
agilent_data <- backgroundCorrect(agilent_data, method = "normexp", offset = 50)

# check if it is already log2 transformed
summary(as.vector(agilent_data$E))


# log transform
exprs_log <- log2(agilent_data$E + 1)
# 2. Normalize between arrays
exprs_norm <- normalizeBetweenArrays(exprs_log, method = "quantile")

gene_symbols <- toupper(agilent_data$genes$GeneName)

keep <- !is.na(gene_symbols) & gene_symbols != ""

exprs_filt <- exprs_norm[keep, ]
genes_filt <- gene_symbols[keep]

exprs_gene <- avereps(exprs_filt, ID = genes_filt)

# load pheno data
gse <- getGEO("GSE33575", GSEMatrix = TRUE)
pheno <- pData(phenoData(gse[[1]]))
head(pheno)

# match GSM IDs between expression data and phenotype data
gsm_ids <- str_extract(colnames(agilent_data$E), "GSM\\d+")
pheno <- pheno[match(gsm_ids, pheno$geo_accession), ]
all(gsm_ids %in% pheno$geo_accession) 

# define experimental groups
group <- factor(pheno$`genotype:ch1`, levels = c("wild type", "Nrf2 KO"))
levels(group) <- make.names(levels(group))  # Converts to "wild.type" and "Nrf2.KO"

# create design matrix
design <- model.matrix(~0 + group)
colnames(design) <- c("wild.type", "Nrf2.KO")


fit <- lmFit(exprs_gene, design)
contrast.matrix <- makeContrasts(Nrf2_vs_WT = Nrf2.KO - wild.type, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# get full results (GENE-LEVEL)
results_gene <- topTable(
  fit2,
  coef = "Nrf2_vs_WT",
  adjust.method = "BH",
  number = Inf
)

# Add GeneSymbol column
results_gene$GeneSymbol <- rownames(results_gene)



# 4. Example: pull individual genes
# Example: pull G6PD
G6PD_results <- results_gene[results_gene$GeneSymbol == "G6PDX", ]
G6PD_results




# 6. Function to extract NRF2 target gene set
extract_nrf2_results <- function(res_df, nrf2_genes) {
  # Subset NRF2 genes
  res_sub <- res_df[res_df$GeneSymbol %in% nrf2_genes, ]
  
  # Order by logFC descending
  res_sub <- res_sub[order(res_sub$logFC, decreasing = TRUE), ]
  
  return(res_sub)
}

nrf2_results <- extract_nrf2_results(results_gene, nrf2_genes)

nrf2_results_ordered <- nrf2_results[match(nrf2_genes, nrf2_results$GeneSymbol), ]
nrf2_results_ordered <- nrf2_results_ordered %>%
  dplyr::select(GeneSymbol, dplyr::everything())

nrf2_results_ordered <- nrf2_results_ordered %>%
  dplyr::select(GeneSymbol, logFC, AveExpr, t, P.Value, adj.P.Val, B, everything())

write.csv(
  nrf2_results_ordered,
  file = "GSE33575_HFD_NRF2_KO_vs_WT_averaged.csv",
  row.names = FALSE
)













#------- enrichment analysis
# Use original gene names (no toupper)
gene_symbols_raw <- agilent_data$genes$GeneName

keep_raw <- !is.na(gene_symbols_raw) & gene_symbols_raw != ""
exprs_filt_raw <- exprs_norm[keep_raw, ]
genes_filt_raw <- gene_symbols_raw[keep_raw]

# Average probes by gene symbol (proper mouse case)
exprs_gene_raw <- avereps(exprs_filt_raw, ID = genes_filt_raw)


fit_raw <- lmFit(exprs_gene_raw, design)
fit2_raw <- contrasts.fit(fit_raw, contrast.matrix)
fit2_raw <- eBayes(fit2_raw)


results_gene_raw <- topTable(
  fit2_raw,
  coef = "Nrf2_vs_WT",
  adjust.method = "BH",
  number = Inf
)

# Add proper mouse gene symbols
results_gene_raw$GeneSymbol <- rownames(results_gene_raw)

DEGs_up_raw <- subset(results_gene_raw, adj.P.Val < 0.05 & logFC > 1)
DEGs_down_raw <- subset(results_gene_raw, adj.P.Val < 0.05 & logFC < -1)


library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)


symbols_up_raw   <- DEGs_up_raw$GeneSymbol
symbols_down_raw <- DEGs_down_raw$GeneSymbol



entrez_up_raw <- mapIds(
  org.Mm.eg.db,
  keys = symbols_up_raw,
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)

entrez_down_raw <- mapIds(
  org.Mm.eg.db,
  keys = symbols_down_raw,
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)

# Remove NAs
entrez_up_raw   <- na.omit(entrez_up_raw)
entrez_down_raw <- na.omit(entrez_down_raw)

nfe2l2_entrez <- mapIds(
  org.Mm.eg.db,
  keys = "Nfe2l2",
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)

entrez_down_noNrf2 <- setdiff(entrez_down_raw, nfe2l2_entrez)


ego_up_raw <- enrichGO(
  gene          = entrez_up_raw,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

ego_down_noNrf2_raw <- enrichGO(
  gene          = entrez_down_noNrf2,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

ekegg_up_raw <- enrichKEGG(
  gene          = entrez_up_raw,
  organism      = "mmu",
  keyType       = "ncbi-geneid",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

ekegg_down_noNrf2_raw <- enrichKEGG(
  gene          = entrez_down_noNrf2,
  organism      = "mmu",
  keyType       = "ncbi-geneid",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

ego_up_raw_simpl <- simplify(ego_up_raw, cutoff = 0.7, by = "p.adjust", select_fun = min)
ego_down_noNrf2_raw_simpl <- simplify(ego_down_noNrf2_raw, cutoff = 0.7, by = "p.adjust", select_fun = min)

df_ego_up_raw <- as.data.frame(ego_up_raw_simpl@result)
df_ego_down_noNrf2_raw <- as.data.frame(ego_down_noNrf2_raw_simpl@result)

df_kegg_up_raw <- as.data.frame(ekegg_up_raw@result)
df_kegg_down_noNrf2_raw <- as.data.frame(ekegg_down_noNrf2_raw@result)

# ✅ Filter KEGG results for adjusted p-value <= 0.05
df_kegg_up_raw <- df_kegg_up_raw %>% filter(p.adjust <= 0.05)
df_kegg_down_noNrf2_raw <- df_kegg_down_noNrf2_raw %>% filter(p.adjust <= 0.05)


# Add GeneSymbols column for upregulated GO
df_ego_up_raw$GeneSymbols <- sapply(strsplit(df_ego_up_raw$geneID, "/"), function(ids) {
  paste(mapIds(org.Mm.eg.db, ids, column="SYMBOL", keytype="ENTREZID", multiVals="first"), collapse="/")
})

# Add GeneSymbols column for downregulated GO (Nfe2l2 excluded)
df_ego_down_noNrf2_raw$GeneSymbols <- sapply(strsplit(df_ego_down_noNrf2_raw$geneID, "/"), function(ids) {
  paste(mapIds(org.Mm.eg.db, ids, column="SYMBOL", keytype="ENTREZID", multiVals="first"), collapse="/")
})

# Add GeneSymbols column for upregulated KEGG
df_kegg_up_raw$GeneSymbols <- sapply(strsplit(df_kegg_up_raw$geneID, "/"), function(ids) {
  paste(mapIds(org.Mm.eg.db, ids, column="SYMBOL", keytype="ENTREZID", multiVals="first"), collapse="/")
})

# Add GeneSymbols column for downregulated KEGG (Nfe2l2 excluded)
df_kegg_down_noNrf2_raw$GeneSymbols <- sapply(strsplit(df_kegg_down_noNrf2_raw$geneID, "/"), function(ids) {
  paste(mapIds(org.Mm.eg.db, ids, column="SYMBOL", keytype="ENTREZID", multiVals="first"), collapse="/")
})

nrow(df_ego_up_raw)
nrow(df_ego_down_noNrf2_raw)

nrow(df_kegg_up_raw)
nrow(df_kegg_down_noNrf2_raw)

top20_go_down_raw <- df_ego_down_noNrf2_raw %>%
  arrange(p.adjust) %>%
  slice(1:20) %>%
  mutate(Description = factor(Description, levels = rev(Description)))

top20_go_down_raw$GeneRatio_num <- sapply(top20_go_down_raw$GeneRatio, function(x) {
  parts <- strsplit(x, "/")[[1]]
  as.numeric(parts[1]) / as.numeric(parts[2])
})


go_down_plot <- ggplot(top20_go_down_raw, aes(x = GeneRatio_num, y = Description)) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_gradient(low = "red", high = "blue") +
  labs(
    x = "Gene Ratio",
    y = "GO Biological Process",
    color = "Adjusted p-value",
    size = "Gene Count",
    title = "Top 20 Downregulated GO BP (Nfe2l2 Excluded)"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  )

# KEG graph
top20_kegg_down_raw <- df_kegg_down_noNrf2_raw %>%
  arrange(p.adjust) %>%
  slice(1:20) %>%
  mutate(Description = factor(Description, levels = rev(Description)))

top20_kegg_down_raw$GeneRatio_num <- sapply(top20_kegg_down_raw$GeneRatio, function(x) {
  parts <- strsplit(x, "/")[[1]]
  as.numeric(parts[1]) / as.numeric(parts[2])
})

top20_kegg_down_raw$Description <- gsub(
  " - Mus musculus \\(house mouse\\)",
  "",
  top20_kegg_down_raw$Description
)


kegg_plot <- ggplot(top20_kegg_down_raw, aes(x = GeneRatio_num, y = Description)) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_gradient(low = "red", high = "blue") +
  labs(
    x = "Gene Ratio",
    y = "KEGG Pathway",
    color = "Adjusted p-value",
    size = "Gene Count",
    title = "Top 20 Downregulated KEGG Pathways (Nfe2l2 Excluded)"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  )








top20_go_up_raw <- df_ego_up_raw %>%
  arrange(p.adjust) %>%
  slice(1:20) %>%
  mutate(Description = factor(Description, levels = rev(Description)))

top20_go_up_raw$GeneRatio_num <- sapply(top20_go_up_raw$GeneRatio, function(x) {
  parts <- strsplit(x, "/")[[1]]
  as.numeric(parts[1]) / as.numeric(parts[2])
})

go_up_plot <- ggplot(top20_go_up_raw, aes(x = GeneRatio_num, y = Description)) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_gradient(low = "red", high = "blue") +
  labs(
    x = "Gene Ratio",
    y = "GO Biological Process",
    color = "Adjusted p-value",
    size = "Gene Count",
    title = "Top 20 Upregulated GO BP"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  )








top20_kegg_up_raw <- df_kegg_up_raw %>%
  arrange(p.adjust) %>%
  slice(1:20) %>%
  mutate(Description = factor(Description, levels = rev(Description)))

top20_kegg_up_raw$GeneRatio_num <- sapply(top20_kegg_up_raw$GeneRatio, function(x) {
  parts <- strsplit(x, "/")[[1]]
  as.numeric(parts[1]) / as.numeric(parts[2])
})

top20_kegg_up_raw$Description <- gsub(
  " - Mus musculus \\(house mouse\\)",
  "",
  top20_kegg_up_raw$Description
)

kegg_up_plot <- ggplot(top20_kegg_up_raw, aes(x = GeneRatio_num, y = Description)) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_gradient(low = "red", high = "blue") +
  labs(
    x = "Gene Ratio",
    y = "KEGG Pathway",
    color = "Adjusted p-value",
    size = "Gene Count",
    title = "Top 20 Upregulated KEGG Pathways"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  )


kegg_up_plot




# Export GO enrichment results
write.csv(df_ego_up_raw, "GSE33575_GO_up_raw.csv", row.names = FALSE)
write.csv(df_ego_down_noNrf2_raw, "GSE33575_GO_down_raw.csv", row.names = FALSE)

# Export KEGG enrichment results
write.csv(df_kegg_up_raw, "GSE33575_KEGG_up_raw.csv", row.names = FALSE)
write.csv(df_kegg_down_noNrf2_raw, "GSE33575_KEGG_down_raw.csv", row.names = FALSE)

go_up_df_GSE33575 <- read.csv(
  "C:/Users/chukw/OneDrive/Desktop/capstone project/Rstudio capstone/capstone/GSE33575_GO_up_raw.csv",
  stringsAsFactors = FALSE
)

go_down_df_GSE33575 <- read.csv(
  "C:/Users/chukw/OneDrive/Desktop/capstone project/Rstudio capstone/capstone/GSE33575_GO_down_raw.csv",
  stringsAsFactors = FALSE
)

kegg_up_df_GSE33575 <- read.csv(
  "C:/Users/chukw/OneDrive/Desktop/capstone project/Rstudio capstone/capstone/GSE33575_KEGG_up_raw.csv",
  stringsAsFactors = FALSE
)

kegg_down_df_GSE33575 <- read.csv(
  "C:/Users/chukw/OneDrive/Desktop/capstone project/Rstudio capstone/capstone/GSE33575_KEGG_down_raw.csv",
  stringsAsFactors = FALSE
)

go_up_df_GSE33575
go_down_df_GSE33575

kegg_up_df_GSE33575
kegg_down_df_GSE33575

nrow(go_up_df_GSE33575)
nrow(go_down_df_GSE33575)
nrow(kegg_up_df_GSE33575)
nrow(kegg_down_df_GSE33575)

jjk



library(patchwork)  # for combining ggplots

# ----------------------
# Combine GO plots vertically
go_combined <- go_up_plot / go_down_plot + 
  plot_annotation(title = "GO Enrichment (Up vs Down)")

# Combine KEGG plots vertically
kegg_combined <- kegg_up_plot / kegg_plot + 
  plot_annotation(title = "KEGG Enrichment (Up vs Down)")

# Combine GO and KEGG side by side
combined_all <- go_combined | kegg_combined



# ----------------------
# Export combined image
ggsave(
  filename = "GO_KEGG_combined.png",
  plot = combined_all,
  width = 16, height = 12, dpi = 300
)

# Export individual plots
ggsave("GO_up_plot.png", plot = go_up_plot, width = 8, height = 6, dpi = 300)
ggsave("GO_down_plot.png", plot = go_down_plot, width = 8, height = 6, dpi = 300)
ggsave("KEGG_up_plot.png", plot = kegg_up_plot, width = 8, height = 6, dpi = 300)
ggsave("KEGG_down_plot.png", plot = kegg_plot, width = 8, height = 6, dpi = 300)


jjk



DEGs <- all_genes_annot %>%
  filter(adj.P.Val < 0.05, abs(logFC) >= 1)

head(DEGs)


# map the mouose symbols also use version 105 as later versions run into issues
mouse <- useEnsembl(
  biomart = "ensembl",
  dataset = "mmusculus_gene_ensembl",
  version = 105
)

human <- useEnsembl(
  biomart = "ensembl",
  dataset = "hsapiens_gene_ensembl",
  version = 105
)

# Map mouse symbols to human symbols
mapping_all <- getLDS(
  attributes = "mgi_symbol",
  filters = "mgi_symbol",
  values = all_genes_annot$GeneName,   # <- all genes, not DEGs
  mart = mouse,
  attributesL = "hgnc_symbol",
  martL = human,
  uniqueRows = TRUE
)

# Merge DEGs with mapping
all_genes_human <- merge(all_genes_annot, mapping_all,
                         by.x = "GeneName", by.y = "MGI.symbol", all.x = TRUE)

# keep only mapped DEGs
all_genes_human_mapped <- all_genes_human[!is.na(all_genes_human$HGNC.symbol), ]

# --- Create ranked gene list for GSEA
gene_ranks <- all_genes_human_mapped$logFC   # or t
names(gene_ranks) <- all_genes_human_mapped$HGNC.symbol

# Handle duplicates by keeping the probe with the strongest (max) t
gene_ranks <- tapply(gene_ranks, names(gene_ranks), max)
gene_ranks <- sort(gene_ranks, decreasing = TRUE)

# load nrf2 gene set
nrf2_sets <- read.gmt("C:/Users/chukw/OneDrive/Desktop/capstone project/NRF2_pathway_from_KEG.gmt")


# Convert names to uppercase to match the gene set
names(gene_ranks) <- toupper(names(gene_ranks))
nrf2_sets$gene <- toupper(nrf2_sets$gene)

# Sort descending (required for GSEA)
gene_ranks <- sort(gene_ranks, decreasing = TRUE)

# (Optional) Check overlap before GSEA
nrf2_genes <- unique(nrf2_sets$gene)
overlap_genes <- intersect(names(gene_ranks), nrf2_genes)
length(overlap_genes)
overlap_genes

# --- Run GSEA ---
gsea_result <- GSEA(
  geneList = gene_ranks,
  TERM2GENE = nrf2_sets,
  minGSSize = 5,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

# View full result table
as.data.frame(gsea_result)




nrf2_ranks <- gene_ranks[names(gene_ranks) %in% nrf2_sets$gene]
nrf2_ranks <- sort(nrf2_ranks, decreasing = TRUE)
nrf2_ranks


# rerunning th egsea withotu NRF2 to see if NRF2 deletion isnt just skewing the result downward
# Remove NFE2L2 from the Nrf2 gene set
nrf2_sets_noNFE2L2 <- nrf2_sets[nrf2_sets$gene != "NFE2L2", ]

# Rerun GSEA
gsea_result_noNFE2L2 <- GSEA(
  geneList = gene_ranks,
  TERM2GENE = nrf2_sets_noNFE2L2,
  minGSSize = 5,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

as.data.frame(gsea_result_noNFE2L2)

DEGs <- all_genes_annot %>%
  filter(adj.P.Val < 0.05, abs(logFC) >= 1)


DEGs_entrez <- bitr(DEGs$GeneName,
                    fromType = "SYMBOL",
                    toType   = "ENTREZID",
                    OrgDb    = org.Mm.eg.db)


gene_universe_entrez <- bitr(all_genes_annot$GeneName,
                             fromType = "SYMBOL",
                             toType   = "ENTREZID",
                             OrgDb    = org.Mm.eg.db)



ora_bp <- enrichGO(
  gene          = DEGs_entrez$ENTREZID,
  universe      = gene_universe_entrez$ENTREZID,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

ora_mf <- enrichGO(
  gene          = DEGs_entrez$ENTREZID,
  universe      = gene_universe_entrez$ENTREZID,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

ora_cc <- enrichGO(
  gene          = DEGs_entrez$ENTREZID,
  universe      = gene_universe_entrez$ENTREZID,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# --- Simplify results (remove redundant terms) ---
ora_bp_simple <- simplify(ora_bp)
ora_mf_simple <- simplify(ora_mf)
ora_cc_simple <- simplify(ora_cc)

nrow(ora_bp_simple@result)
nrow(ora_mf_simple@result)
nrow(ora_cc_simple@result)

head(ora_bp_simple@result, 10)
head(ora_mf_simple@result, 10)
head(ora_cc_simple@result, 10)


# Split DEGs into upregulated and downregulated
DEGs_up <- DEGs %>% filter(logFC > 1)
DEGs_down <- DEGs %>% filter(logFC < -1)

# map to entrezIDs
DEGs_up_entrez <- bitr(DEGs_up$Gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Mm.eg.db)
DEGs_down_entrez <- bitr(DEGs_down$Gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Mm.eg.db)

# Upregulated pathways
ora_bp_up <- enrichGO(
  gene = DEGs_up_entrez$ENTREZID,
  universe = gene_universe_entrez$ENTREZID,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

# Downregulated pathways
ora_bp_down <- enrichGO(
  gene = DEGs_down_entrez$ENTREZID,
  universe = gene_universe_entrez$ENTREZID,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

df_up <- as.data.frame(ora_bp_up)
df_down <- as.data.frame(ora_bp_down)


write.csv(df_up, "ORA_BP_upregulated_GO_results.csv", row.names = FALSE)
write.csv(df_down, "ORA_BP_downregulated_GO_results.csv", row.names = FALSE)

df_up <- read.csv("C:/Users/chukw/OneDrive/Desktop/capstone project/Rstudio capstone/capstone/ORA_BP_upregulated_GO_results.csv")
df_down <- read.csv("C:/Users/chukw/OneDrive/Desktop/capstone project/Rstudio capstone/capstone/ORA_BP_downregulated_GO_results.csv")

df_down[df_down$ID == "GO:0000303", ]




genes_of_interest <- c("GCLC", "GCLM", "GSR", "SLC7A11")


# pllot simplified results
# simplify
ora_bp_up_simple <- simplify(ora_bp_up)
ora_bp_down_simple <- simplify(ora_bp_down)

ora_bp_up_simple_data_frame <- as.data.frame(ora_bp_up_simple)
ora_bp_down_simple_data_frame <- as.data.frame(ora_bp_down_simple)

write.csv(ora_bp_up_simple_data_frame, "ORA_BP_upregulated_GO_results_simple.csv", row.names = FALSE)
write.csv(ora_bp_down_simple_data_frame, "ORA_BP_downregulated_GO_results_simple.csv", row.names = FALSE)


# graph of simplified
p_up <- dotplot(
  ora_bp_up_simple,
  showCategory = 15,
  title = "(A) Upregulated GO Biological Processes"
) +
  scale_size(range = c(1, 4))

p_down <- dotplot(
  ora_bp_down_simple,
  showCategory = 15,
  title = "(B) Downregulated GO Biological Processes"
) +
  scale_size(range = c(1, 4))

p_up
p_down

# put images together
p_up / p_down






# convert ot data frame
up_results <- as.data.frame(ora_bp_up_simple@result)
down_results <- as.data.frame(ora_bp_down_simple@result)

down_results$Description





# make boxplots
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
# 2. Prepare expression data
# -------------------------------
expr_df <- as.data.frame(exprs_norm)
expr_df$ProbeID <- rownames(expr_df)

# Join GeneSymbol from annotation
expr_df <- expr_df %>%
  left_join(
    results_annot_noNA %>% dplyr::select(ProbeID, GeneSymbol),
    by = "ProbeID"
  ) %>%
  filter(!is.na(GeneSymbol) & GeneSymbol != "")

# -------------------------------
# 3. Pivot to long format
# -------------------------------
expr_long <- expr_df %>%
  pivot_longer(
    cols = -c(ProbeID, GeneSymbol),
    names_to = "Sample",
    values_to = "Expression"
  )

# -------------------------------
# 4. Extract GSM IDs from file paths
# -------------------------------
expr_long$Sample_GSM <- str_extract(expr_long$Sample, "GSM\\d+")

# Make sure pheno rownames are GSM IDs
rownames(pheno) <- pheno$geo_accession

# Assign group based on GSM ID
expr_long$Group <- factor(
  pheno$`genotype:ch1`[match(expr_long$Sample_GSM, rownames(pheno))],
  levels = c("wild type", "Nrf2 KO")
)

# -------------------------------
# 5. Collapse multiple probes to mean per gene per sample
# -------------------------------
expr_gene <- expr_long %>%
  group_by(GeneSymbol, Group, Sample_GSM) %>%
  summarise(Expression = mean(Expression, na.rm = TRUE), .groups = "drop")

# -------------------------------
# 6. Filter for genes of interest
# -------------------------------
plot_df <- expr_gene %>%
  filter(GeneSymbol %in% genes_of_interest)

# -------------------------------
# 7. Create boxplots
# -------------------------------
library(ggplot2)

ggplot(plot_df, aes(x = Group, y = Expression)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", size = 4, shape = 18) +
  stat_summary(
    fun = mean,
    geom = "text",
    aes(label = round(after_stat(y), 2)),
    vjust = -1.2,
    color = "red",
    size = 3
  ) +
  facet_wrap(~ GeneSymbol, scales = "free_y") +
  theme_bw() +
  labs(
    title = "Expression of NADPH-Regenerating Genes (Collapsed by Gene)",
    y = "Log2 Normalized Expression",
    x = "Genotype"
  )

