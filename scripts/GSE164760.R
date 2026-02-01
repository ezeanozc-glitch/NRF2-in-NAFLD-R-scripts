library(oligo)
library(hugene11sttranscriptcluster.db)
library(GEOquery)
library(limma)
library(dplyr)
library(AnnotationDbi)
library(ggplot2)
library(tidyr)


GSE48452
GSE164760  
GSE89632
GSE66676 


gse <- getGEO("GSE48452", GSEMatrix = TRUE)
pheno <- pData(phenoData(gse[[1]]))
head(pheno)


# load data from folder
folder_path <- "C:/Users/chukw/OneDrive/Desktop/capstone project/Rstudio capstone/capstone/GSE48452_RAW"
list.files(folder_path)

# load files
cel_files <- list.files(folder_path, pattern = "\\.CEL$", full.names = TRUE)
raw_data <- read.celfiles(cel_files)

# normalisation and background correction
norm_data <- rma(raw_data)

# extract expression matrix
exprs_mat <- exprs(norm_data)

# extract GSM IDs from expression matrix column names
exprs_samples <- sub("_.*", "", colnames(exprs_mat))
# phenotype sample names
pheno_samples <- rownames(pheno)

# check if they are identical and in the same order
all(exprs_samples == pheno_samples)

colnames(pheno)[colnames(pheno) == "group:ch1"] <- "group"
colnames(pheno)[colnames(pheno) == "fat:ch1"] <- "fat"
colnames(pheno)[colnames(pheno) == "inflammation:ch1"] <- "inflammation"
colnames(pheno)[colnames(pheno) == "fibrosis:ch1"] <- "fibrosis"
colnames(pheno)[colnames(pheno) == "nas:ch1"] <- "nas"

unique(pheno$group)
unique(pheno$fat)
unique(pheno$inflammation)
unique(pheno$fibrosis)
unique(pheno$nas)

group <- factor(pheno$group,
                levels = c("Control", "Healthy obese", "Steatosis", "Nash"))

# rename levels to syntactically valid names
levels(group) <- c("Control", "Healthy_Obese", "Steatosis", "Nash")

# make design matrix
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# fit model
fit <- lmFit(exprs_mat, design)
# create contrasts
contrasts <- makeContrasts(
  HealthyObese_vs_Control = Healthy_Obese - Control,
  Steatosis_vs_Control    = Steatosis - Control,
  Nash_vs_Control         = Nash - Control,
  levels = design
)

# Apply contrasts + empirical Bayes
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)

# Extract results
res_HealthyObese <- topTable(fit2, coef = "HealthyObese_vs_Control", number = Inf)
res_Steatosis    <- topTable(fit2, coef = "Steatosis_vs_Control", number = Inf)
res_Nash         <- topTable(fit2, coef = "Nash_vs_Control", number = Inf)


# add gene symbols
add_symbols <- function(res) {
  res$SYMBOL <- mapIds(
    hugene11sttranscriptcluster.db,
    keys = rownames(res),
    column = "SYMBOL",
    keytype = "PROBEID",
    multiVals = "first"
  )
  res
}

res_HealthyObese <- add_symbols(res_HealthyObese)
res_Steatosis <- add_symbols(res_Steatosis)
res_Nash  <- add_symbols(res_Nash)

res_HealthyObese_clean <- res_HealthyObese[!is.na(res_HealthyObese$SYMBOL), ]
res_Steatosis_clean    <- res_Steatosis[!is.na(res_Steatosis$SYMBOL), ]
res_Nash_clean <- res_Nash[!is.na(res_Nash$SYMBOL), ]

res_HealthyObese_clean[res_HealthyObese_clean$SYMBOL == "NFE2L2", ]
res_Steatosis_clean[res_Steatosis_clean$SYMBOL == "NFE2L2", ]
res_Nash_clean[res_Nash_clean$SYMBOL == "NFE2L2", ]

# SREBF1
res_HealthyObese_clean[res_HealthyObese_clean$SYMBOL == "SREBF1", ]
res_Steatosis_clean[res_Steatosis_clean$SYMBOL == "SREBF1", ]
res_Nash_clean[res_Nash_clean$SYMBOL == "SREBF1", ]

srebp1_targets <- c("FASN", "ACACA", "SCD", "ELOVL6", "DGAT1", "DGAT2")

res_HealthyObese_clean[res_HealthyObese_clean$SYMBOL %in% srebp1_targets, ]
res_Steatosis_clean[res_Steatosis_clean$SYMBOL %in% srebp1_targets, ]
res_Nash_clean[res_Nash_clean$SYMBOL %in% srebp1_targets, ]

# check NRF2 changes
res_HealthyObese_clean[res_HealthyObese_clean$SYMBOL == "NFE2L2", ]
res_Steatosis_clean[res_Steatosis_clean$SYMBOL == "NFE2L2", ]
res_Nash_clean[res_Nash_clean$SYMBOL == "NFE2L2", ]


# check results for individual genes
extract_nrf2_results <- function(res_df, nrf2_genes) {
  # Subset rows matching NRF2-related genes
  res_sub <- res_df[res_df$SYMBOL %in% nrf2_genes, ]
  
  # If no genes are found, return empty data frame safely
  if (nrow(res_sub) == 0) return(res_sub)
  
  # Collapse multiple probes per gene:
  # keep probe with largest absolute logFC
  res_sub <- res_sub %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::slice_max(order_by = abs(logFC), n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  # Order by logFC
  res_sub <- res_sub[order(res_sub$logFC, decreasing = TRUE), ]
  
  return(res_sub)
}

nrf2_HealthyObese <- extract_nrf2_results(res_HealthyObese_clean, nrf2_genes)
nrf2_Steatosis    <- extract_nrf2_results(res_Steatosis_clean, nrf2_genes)
nrf2_Nash  <- extract_nrf2_results(res_Nash_clean, nrf2_genes)

#reorder NRF2 results
nrf2_HealthyObese_ordered <-
  nrf2_HealthyObese[match(nrf2_genes, nrf2_HealthyObese$SYMBOL), ]

nrf2_HealthyObese_ordered <-
  nrf2_HealthyObese_ordered[
    , c("SYMBOL", setdiff(colnames(nrf2_HealthyObese_ordered), "SYMBOL"))
  ]

nrf2_Steatosis_ordered <-
  nrf2_Steatosis[match(nrf2_genes, nrf2_Steatosis$SYMBOL), ]

nrf2_Steatosis_ordered <-
  nrf2_Steatosis_ordered[
    , c("SYMBOL", setdiff(colnames(nrf2_Steatosis_ordered), "SYMBOL"))
  ]

nrf2_Nash_ordered <-
  nrf2_Nash[match(nrf2_genes, nrf2_Nash$SYMBOL), ]

nrf2_Nash_ordered <-
  nrf2_Nash_ordered[
    , c("SYMBOL", setdiff(colnames(nrf2_Nash_ordered), "SYMBOL"))
  ]

# export
write.csv(
  nrf2_HealthyObese_ordered,
  file = "GSE48452_NRF2_HealthyObese_vs_Control.csv",
  row.names = FALSE
)

write.csv(
  nrf2_Steatosis_ordered,
  file = "GSE48452_NRF2_Steatosis_vs_Control.csv",
  row.names = FALSE
)

write.csv(
  nrf2_Nash_ordered,
  file = "GSE48452_NRF2_Nash_vs_Control.csv",
  row.names = FALSE
)










#------ correlation of genes with histological scores
# Convert clinical scores to numeric
pheno$fat <- as.numeric(pheno$fat)
pheno$inflammation <- as.numeric(pheno$inflammation)
pheno$fibrosis <- as.numeric(pheno$fibrosis)
pheno$nas <- as.numeric(pheno$nas)

# Map probe IDs to gene symbols
probe2symbol <- mapIds(
  hugene11sttranscriptcluster.db,
  keys = rownames(exprs_mat),
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

# Add symbols to expression matrix
exprs_df <- as.data.frame(exprs_mat)
exprs_df$SYMBOL <- probe2symbol

# Keep only NRF2 genes
exprs_nrf2 <- exprs_df[exprs_df$SYMBOL %in% nrf2_genes, ]

# Remove probes without symbols
exprs_nrf2 <- exprs_nrf2[!is.na(exprs_nrf2$SYMBOL), ]

exprs_nrf2_gene <- exprs_nrf2 %>%
  group_by(SYMBOL) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

# Convert to matrix
exprs_gene_mat <- as.matrix(exprs_nrf2_gene[, -1])
rownames(exprs_gene_mat) <- exprs_nrf2_gene$SYMBOL

# Transpose so rows = samples, cols = genes
exprs_gene_mat_t <- t(exprs_gene_mat)


exprs_samples <- sub("_.*", "", colnames(exprs_mat))

# Strip CEL suffix from expression matrix sample names
rownames(exprs_gene_mat_t) <- sub("_.*", "", rownames(exprs_gene_mat_t))

pheno <- pheno[rownames(exprs_gene_mat_t), ]

# Check alignment
all(rownames(exprs_gene_mat_t) == rownames(pheno))


clinical_vars <- c("fat", "inflammation", "fibrosis", "nas")

cor_results <- data.frame()

for (gene in colnames(exprs_gene_mat_t)) {
  for (var in clinical_vars) {
    
    x <- exprs_gene_mat_t[, gene]
    y <- pheno[[var]]
    
    # Remove NAs
    valid <- complete.cases(x, y)
    
    if (sum(valid) >= 10) {  # minimum samples safeguard
      
      cor_test <- cor.test(
        x[valid],
        y[valid],
        method = "spearman"
      )
      
      cor_results <- rbind(
        cor_results,
        data.frame(
          Gene = gene,
          Clinical_Score = var,
          Spearman_rho = cor_test$estimate,
          P_value = cor_test$p.value,
          N = sum(valid)
        )
      )
    }
  }
}

cor_results$Adj_P_value <- p.adjust(cor_results$P_value, method = "BH")

head(cor_results)


write.csv(
  cor_results,
  file = "NRF2_gene_clinical_score_correlations.csv",
  row.names = FALSE
)







#--------- create graphs with single gene multiple pannels and dots
gene <- "ME1"

gene <- "NFE2L2"


df <- data.frame(
  Expression = exprs_gene_mat_t[, gene],
  fat = pheno$fat,
  inflammation = pheno$inflammation,
  fibrosis = pheno$fibrosis,
  nas = pheno$nas
)

df_long <- tidyr::pivot_longer(
  df,
  cols = -Expression,
  names_to = "Score",
  values_to = "Value"
)

ggplot(df_long, aes(x = Value, y = Expression)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ Score, scales = "free_x") +
  theme_bw() +
  labs(title = paste(gene, "expression vs clinical scores"))








#------ create correlation graph with multiple pannels and multiple genes but no dots
genes_to_plot <- c("GSTA1", "GSTA2", "GSTA3", "GSTA5",
                    "GSTM1", "GSTM3")

genes_to_plot <- c("FTL", "FTH1", "HMOX1")

genes_to_plot <- c("G6PD", "PGD", "TKT", "TALDO1",
                   "ME1", "IDH1")

genes_to_plot <- c("NFE2L2")
nrf2_genes <- c("NFE2L2")

df <- data.frame(
  Sample = rownames(exprs_gene_mat_t),
  pheno,
  exprs_gene_mat_t[, genes_to_plot]
)

df_long <- df %>%
  pivot_longer(
    cols = all_of(genes_to_plot),
    names_to = "Gene",
    values_to = "Expression"
  ) %>%
  pivot_longer(
    cols = c(fat, inflammation, fibrosis, nas),
    names_to = "Score",
    values_to = "Value"
  )


ggplot(df_long, aes(x = Value, y = Expression, color = Gene)) +
  geom_point(alpha = 0.6, size = 1) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  facet_wrap(~ Score, scales = "free_x") +
  theme_bw() +
  labs(
    title = "NRF2 target gene expression vs histological scores",
    x = "Histological score",
    y = "Gene expression",
    color = "Gene"
  )











