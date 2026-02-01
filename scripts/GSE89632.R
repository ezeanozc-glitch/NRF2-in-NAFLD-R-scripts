library(tidyverse)
library(GEOquery)
library(limma)
library(illuminaHumanv4.db)
library(AnnotationDbi) 

# load metadata
gse <- getGEO("GSE89632", GSEMatrix = TRUE)
expr <- exprs(gse[[1]])
pheno <- pData(gse[[1]])



head(expr)
head(pheno)

# remove symbols
colnames(pheno)[colnames(pheno) == "steatosis (%):ch1"] <- "steatosis"
colnames(pheno)[colnames(pheno) == "ballooning (intensity):ch1"] <- "ballooning"
colnames(pheno)[colnames(pheno) == "fibrosis (stage):ch1"] <- "fibrosis_stage"
colnames(pheno)[colnames(pheno) == "lobular inflammation (severity):ch1"] <- "inflammation"
colnames(pheno)[colnames(pheno) == "nafld activity score:ch1"] <- "NAS"
colnames(pheno)[colnames(pheno) == "diagnosis:ch1"] <- "diagnosis"
colnames(pheno)[colnames(pheno) == "NAS:ch1"] <- "NAS"

# Frequency tables for each histology variable
table(pheno$steatosis)
table(pheno$ballooning)
table(pheno$fibrosis_stage)
table(pheno$inflammation)
table(pheno$NAS)
table(pheno$diagnosis)






# Keep only matching columns in expression matrix
expr_subset <- expr[, rownames(pheno)]
# Make sure factors are correct AFTER subsetting
pheno_subset$diagnosis <- factor(pheno_subset$diagnosis)  # must be HC, SS, NASH
# make design matrix
design <- model.matrix(~0 + diagnosis, data = pheno)
colnames(design) <- gsub("diagnosis", "", colnames(design))


contrast_matrix <- makeContrasts(
  NASH_vs_HC = NASH - HC,
  SS_vs_HC   = SS - HC,
  NASH_vs_SS = NASH - SS,
  levels = design
)

fit  <- lmFit(expr_subset, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

deg_NASH_HC <- topTable(fit2, coef="NASH_vs_HC", number=Inf, adjust.method="fdr")
deg_SS_HC   <- topTable(fit2, coef="SS_vs_HC",   number=Inf, adjust.method="fdr")
deg_NASH_SS <- topTable(fit2, coef="NASH_vs_SS", number=Inf, adjust.method="fdr")

deg_NASH_HC$GeneSymbol <- mapIds(
  illuminaHumanv4.db,
  keys = rownames(deg_NASH_HC),
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

deg_SS_HC$GeneSymbol <- mapIds(
  illuminaHumanv4.db,
  keys = rownames(deg_SS_HC),
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

deg_NASH_SS$GeneSymbol <- mapIds(
  illuminaHumanv4.db,
  keys = rownames(deg_NASH_SS),
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)


# create function to subset nrf2 target genes
extract_nrf2_results <- function(res_df, nrf2_genes) {
  res_sub <- res_df[res_df$GeneSymbol %in% nrf2_genes, ]
  res_sub <- res_sub %>%
    dplyr::group_by(GeneSymbol) %>%
    dplyr::slice_max(order_by = abs(logFC), n = 1) %>%
    dplyr::ungroup()
  res_sub <- res_sub[order(res_sub$logFC, decreasing = TRUE), ]
  return(res_sub)
}

# Run
nrf2_nafl <- extract_nrf2_results(deg_SS_HC,   nrf2_genes)
nrf2_nash <- extract_nrf2_results(deg_NASH_HC, nrf2_genes)

# Re-order to match original gene list
nrf2_nafl_ordered <- nrf2_nafl[match(nrf2_genes, nrf2_nafl$GeneSymbol), ]
nrf2_nash_ordered <- nrf2_nash[match(nrf2_genes, nrf2_nash$GeneSymbol), ]

# optional obese
write.csv(nrf2_nafl_ordered,
          file = "GSE89632_cov_obese_nrf2_nafld.csv",
          row.names = FALSE)

# export
write.csv(nrf2_nash_ordered,
          file = "GSE89632_cov_obese_nrf2_nash.csv",
          row.names = FALSE)









#------- correlation
histology_vars <- c(
  "steatosis",
  "ballooning",
  "fibrosis_stage",
  "inflammation",
  "NAS"
)

pheno[histology_vars] <- lapply(pheno[histology_vars], as.numeric)

probe2symbol <- mapIds(
  illuminaHumanv4.db,
  keys = rownames(expr),
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

expr_df <- as.data.frame(expr)
expr_df$SYMBOL <- probe2symbol

# Remove probes without symbols
expr_df <- expr_df[!is.na(expr_df$SYMBOL), ]

expr_nrf2 <- expr_df[expr_df$SYMBOL %in% nrf2_genes, ]

expr_nrf2_gene <- expr_nrf2 %>%
  group_by(SYMBOL) %>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))

expr_gene_mat <- as.matrix(expr_nrf2_gene[, -1])
rownames(expr_gene_mat) <- expr_nrf2_gene$SYMBOL

# Transpose: rows = samples, cols = genes
expr_gene_mat_t <- t(expr_gene_mat)

# Ensure alignment
all(rownames(expr_gene_mat_t) == rownames(pheno))

pheno <- pheno[rownames(expr_gene_mat_t), ]

clinical_vars <- c(
  "steatosis",
  "ballooning",
  "fibrosis_stage",
  "inflammation",
  "NAS"
)


cor_results <- data.frame()

for (gene in colnames(expr_gene_mat_t)) {
  for (var in clinical_vars) {
    
    x <- expr_gene_mat_t[, gene]
    y <- pheno[[var]]
    
    valid <- complete.cases(x, y)
    
    if (sum(valid) >= 10) {   # safeguard
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
          Spearman_rho = unname(cor_test$estimate),
          P_value = cor_test$p.value,
          N = sum(valid)
        )
      )
    }
  }
}


cor_results$Adj_P_value <- p.adjust(
  cor_results$P_value,
  method = "BH"
)

head(cor_results)

write.csv(
  cor_results,
  file = "GSE89632_NRF2_gene_histology_correlations.csv",
  row.names = FALSE
)


cor_results








#------- crate a graph for single genes
genes_to_plot <- c("NFE2L2")
nrf2_genes <- c("NFE2L2")

plot_df <- expr_gene_mat_t[, genes_to_plot, drop = FALSE] %>%
  as.data.frame() %>%
  mutate(Sample = rownames(.)) %>%
  pivot_longer(
    cols = all_of(genes_to_plot),
    names_to = "Gene",
    values_to = "Expression"
  ) %>%
  left_join(
    pheno %>%
      mutate(Sample = rownames(pheno)) %>%
      dplyr::select(
        Sample,
        steatosis,
        ballooning,
        fibrosis_stage,
        inflammation,
        NAS
      ),
    by = "Sample"
  ) %>%
  pivot_longer(
    cols = c(steatosis, ballooning, fibrosis_stage, inflammation, NAS),
    names_to = "Clinical_Score",
    values_to = "Score_Value"
  )

ggplot(plot_df, aes(x = Score_Value, y = Expression, colour = Gene)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ Clinical_Score, scales = "free_x") +
  theme_bw() +
  labs(
    title = "NRF2 gene expression vs histological severity (GSE89632)",
    x = "Histological score",
    y = "Gene expression",
    colour = "Gene"
  )










#---------- create a grpah fro multiple genes
genes_to_plot <- c(
  "GCLC", "GCLM", "GSR", "SLC7A11"
)

genes_to_plot <- c("NFE2L2")
nrf2_genes <- c("NFE2L2")


plot_df <- expr_gene_mat_t[, genes_to_plot] %>%
  as.data.frame() %>%
  mutate(Sample = rownames(.)) %>%
  pivot_longer(
    cols = all_of(genes_to_plot),
    names_to = "Gene",
    values_to = "Expression"
  ) %>%
  left_join(
    pheno %>%
      mutate(Sample = rownames(pheno)) %>%
      dplyr::select(Sample, steatosis, ballooning, fibrosis_stage, inflammation, NAS),
    by = "Sample"
  ) %>%
  pivot_longer(
    cols = c(steatosis, ballooning, fibrosis_stage, inflammation, NAS),
    names_to = "Clinical_Score",
    values_to = "Score_Value"
  )


ggplot(plot_df, aes(x = Score_Value, y = Expression, colour = Gene)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ Clinical_Score, scales = "free_x") +
  theme_bw() +
  labs(
    title = "NRF2 gene expression vs histological severity (GSE89632)",
    x = "Histological score",
    y = "Gene expression",
    colour = "Gene"
  )

