library(dplyr)
library(GEOquery)
library(limma)
library(illuminaHumanv4.db)
library(AnnotationDbi)  # required for mapIds
library(illuminaHumanv4.db)
library(tidyr)
library(ggplot2)


gse <- getGEO("GSE89632", GSEMatrix = TRUE)
expr <- exprs(gse[[1]])
pheno <- pData(gse[[1]])



head(expr)
head(pheno)


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




#---------balloonign comparison
# Identify samples with non-NA ballooning
keep_samples <- !is.na(pheno$ballooning)

# Subset both phenotype and expression
pheno_clean <- pheno[keep_samples, ]
expr_clean <- expr[, keep_samples]

# Recreate design matrix
pheno_clean$ballooning <- factor(pheno_clean$ballooning, levels = c(0, 1, 2))
design <- model.matrix(~ ballooning, data = pheno_clean)

# Fit the model
fit <- lmFit(expr_clean, design)

# Define contrasts
contrast.matrix <- makeContrasts(
  "1_vs_0" = ballooning1,
  "2_vs_0" = ballooning2,
  "2_vs_1" = ballooning2 - ballooning1,
  levels = design
)

# Fit contrasts
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# For 1 vs 0
top_1_vs_0 <- topTable(fit2, coef = "1_vs_0", number = Inf, sort.by = "P")

# For 2 vs 0
top_2_vs_0 <- topTable(fit2, coef = "2_vs_0", number = Inf, sort.by = "P")

# For 2 vs 1
top_2_vs_1 <- topTable(fit2, coef = "2_vs_1", number = Inf, sort.by = "P")


# Function to annotate a topTable
annotate_results <- function(tt) {
  tt$GeneSymbol <- mapIds(
    illuminaHumanv4.db,
    keys = rownames(tt),
    column = "SYMBOL",
    keytype = "PROBEID",
    multiVals = "first"
  )
  
  tt$GeneName <- mapIds(
    illuminaHumanv4.db,
    keys = rownames(tt),
    column = "GENENAME",
    keytype = "PROBEID",
    multiVals = "first"
  )
  
  return(tt)
}

# Apply to each comparison
top_1_vs_0_annot <- annotate_results(top_1_vs_0)
top_2_vs_0_annot <- annotate_results(top_2_vs_0)
top_2_vs_1_annot <- annotate_results(top_2_vs_1)

# check individual genes
extract_nrf2_results <- function(res_df, nrf2_genes) {
  # Subset rows that match your nrf2_genes via the GeneSymbol column
  res_sub <- res_df[res_df$GeneSymbol %in% nrf2_genes, ]
  
  # Collapse multiple probes: keep the one with largest absolute logFC
  res_sub <- res_sub %>%
    dplyr::group_by(GeneSymbol) %>%
    dplyr::slice_max(order_by = abs(logFC), n = 1) %>%
    dplyr::ungroup()
  
  # Order by logFC
  res_sub <- res_sub[order(res_sub$logFC, decreasing = TRUE), ]
  
  return(res_sub)
}



nrf2_1_vs_0_annot  <- extract_nrf2_results(top_1_vs_0_annot, nrf2_genes)
nrf2_2_vs_0_annot   <- extract_nrf2_results(top_2_vs_0_annot,  nrf2_genes)

# Reorder nrf2_1_vs_0_annot to match nrf2_genes
nrf2_1_vs_0_ordered <- nrf2_1_vs_0_annot[match(nrf2_genes,
                                               nrf2_1_vs_0_annot$GeneSymbol), ]

# Reorder nrf2_2_vs_0_annot to match nrf2_genes
nrf2_2_vs_0_ordered <- nrf2_2_vs_0_annot[match(nrf2_genes,
                                               nrf2_2_vs_0_annot$GeneSymbol), ]

# For 1 vs 0
nrf2_1_vs_0_ordered <- nrf2_1_vs_0_ordered[
  , c("GeneSymbol", setdiff(colnames(nrf2_1_vs_0_ordered), "GeneSymbol"))
]

# For 2 vs 0
nrf2_2_vs_0_ordered <- nrf2_2_vs_0_ordered[
  , c("GeneSymbol", setdiff(colnames(nrf2_2_vs_0_ordered), "GeneSymbol"))
]


write.csv(nrf2_1_vs_0_ordered,
          file = "GSE89632_nrf2_1_vs_0_results.csv",
          row.names = FALSE)

write.csv(nrf2_2_vs_0_ordered,
          file = "GSE89632_nrf2_2_vs_0_results.csv",
          row.names = FALSE)







#-------- comparison for fibrosis stages
# Keep samples with non-missing fibrosis stage
keep_samples_fib <- !is.na(pheno$fibrosis_stage)

# Subset phenotype & expression
pheno_fib <- pheno[keep_samples_fib, ]
expr_fib <- expr[, keep_samples_fib]

# Make a collapsed fibrosis variable (keep original intact)
pheno_fib$fibrosis_collapsed <- as.character(pheno_fib$fibrosis_stage)
pheno_fib$fibrosis_collapsed[ pheno_fib$fibrosis_collapsed %in% c("2","3") ] <- "2_3"

# Convert to factor with the order you want (important for contrast names)
pheno_fib$fibrosis_collapsed <- factor(
  pheno_fib$fibrosis_collapsed,
  levels = c("0","1","2_3","4")
)

# Quick table to confirm group sizes
table(pheno_fib$fibrosis_collapsed)

# Build design matrix using the collapsed factor
design_fib_coll <- model.matrix(~ fibrosis_collapsed, data = pheno_fib)
colnames(design_fib_coll)  # inspect to see how R named the coefficients

# Fit linear model
fit_fib_coll <- lmFit(expr_fib, design_fib_coll)

# Define contrasts
contrast_fib_coll <- makeContrasts(
  "1_vs_0"     = fibrosis_collapsed1,
  "2_3_vs_0"   = fibrosis_collapsed2_3,
  "4_vs_0"     = fibrosis_collapsed4,
  "4_vs_1"     = fibrosis_collapsed4 - fibrosis_collapsed1,
  "4_vs_2_3"   = fibrosis_collapsed4 - fibrosis_collapsed2_3,
  levels = design_fib_coll
)

# Fit contrasts & empirical Bayes
fit_fib_coll2 <- contrasts.fit(fit_fib_coll, contrast_fib_coll)
fit_fib_coll2 <- eBayes(fit_fib_coll2)

# Extract topTables for comparisons
fib_1_vs_0   <- topTable(fit_fib_coll2, coef = "1_vs_0",   number = Inf, sort.by = "P")
fib_2_3_vs_0 <- topTable(fit_fib_coll2, coef = "2_3_vs_0", number = Inf, sort.by = "P")
fib_4_vs_0   <- topTable(fit_fib_coll2, coef = "4_vs_0",   number = Inf, sort.by = "P")
fib_4_vs_1   <- topTable(fit_fib_coll2, coef = "4_vs_1",   number = Inf, sort.by = "P")
fib_4_vs_2_3 <- topTable(fit_fib_coll2, coef = "4_vs_2_3", number = Inf, sort.by = "P")

# Annotate results using your annotation function
fib_1_vs_0_annot   <- annotate_results(fib_1_vs_0)
fib_2_3_vs_0_annot <- annotate_results(fib_2_3_vs_0)
fib_4_vs_0_annot   <- annotate_results(fib_4_vs_0)
fib_4_vs_1_annot   <- annotate_results(fib_4_vs_1)
fib_4_vs_2_3_annot <- annotate_results(fib_4_vs_2_3)

# 1 vs 0
nrf2_fib_1_vs_0 <- extract_nrf2_results(fib_1_vs_0_annot, nrf2_genes)
nrf2_fib_1_vs_0_ordered <- nrf2_fib_1_vs_0[
  match(nrf2_genes, nrf2_fib_1_vs_0$GeneSymbol),
]
nrf2_fib_1_vs_0_ordered <- nrf2_fib_1_vs_0_ordered[
  , c("GeneSymbol", setdiff(colnames(nrf2_fib_1_vs_0_ordered), "GeneSymbol"))
]

# 2_3 vs 0
nrf2_fib_2_3_vs_0 <- extract_nrf2_results(fib_2_3_vs_0_annot, nrf2_genes)
nrf2_fib_2_3_vs_0_ordered <- nrf2_fib_2_3_vs_0[
  match(nrf2_genes, nrf2_fib_2_3_vs_0$GeneSymbol),
]
nrf2_fib_2_3_vs_0_ordered <- nrf2_fib_2_3_vs_0_ordered[
  , c("GeneSymbol", setdiff(colnames(nrf2_fib_2_3_vs_0_ordered), "GeneSymbol"))
]

# 4 vs 0
nrf2_fib_4_vs_0 <- extract_nrf2_results(fib_4_vs_0_annot, nrf2_genes)
nrf2_fib_4_vs_0_ordered <- nrf2_fib_4_vs_0[
  match(nrf2_genes, nrf2_fib_4_vs_0$GeneSymbol),
]
nrf2_fib_4_vs_0_ordered <- nrf2_fib_4_vs_0_ordered[
  , c("GeneSymbol", setdiff(colnames(nrf2_fib_4_vs_0_ordered), "GeneSymbol"))
]

write.csv(nrf2_fib_1_vs_0_ordered,
          "GSE89632_nrf2_fibrosis_1_vs_0.csv",
          row.names = FALSE)

write.csv(nrf2_fib_2_3_vs_0_ordered,
          "GSE89632_nrf2_fibrosis_2_3_vs_0.csv",
          row.names = FALSE)

write.csv(nrf2_fib_4_vs_0_ordered,
          "GSE89632_nrf2_fibrosis_4_vs_0.csv",
          row.names = FALSE)





# -------- correlation with inflammation
# Keep only non-missing inflammation samples
keep_inf <- !is.na(pheno$inflammation)

pheno_inf <- pheno[keep_inf, ]
expr_inf  <- expr[, keep_inf]

# Collapse inflammation stages
pheno_inf$inflam_collapsed <- as.character(pheno_inf$inflammation)
pheno_inf$inflam_collapsed[ pheno_inf$inflam_collapsed %in% c("2","3") ] <- "2_3"

pheno_inf$inflam_collapsed <- factor(
  pheno_inf$inflam_collapsed,
  levels = c("0","1","2_3")
)

# Check group counts
table(pheno_inf$inflam_collapsed)

# Design matrix
design_inf <- model.matrix(~ inflam_collapsed, data = pheno_inf)
colnames(design_inf)

# Fit linear model
fit_inf <- lmFit(expr_inf, design_inf)


# Contrasts
contrast_inf <- makeContrasts(
  "1_vs_0"     = inflam_collapsed1,
  "2_3_vs_0"   = inflam_collapsed2_3,
  "2_3_vs_1"   = inflam_collapsed2_3 - inflam_collapsed1,
  levels = design_inf
)

# Fit contrasts
fit_inf2 <- contrasts.fit(fit_inf, contrast_inf)
fit_inf2 <- eBayes(fit_inf2)

# Extract tables
inf_1_vs_0   <- topTable(fit_inf2, coef = "1_vs_0", number = Inf, sort.by = "P")
inf_2_3_vs_0 <- topTable(fit_inf2, coef = "2_3_vs_0", number = Inf, sort.by = "P")
inf_2_3_vs_1 <- topTable(fit_inf2, coef = "2_3_vs_1", number = Inf, sort.by = "P")

# Annotate results
inf_1_vs_0_annot   <- annotate_results(inf_1_vs_0)
inf_2_3_vs_0_annot <- annotate_results(inf_2_3_vs_0)
inf_2_3_vs_1_annot <- annotate_results(inf_2_3_vs_1)

# Extract NRF2 gene results
nrf2_inf_1_vs_0   <- extract_nrf2_results(inf_1_vs_0_annot, nrf2_genes)
nrf2_inf_2_3_vs_0 <- extract_nrf2_results(inf_2_3_vs_0_annot, nrf2_genes)

# Reorder to match nrf2 list
order_to <- function(df, genes) {
  df2 <- df[match(genes, df$GeneSymbol), ]
  df2 <- df2[, c("GeneSymbol", setdiff(colnames(df2), "GeneSymbol"))]
  return(df2)
}

nrf2_inf_1_vs_0_ordered   <- order_to(nrf2_inf_1_vs_0,   nrf2_genes)
nrf2_inf_2_3_vs_0_ordered <- order_to(nrf2_inf_2_3_vs_0, nrf2_genes)

# Export CSVs
write.csv(nrf2_inf_1_vs_0_ordered,
          "GSE89632_nrf2_inflammation_1_vs_0.csv",
          row.names = FALSE)

write.csv(nrf2_inf_2_3_vs_0_ordered,
          "GSE89632_nrf2_inflammation_2_3_vs_0.csv",
          row.names = FALSE)





#-------------- steatosis compariosn as continuous varriable
# Make a clean numeric steatosis variable
pheno$steatosis_num <- as.numeric(as.character(pheno$steatosis))

# Check for NAs after conversion
summary(pheno$steatosis_num)

keep_samples_stea <- !is.na(pheno$steatosis_num)

expr_stea <- expr[, keep_samples_stea]
pheno_stea <- pheno[keep_samples_stea, ]

design_stea <- model.matrix(~ steatosis_num, data = pheno_stea)
colnames(design_stea)

fit_stea <- lmFit(expr_stea, design_stea)
fit_stea <- eBayes(fit_stea)

stea_results <- topTable(fit_stea, coef = "steatosis_num", number = Inf, sort.by = "P")

stea_results_annot <- annotate_results(stea_results)

nrf2_steatosis <- extract_nrf2_results(stea_results_annot, nrf2_genes)

nrf2_steatosis_ordered <- nrf2_steatosis[
  match(nrf2_genes, nrf2_steatosis$GeneSymbol),
]

as.data.frame(nrf2_steatosis_ordered)

nrf2_steatosis_ordered <- nrf2_steatosis_ordered[
  , c("GeneSymbol", setdiff(colnames(nrf2_steatosis_ordered), "GeneSymbol"))
]

write.csv(nrf2_steatosis_ordered,
          file = "nrf2_steatosis_results.csv",
          row.names = FALSE)





#---------- steatosis comparisons category
# Convert steatosis into numeric
pheno$steatosis_num <- as.numeric(as.character(pheno$steatosis))

# Categorize based on NASH CRN standard thresholds
pheno$steatosis_grade <- cut(
  pheno$steatosis_num,
  breaks = c(-Inf, 5, 33, 66, Inf),
  labels = c("0_Normal", "1_Mild", "2_Moderate", "3_Severe"),
  right = FALSE
)

# View counts per grade
table(pheno$steatosis_grade)

keep_stea <- !is.na(pheno$steatosis_grade)

pheno_stea <- pheno[keep_stea, ]
expr_stea  <- expr[, keep_stea]

pheno_stea$steatosis_grade <- factor(
  pheno_stea$steatosis_grade,
  levels = c("0_Normal", "1_Mild", "2_Moderate", "3_Severe")
)

design_stea_cat <- model.matrix(~ steatosis_grade, data = pheno_stea)
colnames(design_stea_cat)

fit_stea_cat <- lmFit(expr_stea, design_stea_cat)

contrast_stea_cat <- makeContrasts(
  Mild_vs_Normal      = steatosis_grade1_Mild,
  Moderate_vs_Normal  = steatosis_grade2_Moderate,
  Severe_vs_Normal    = steatosis_grade3_Severe,
  Moderate_vs_Mild    = steatosis_grade2_Moderate - steatosis_grade1_Mild,
  Severe_vs_Mild      = steatosis_grade3_Severe - steatosis_grade1_Mild,
  Severe_vs_Moderate  = steatosis_grade3_Severe - steatosis_grade2_Moderate,
  levels = design_stea_cat
)

fit_stea_cat2 <- contrasts.fit(fit_stea_cat, contrast_stea_cat)
fit_stea_cat2 <- eBayes(fit_stea_cat2)


stea_Mild_vs_Normal      <- topTable(fit_stea_cat2, coef="Mild_vs_Normal",     number=Inf)
stea_Moderate_vs_Normal  <- topTable(fit_stea_cat2, coef="Moderate_vs_Normal", number=Inf)
stea_Severe_vs_Normal    <- topTable(fit_stea_cat2, coef="Severe_vs_Normal",   number=Inf)

stea_Mild_vs_Normal_annot      <- annotate_results(stea_Mild_vs_Normal)
stea_Moderate_vs_Normal_annot  <- annotate_results(stea_Moderate_vs_Normal)
stea_Severe_vs_Normal_annot    <- annotate_results(stea_Severe_vs_Normal)

nrf2_Mild_vs_Normal     <- extract_nrf2_results(stea_Mild_vs_Normal_annot, nrf2_genes)
nrf2_Moderate_vs_Normal <- extract_nrf2_results(stea_Moderate_vs_Normal_annot, nrf2_genes)
nrf2_Severe_vs_Normal   <- extract_nrf2_results(stea_Severe_vs_Normal_annot, nrf2_genes)


order_to <- function(df, genes) {
  df2 <- df[match(genes, df$GeneSymbol), ]
  df2 <- df2[, c("GeneSymbol", setdiff(colnames(df2), "GeneSymbol"))]
  return(df2)
}

nrf2_Mild_vs_Normal_ordered <- order_to(nrf2_Mild_vs_Normal, nrf2_genes)
nrf2_Moderate_vs_Normal_ordered <- order_to(nrf2_Moderate_vs_Normal, nrf2_genes)
nrf2_Severe_vs_Normal_ordered <- order_to(nrf2_Severe_vs_Normal, nrf2_genes)


write.csv(
  nrf2_Mild_vs_Normal_ordered,
  "GSE89632_nrf2_steatosis_Mild_vs_Normal.csv",
  row.names = FALSE
)

write.csv(
  nrf2_Moderate_vs_Normal_ordered,
  "GSE89632_nrf2_steatosis_Moderate_vs_Normal.csv",
  row.names = FALSE
)

write.csv(
  nrf2_Severe_vs_Normal_ordered,
  "GSE89632_nrf2_steatosis_Severe_vs_Normal.csv",
  row.names = FALSE
)




















pheno_subset <- pheno

# Keep only matching columns in expression matrix
expr_subset <- expr[, rownames(pheno_subset)]
# Make sure factors are correct AFTER subsetting
pheno_subset$diagnosis <- factor(pheno_subset$diagnosis)  # must be HC, SS, NASH
# make design matrix
design <- model.matrix(~0 + diagnosis, data = pheno_subset)
colnames(design) <- gsub("diagnosis", "", colnames(design))
design

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

# NFRE2L2
deg_SS_HC[ !is.na(deg_SS_HC$GeneSymbol) & deg_SS_HC$GeneSymbol == "NFE2L2", ]
deg_NASH_HC[ !is.na(deg_NASH_HC$GeneSymbol) & deg_NASH_HC$GeneSymbol == "NFE2L2", ]


# SREBF check
deg_SS_HC[ !is.na(deg_SS_HC$GeneSymbol) & deg_SS_HC$GeneSymbol == "SREBF1", ]
deg_NASH_HC[ !is.na(deg_NASH_HC$GeneSymbol) & deg_NASH_HC$GeneSymbol == "SREBF1", ]

"HMOX1" %in% deg_NASH_HC$GeneSymbol


extract_nrf2_results <- function(res_df, nrf2_genes) {
  res_sub <- res_df[res_df$GeneSymbol %in% nrf2_genes, ]
  res_sub <- res_sub %>%
    dplyr::group_by(GeneSymbol) %>%
    dplyr::slice_max(order_by = abs(logFC), n = 1) %>%
    dplyr::ungroup()
  res_sub <- res_sub[order(res_sub$logFC, decreasing = TRUE), ]
  return(res_sub)
}

nrf2_nafl <- extract_nrf2_results(deg_SS_HC,   nrf2_genes)
nrf2_nash <- extract_nrf2_results(deg_NASH_HC, nrf2_genes)

# Re-order to match original gene list
nrf2_nafl_ordered <- nrf2_nafl[match(nrf2_genes, nrf2_nafl$GeneSymbol), ]
nrf2_nash_ordered <- nrf2_nash[match(nrf2_genes, nrf2_nash$GeneSymbol), ]

# optional obese
write.csv(nrf2_nafl_ordered,
          file = "GSE89632_cov_obese_nrf2_nafld.csv",
          row.names = FALSE)

write.csv(nrf2_nash_ordered,
          file = "GSE89632_cov_obese_nrf2_nash.csv",
          row.names = FALSE)









contrast_matrix <- makeContrasts(
  NASH_vs_HC = diagnosisNASH,
  SS_vs_HC   = diagnosisSS,
  NASH_vs_SS = diagnosisNASH - diagnosisSS,
  levels = design
)
contrast_matrix

fit <- lmFit(expr, design)         # expr is your expression matrix
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# ------- continue optional from here
deg_NASH_HC <- topTable(fit2, coef="NASH_vs_HC", number=Inf, adjust.method="fdr")
deg_SS_HC <- topTable(fit2, coef="SS_vs_HC", number=Inf, adjust.method="fdr")
deg_NASH_SS <- topTable(fit2, coef="NASH_vs_SS", number=Inf, adjust.method="fdr")

head(deg_NASH_HC)


deg_NASH_HC$GeneSymbol <- mapIds(
  illuminaHumanv4.db,         # annotation database
  keys = rownames(deg_NASH_HC), # the probe IDs
  column = "SYMBOL",           # what you want: gene symbol
  keytype = "PROBEID",         # the type of keys you have
  multiVals = "first"          # if multiple symbols exist, take the first
)

deg_SS_HC$GeneSymbol <- mapIds(
  illuminaHumanv4.db,
  keys = rownames(deg_SS_HC),
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)


# check individual genes
extract_nrf2_results <- function(res_df, nrf2_genes) {
  # Subset rows that match your nrf2_genes via the GeneSymbol column
  res_sub <- res_df[res_df$GeneSymbol %in% nrf2_genes, ]
  
  # Collapse multiple probes: keep the one with largest absolute logFC
  res_sub <- res_sub %>%
    dplyr::group_by(GeneSymbol) %>%
    dplyr::slice_max(order_by = abs(logFC), n = 1) %>%
    dplyr::ungroup()
  
  # Order by logFC
  res_sub <- res_sub[order(res_sub$logFC, decreasing = TRUE), ]
  
  return(res_sub)
}


nrf2_nafl  <- extract_nrf2_results(deg_SS_HC, nrf2_genes)
nrf2_nash   <- extract_nrf2_results(deg_NASH_HC,  nrf2_genes)


# Reorder nrf2_nafl to match the order in nrf2_genes
nrf2_nafl_ordered <- nrf2_nafl[match(nrf2_genes, nrf2_nafl$GeneSymbol), ]

# Reorder nrf2_nash similarly
nrf2_nash_ordered <- nrf2_nash[match(nrf2_genes, nrf2_nash$GeneSymbol), ]

nrf2_nafl_ordered <- nrf2_nafl_ordered %>%
  dplyr::select(GeneSymbol, dplyr::everything())

nrf2_nash_ordered <- nrf2_nash_ordered %>%
  dplyr::select(GeneSymbol, dplyr::everything())

write.csv(nrf2_nafl_ordered,
          file = "GSE89632_cov_nrf2_nafld.csv",
          row.names = FALSE)

write.csv(nrf2_nash_ordered,
          file = "GSE89632_cov_nrf2_nash.csv",
          row.names = FALSE)


all(rownames(pheno_subset) %in% colnames(expr))  # should be TRUE

j

















library(GEOquery)
library(limma)
library(illuminaio)
library(GEOquery)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)


folder_path <- "C:/Users/chukw/OneDrive/Desktop/capstone project/Rstudio capstone/capstone/GSE89632_RAW"
files <- list.files(path = folder_path, pattern = "\\.idat$", full.names = TRUE)
head(files)



# --- 3. Read all IDATs ---
idat_list <- lapply(files, readIDAT)



# --- 2. Extract raw expression matrix ---
expr_matrix <- do.call(cbind, lapply(idat_list, function(x) x$Quants[, "MeanBinData"]))
rownames(expr_matrix) <- idat_list[[1]]$Quants[, "CodesBinData"]
colnames(expr_matrix) <- gsub("_Grn.idat$", "", basename(files))


# load pheno data
gse <- getGEO("GSE89632", GSEMatrix = TRUE)
pheno <- pData(phenoData(gse[[1]]))
head(pheno)
colnames(pheno)[colnames(pheno) == "diagnosis:ch1"] <- "diagnosis"
colnames(pheno)[colnames(pheno) ==  "fibrosis (stage):ch1"] <- "fibrosis_stage"
colnames(pheno)[colnames(pheno) ==  "lobular inflammation (severity):ch1"] <- "lobular_inflammation"
colnames(pheno)[colnames(pheno) ==  "ballooning (intensity):ch1"] <- "ballooning"
colnames(pheno)[colnames(pheno) ==  "nafld activity score:ch1"] <- "nas_score"


unique(pheno$diagnosis)
unique(pheno$fibrosis)
unique(pheno$lobular_inflammation)
unique(pheno$ballooning)
unique(pheno$nas_score)

# --- 4. Match GSM IDs and filter phenotype at the same time ---
gsm_ids <- sub("_.*", "", colnames(expr_matrix))
pheno <- pheno[match(gsm_ids, pheno$geo_accession), ]
colnames(expr_matrix) <- gsm_ids

# --- 5. Log2 transform ---
expr_matrix_log2 <- log2(expr_matrix + 1)

# --- 6. Quantile normalization ---
norm_matrix <- normalizeBetweenArrays(expr_matrix_log2, method = "quantile")

# --- 7. Create design matrix ---
group <- pheno$diagnosis
design <- model.matrix(~ 0 + factor(group))
colnames(design) <- levels(factor(group))

# --- 8. Fit limma model ---
fit <- lmFit(norm_matrix, design)

# --- 9. Define contrasts ---
contrast.matrix <- makeContrasts(NASHvsSS = NASH - SS, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# All genes
all_results <- topTable(fit2, coef = "NASHvsSS", number = Inf, sort.by = "none")

# --- 9. Add probe IDs ---
all_results$ProbeID <- rownames(all_results)

# --- 10. Get GPL annotation ---
gpl <- getGEO("GPL14951", destdir = tempdir())
gpl_table <- Table(gpl)

# --- 11. Map numeric probe IDs to GPL IDs ---
probe_map <- data.frame(
  NumericID = idat_list[[1]]$Quants$IllumicodeBinData,
  GPL_ID = gpl_table$ID[match(idat_list[[1]]$Quants$IllumicodeBinData,
                              idat_list[[1]]$Quants$IllumicodeBinData)]
)

# Add GPL_ID to all_results
all_results$ProbeID <- probe_map$GPL_ID[match(rownames(all_results), probe_map$NumericID)]


# --- 12. Merge with GPL table to get gene symbols ---
all_results_annot <- merge(
  all_results,
  gpl_table[, c("ID", "ILMN_Gene")],
  by.x = "ProbeID",
  by.y = "ID",
  all.x = TRUE
)

# --- 13. Check Nrf2 (NFE2L2) ---
all_results_annot[all_results_annot$ILMN_Gene == "NFE2L2", ]


# comparison of HC vs NASH
# Define a new contrast for NASH vs HC
contrast.matrix <- makeContrasts(NASHvsHC = NASH - HC, levels = design)

# Fit contrasts and compute statistics
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Get results for all genes
all_results <- topTable(fit2, coef = "NASHvsHC", number = Inf, sort.by = "none")

# Filter for significant DEGs: adjusted p < 0.05 and |logFC| ≥ 1
DEGs <- all_results_annot %>%
  filter(adj.P.Val < 0.05, abs(logFC) >= 1)


# Remove rows without gene symbols
all_results_human <- all_results_annot %>%
  filter(!is.na(ILMN_Gene), ILMN_Gene != "")

# Create a named numeric vector: gene name → logFC
gene_ranks <- all_results_human$logFC
names(gene_ranks) <- all_results_human$ILMN_Gene

# Collapse duplicates by keeping the strongest |logFC|
gene_ranks <- tapply(gene_ranks, names(gene_ranks), function(x) x[which.max(abs(x))])

# Convert from list to numeric vector while preserving names
gene_ranks <- unlist(gene_ranks)
gene_ranks <- sort(gene_ranks, decreasing = TRUE)

nrf2_sets <- read.gmt("C:/Users/chukw/OneDrive/Desktop/capstone project/NRF2_pathway_from_KEG.gmt")


# Convert gene names to uppercase for consistency
names(gene_ranks) <- toupper(names(gene_ranks))
nrf2_sets$gene <- toupper(nrf2_sets$gene)

# Check overlap between ranked genes and NRF2 gene set
overlap_genes <- intersect(names(gene_ranks), nrf2_sets$gene)
cat("Number of overlapping genes:", length(overlap_genes), "\n")
print(head(overlap_genes))


# Run GSEA
gsea_result <- GSEA(
  geneList = gene_ranks,
  TERM2GENE = nrf2_sets,
  minGSSize = 5,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

# Get all NRF2 genes present in your data
nrf2_gene_data <- all_results_human %>%
  mutate(ILMN_Gene = toupper(ILMN_Gene)) %>%
  filter(ILMN_Gene %in% nrf2_sets$gene) %>%
  arrange(desc(logFC))

nrf2_gene_data



# correlation with fibrosis score
# --- 1. Map probes to GPL IDs for normalized matrix ---
probe_map <- data.frame(
  NumericID = idat_list[[1]]$Quants$IllumicodeBinData,
  GPL_ID = gpl_table$ID[match(idat_list[[1]]$Quants$IllumicodeBinData,
                              idat_list[[1]]$Quants$IllumicodeBinData)],
  stringsAsFactors = FALSE
)

rownames_before <- rownames(norm_matrix)
match_idx <- match(rownames_before, as.character(probe_map$NumericID))
mapped_probe_ids <- probe_map$GPL_ID[match_idx]

norm_mat_gpl <- norm_matrix
rownames(norm_mat_gpl) <- mapped_probe_ids
norm_mat_gpl <- norm_mat_gpl[!is.na(rownames(norm_mat_gpl)), , drop = FALSE]


# --- 2. Identify NFE2L2 probe(s) ---
nfe2l2_probes <- gpl_table$ID[toupper(gpl_table$ILMN_Gene) == "NFE2L2"]
probes_in_data <- nfe2l2_probes[nfe2l2_probes %in% rownames(norm_mat_gpl)]
cat("NFE2L2 probes present in data:", paste(probes_in_data, collapse = ", "), "\n")


# --- 3. Extract NFE2L2 expression across samples ---
nfe2l2_expr_by_sample <- as.numeric(norm_mat_gpl["ILMN_1790909", ])
names(nfe2l2_expr_by_sample) <- colnames(norm_mat_gpl)

common_samples <- intersect(names(nfe2l2_expr_by_sample), rownames(pheno))
df <- data.frame(
  sample = common_samples,
  nfe2l2_expr = nfe2l2_expr_by_sample[common_samples],
  fibrosis_stage = pheno$fibrosis_stage[match(common_samples, rownames(pheno))]
)

# Remove any NA values
df <- df[!is.na(df$nfe2l2_expr) & !is.na(df$fibrosis_stage), ]

# Convert fibrosis_stage to numeric
df$fibrosis_stage <- as.numeric(as.character(df$fibrosis_stage))


# correlation test
cor_spearman <- cor.test(df$nfe2l2_expr, df$fibrosis_stage, method = "spearman")
cat("Spearman rho =", round(cor_spearman$estimate, 3),
    "p =", signif(cor_spearman$p.value, 3), "\n")








# correlation with lobular inflammation
# --- Columns to correlate with NFE2L2 ---
histology_cols <- c("lobular_inflammation", "ballooning", "nas_score")

# --- Prepare data frame with NFE2L2 expression ---
common_samples <- intersect(names(nfe2l2_expr_by_sample), rownames(pheno))
nfe2l2_df <- data.frame(
  sample = common_samples,
  nfe2l2_expr = nfe2l2_expr_by_sample[common_samples]
)

# --- Loop over histology columns and compute correlations ---
cor_results <- data.frame(
  feature = character(),
  spearman_rho = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (col in histology_cols) {
  # Extract the histology values
  histo_vals <- pheno[[col]][match(common_samples, rownames(pheno))]
  # Remove NAs
  valid_idx <- !is.na(histo_vals) & !is.na(nfe2l2_df$nfe2l2_expr)
  # Convert to numeric if not already
  histo_vals <- as.numeric(as.character(histo_vals[valid_idx]))
  nfe2l2_vals <- nfe2l2_df$nfe2l2_expr[valid_idx]
  
  # Compute Spearman correlation
  cor_test <- cor.test(nfe2l2_vals, histo_vals, method = "spearman")
  
  # Store results
  cor_results <- rbind(cor_results, data.frame(
    feature = col,
    spearman_rho = round(cor_test$estimate, 3),
    p_value = signif(cor_test$p.value, 3)
  ))
}

# --- View results ---
print(cor_results)







jjk

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

genes_to_plot <- c(
  "GSTA1", "GSTA2", "GSTA3", "GSTA5",
  "GSTM1", "GSTM3"
)

genes_to_plot <- c(
  "G6PD", "PGD", "TKT", "TALDO1", "ME1", "IDH1"
)

genes_to_plot <- c(
  "NFE2L2"
)


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
