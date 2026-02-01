
library(affy)
library(limma)
library(GEOquery)
library(mouse4302.db)
library(annotate)
library(dplyr)
library(biomaRt)


folder_path <- "C:\\Users\\chukw\\OneDrive\\Desktop\\capstone project\\Rstudio capstone\\capstone\\GSE15633_RAW"
list.files(folder_path)

# Read CEL files from folder
data <- ReadAffy(celfile.path = folder_path)

# Check sample names
sampleNames(data)

# RMA preprocessing
eset <- rma(data)

# Extract expression matrix
exprs_mat <- exprs(eset)
colnames(exprs_mat) <- sub("\\.CEL$", "", colnames(exprs_mat))

gse <- getGEO("GSE15633", GSEMatrix = TRUE)
pheno <- pData(phenoData(gse[[1]]))
head(pheno)


all(colnames(exprs_mat) == rownames(pheno))


table(pheno$description)

# create a group
pheno$Group <- NA
pheno$Group[pheno$description == "Basal genetic control wild-type gene expression"] <- "WT"
pheno$Group[pheno$description == "Basal conditional Keap1 knockout gene expression"] <- "KO"

table(pheno$Group)

# keep only WT and KO samples
keep_samples <- which(pheno$Group %in% c("WT", "KO"))
pheno_sub <- pheno[keep_samples, ]
exprs_sub <- exprs_mat[, keep_samples]

table(pheno_sub$Group)

# create design matrix
pheno_sub$Group <- factor(pheno_sub$Group, levels = c("WT", "KO"))

design <- model.matrix(~ Group, data = pheno_sub)

# fit model
fit <- lmFit(exprs_sub, design)
fit <- eBayes(fit)
results <- topTable(fit, coef = "GroupKO", number = Inf)


results_annot <- results %>%
  mutate(ProbeID = rownames(.))

results_annot <- results_annot %>%
  mutate(
    GeneSymbol = getSYMBOL(ProbeID, "mouse4302.db")
  )

results_basal_annot_noNA <- results_annot[!is.na(results_annot$GeneSymbol), ]

results_basal_annot_noNA$GeneSymbol <- toupper(results_basal_annot_noNA$GeneSymbol)

head(results_basal_annot_noNA )




# CDOO comparison
# Step 1: Create a new group column for CDDO-Im treated samples
pheno$CDDO_Group <- NA
pheno$CDDO_Group[pheno$description == "genetic control wild-type mice treated with CDDO-Im gene expression"] <- "WT"
pheno$CDDO_Group[pheno$description == "conditional Keap1 knockout mice treated with CDDO-Im gene expression"] <- "KO"

# Step 3: Keep only the treated WT and KO samples
keep_cdood <- which(pheno$CDDO_Group %in% c("WT", "KO"))
pheno_cdood <- pheno[keep_cdood, ]
exprs_cdood <- exprs_mat[, keep_cdood]

# Step 3: Keep only the treated WT and KO samples
keep_cdood <- which(pheno$CDDO_Group %in% c("WT", "KO"))
pheno_cdood <- pheno[keep_cdood, ]
exprs_cdood <- exprs_mat[, keep_cdood]

# Step 4: Make sure the group is a factor with WT as baseline
pheno_cdood$CDDO_Group <- factor(pheno_cdood$CDDO_Group, levels = c("WT", "KO"))

# Step 5: Create design matrix
design_cdood <- model.matrix(~ CDDO_Group, data = pheno_cdood)


# Step 6: Fit linear model and compute statistics
fit_cdood <- lmFit(exprs_cdood, design_cdood)
fit_cdood <- eBayes(fit_cdood)
results_cdood <- topTable(fit_cdood, coef = "CDDO_GroupKO", number = Inf)


results_annot <- results_cdood %>%
  mutate(ProbeID = rownames(.))

results_annot <- results_annot %>%
  mutate(
    GeneSymbol = getSYMBOL(ProbeID, "mouse4302.db")
  )

results_cddo_annot_noNA <- results_annot[!is.na(results_annot$GeneSymbol), ]

results_cddo_annot_noNA$GeneSymbol <- toupper(results_cddo_annot_noNA$GeneSymbol)

head(results_cddo_annot_noNA)



# comparison of KO + CDDO and WT
# create group
pheno$KOCDDO_vs_WT <- NA
pheno$KOCDDO_vs_WT[pheno$description == "conditional Keap1 knockout mice treated with CDDO-Im gene expression"] <- "KO_CDDO"
pheno$KOCDDO_vs_WT[pheno$description == "Basal genetic control wild-type gene expression"] <- "WT"

keep_compare <- which(pheno$KOCDDO_vs_WT %in% c("KO_CDDO", "WT"))
pheno_compare <- pheno[keep_compare, ]
exprs_compare <- exprs_mat[, keep_compare]

pheno_compare$KOCDDO_vs_WT <- factor(
  pheno_compare$KOCDDO_vs_WT,
  levels = c("WT", "KO_CDDO")  # baseline WT
)

design_compare <- model.matrix(~ KOCDDO_vs_WT, data = pheno_compare)

fit_compare <- lmFit(exprs_compare, design_compare)
fit_compare <- eBayes(fit_compare)

results_compare <- topTable(
  fit_compare,
  coef = "KOCDDO_vs_WTKO_CDDO",
  number = Inf
)


results_annot <- results_compare %>%
  mutate(ProbeID = rownames(.))

results_annot <- results_annot %>%
  mutate(
    GeneSymbol = getSYMBOL(ProbeID, "mouse4302.db")
  )

results_kocddo_annot_noNA  <- results_annot[!is.na(results_annot$GeneSymbol), ]

results_kocddo_annot_noNA $GeneSymbol <- toupper(results_kocddo_annot_noNA $GeneSymbol)

head(results_kocddo_annot_noNA )

# check individual genes
extract_nrf2_results <- function(res_df, nrf2_genes) {
  res_sub <- res_df[res_df$GeneSymbol %in% nrf2_genes, ]
  
  res_sub <- res_sub %>%
    dplyr::group_by(GeneSymbol) %>%
    dplyr::slice_max(order_by = abs(logFC), n = 1) %>%
    dplyr::ungroup()
  
  res_sub <- res_sub[order(res_sub$logFC, decreasing = TRUE), ]
  
  return(res_sub)
}


nrf2_basal <- extract_nrf2_results(
  results_basal_annot_noNA,
  nrf2_genes
)

nrf2_basal <- nrf2_basal[
  match(nrf2_genes, nrf2_basal$GeneSymbol),
]

nrf2_basal <- nrf2_basal %>%
  dplyr::select(GeneSymbol, dplyr::everything())

write.csv(
  nrf2_basal,
  file = "GSE15633_NRF2_KO_basal.csv",
  row.names = FALSE
)


nrf2_cddo <- extract_nrf2_results(
  results_cddo_annot_noNA,
  nrf2_genes
)

nrf2_cddo <- nrf2_cddo[
  match(nrf2_genes, nrf2_cddo$GeneSymbol),
]

nrf2_cddo <- nrf2_cddo %>%
  dplyr::select(GeneSymbol, dplyr::everything())

write.csv(
  nrf2_cddo,
  file = "GSE15633_NRF2_KO_CDDO.csv",
  row.names = FALSE
)


nrf2_kocddo <- extract_nrf2_results(
  results_kocddo_annot_noNA,
  nrf2_genes
)

nrf2_kocddo <- nrf2_kocddo[
  match(nrf2_genes, nrf2_kocddo$GeneSymbol),
]

nrf2_kocddo <- nrf2_kocddo %>%
  dplyr::select(GeneSymbol, dplyr::everything())

write.csv(
  nrf2_kocddo,
  file = "GSE15633_NRF2_KOCDDO_vs_WT.csv",
  row.names = FALSE
)
