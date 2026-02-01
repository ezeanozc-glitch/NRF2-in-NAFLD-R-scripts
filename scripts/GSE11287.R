
library(affy)
library(limma)
library(GEOquery)
library(mouse4302.db)
library(annotate)
library(dplyr)
library(biomaRt)

folder_path <- "C:\\Users\\chukw\\OneDrive\\Desktop\\capstone project\\Rstudio capstone\\capstone\\GSE11287_RAW"
list.files(folder_path)

data <- ReadAffy(celfile.path = folder_path)

sampleNames(data)

eset <- rma(data)

exprs_mat <- exprs(eset)
colnames(exprs_mat) <- sub("\\.CEL$", "", colnames(exprs_mat))


gse <- getGEO("GSE11287", GSEMatrix = TRUE)
pheno <- pData(phenoData(gse[[1]]))
head(pheno)

all(colnames(exprs_mat) == rownames(pheno))

# Extract the genotype from the characteristics column
genotype <- sub("genotype: ", "", sapply(strsplit(pheno$characteristics_ch1, ";"), `[`, 1))

group <- ifelse(grepl("knockout", genotype, ignore.case = TRUE),
                "KO", "WT")

pheno$group <- group

design <- model.matrix(~ 0 + group, data = pheno)
colnames(design) <- c("KO", "WT")

contrast <- makeContrasts(KOvsWT = KO - WT, levels = design)

fit <- lmFit(exprs_mat, design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

results <- topTable(fit2, number = Inf, adjust.method = "BH")

results_annot <- results %>%
  mutate(ProbeID = rownames(.))

results_annot <- results_annot %>%
  mutate(
    GeneSymbol = getSYMBOL(ProbeID, "mouse4302.db")
  )

results_basal_annot_noNA <- results_annot[!is.na(results_annot$GeneSymbol), ]

results_basal_annot_noNA$GeneSymbol <- toupper(results_basal_annot_noNA$GeneSymbol)

head(results_basal_annot_noNA )




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
  file = "GSE11287_NRF2_KO_basal.csv",
  row.names = FALSE
)













