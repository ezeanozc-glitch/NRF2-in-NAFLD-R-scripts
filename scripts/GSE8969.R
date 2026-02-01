
# microarray with nrf2 KO in liver regeneration
library(affy)
library(limma)
library(GEOquery)
library(dplyr)
library(annotate)
library(mouse4302.db)
library(biomaRt)


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

fit <- lmFit(exprs_mat, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)


results <- topTable(
  fit2,
  coef = "KO_vs_WT",
  adjust.method = "BH",
  number = Inf
)


results_annot <- results %>%
  mutate(ProbeID = rownames(.))

results_annot <- results_annot %>%
  mutate(
    GeneSymbol = getSYMBOL(ProbeID, "mouse4302.db")
  )

results_annot_noNA <- results_annot[!is.na(results_annot$GeneSymbol), ]

results_annot_noNA$GeneSymbol <- toupper(results_annot_noNA$GeneSymbol)


#NFE2L2
nfe2l2_results <- results_annot_noNA[results_annot_noNA$GeneSymbol == "NFE2L2", ]
nfe2l2_results

downregulated_genes_adj <- results_annot_noNA[
  results_annot_noNA$logFC < -1 & results_annot_noNA$adj.P.Val < 0.05, 
]

nrow(downregulated_genes_adj)

downregulated_genes <- results_annot_noNA[
  results_annot_noNA$logFC < -1 & results_annot_noNA$P.Value < 0.05, 
]

nrow(downregulated_genes)

head(downregulated_genes)

downregulated_genes$GeneSymbol

as.vector(downregulated_genes$GeneSymbol)

unique(as.vector(downregulated_genes$GeneSymbol))



# check individual genes
extract_nrf2_results <- function(res_df, nrf2_genes) {
  # Subset rows matching nrf2_genes
  res_sub <- res_df[res_df$GeneSymbol %in% nrf2_genes, ]
  
  # Collapse multiple probes: keep probe with largest absolute logFC
  res_sub <- res_sub %>%
    dplyr::group_by(GeneSymbol) %>%
    dplyr::slice_max(order_by = abs(logFC), n = 1) %>%
    dplyr::ungroup()
  
  # Order by logFC
  res_sub <- res_sub[order(res_sub$logFC, decreasing = TRUE), ]
  
  return(res_sub)
}

nrf2_results <- extract_nrf2_results(results_annot_noNA, nrf2_genes)

nrf2_results_ordered <- nrf2_results[match(nrf2_genes, nrf2_results$GeneSymbol), ]

nrf2_results_ordered <- nrf2_results_ordered %>%
  dplyr::select(GeneSymbol, dplyr::everything())


write.csv(nrf2_results_ordered,
          file = "GSE8969_NRF2_KO_basal.csv",
          row.names = FALSE)

















