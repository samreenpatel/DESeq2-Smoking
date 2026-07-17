###############################################################
# Project : Differential Gene Expression Analysis using DESeq2
# Dataset : GSE220113
# Analysis: Smoking vs Non-Smoking
###############################################################

#==============================================================
# 1. Clear Workspace
#==============================================================

rm(list = ls())

#==============================================================
# 2. Load Libraries
#==============================================================

library(DESeq2)

#==============================================================
# 3. Load Prepared Data
#==============================================================

# Loads:
# counts
# metadata

load("prepared_data.RData")

# Check loaded objects

dim(counts)
dim(metadata)

head(metadata)

#==============================================================
# 4. Create DESeq2 Dataset
#==============================================================

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ Smoking
)

# View DESeq2 object

dds

# Number of genes and samples

dim(dds)

# View experimental design

design(dds)

# Sample information

colData(dds)

###############################################################
# 5. Filter Low Expression Genes
###############################################################

# Remove genes with fewer than 10 total reads
# across all samples

dds <- dds[
  rowSums(counts(dds)) >= 10,
]

# Remaining genes

dim(dds)

###############################################################
# 6. Run Differential Expression Analysis
###############################################################

dds <- DESeq(dds)

###############################################################
# 7. Available Comparisons
###############################################################

resultsNames(dds)

###############################################################
# 8. Extract Differential Expression Results
###############################################################

res <- results(dds)

head(res)

summary(res)

###############################################################
# 9. Order Genes by Adjusted P-value
###############################################################

resOrdered <- res[
  order(res$padj),
]

head(resOrdered)

###############################################################
# 10. Save Results
###############################################################

write.csv(
  as.data.frame(resOrdered),
  "DEG_results.csv"
)

###############################################################
# 11. Significant Genes
###############################################################

sig_genes <- subset(
  res,
  padj < 0.05
)

nrow(sig_genes)

head(sig_genes)

###############################################################
# 12. Save Significant Genes
###############################################################

write.csv(
  as.data.frame(sig_genes),
  "Significant_Genes.csv"
)