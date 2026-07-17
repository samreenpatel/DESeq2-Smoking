#Gene Annotation 

#==============================================================
# 1. Load annotation package
#==============================================================

library(org.Hs.eg.db)
library(AnnotationDbi)

#==============================================================
# 2. Convert DESeq2 results to data frame
#==============================================================

res_df <- as.data.frame(resOrdered)

# Add Entrez Gene IDs as a column
res_df$ENTREZID <- rownames(res_df)

# View first rows
head(res_df)

#==============================================================
# 3. Retrieve Gene Symbols and Gene Names
#==============================================================

gene_annotation <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = res_df$ENTREZID,
  keytype = "ENTREZID",
  columns = c(
    "SYMBOL",
    "GENENAME"
  )
)

# View annotation
head(gene_annotation)

#==============================================================
# 4. Merge annotation with DESeq2 results
#==============================================================

annotated_results <- merge(
  res_df,
  gene_annotation,
  by = "ENTREZID",
  all.x = TRUE
)

#==============================================================
# 5. Order by adjusted p-value
#==============================================================

annotated_results <- annotated_results[
  order(annotated_results$padj),
]

#==============================================================
# 6. View annotated results
#==============================================================

head(annotated_results)

#==============================================================
# 7. Save annotated results
#==============================================================

write.csv(
  annotated_results,
  "DEG_results_Annotated.csv",
  row.names = FALSE
)

###############################################################
# 8. Annotated Significant Genes
###############################################################

sig_annotated <- subset(
  annotated_results,
  padj < 0.05
)

# Number of significant genes
nrow(sig_annotated)

# View top genes
head(sig_annotated)

#==============================================================
# 9. Save Significant Annotated Genes
#==============================================================

write.csv(
  sig_annotated,
  "Significant_Genes_Annotated.csv",
  row.names = FALSE
)

###############################################################
# 10. Top 20 Differentially Expressed Genes
###############################################################

head(
  sig_annotated[, c(
    "ENTREZID",
    "SYMBOL",
    "GENENAME",
    "log2FoldChange",
    "padj"
  )],
  20
)
