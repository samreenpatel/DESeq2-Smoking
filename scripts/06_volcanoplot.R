#Volcano Plot

#==============================================================
# 1. Clear Workspace
#==============================================================

rm(list = ls())

#==============================================================
# 2. Load Libraries
#==============================================================

library(DESeq2)
library(EnhancedVolcano)
library(org.Hs.eg.db)
library(AnnotationDbi)

#==============================================================
# 3. Load DESeq2 Results
#==============================================================

load("DESeq2_results.RData")

# Check objects
ls()

#==============================================================
# 4. Convert Results to Data Frame
#==============================================================

res_df <- as.data.frame(res)

# Add Entrez IDs
res_df$ENTREZID <- rownames(res_df)

#==============================================================
# 5. Retrieve Gene Symbols
#==============================================================

gene_annotation <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = res_df$ENTREZID,
  keytype = "ENTREZID",
  columns = "SYMBOL"
)

#==============================================================
# 6. Merge Annotation
#==============================================================

annotated_results <- merge(
  res_df,
  gene_annotation,
  by = "ENTREZID",
  all.x = TRUE
)

#==============================================================
# 7. Reorder to Match DESeq2 Results
#==============================================================

annotated_results <- annotated_results[
  match(rownames(res), annotated_results$ENTREZID),
]

#==============================================================
# 8. Replace Missing Symbols
#==============================================================

annotated_results$SYMBOL[
  is.na(annotated_results$SYMBOL)
] <- annotated_results$ENTREZID[
  is.na(annotated_results$SYMBOL)
]

#==============================================================
# 9. Volcano Plot
#==============================================================

EnhancedVolcano(
  
  res,
  
  lab = annotated_results$SYMBOL,
  
  x = "log2FoldChange",
  
  y = "padj",
  
  title = "Smoking vs Non-Smoking",
  
  subtitle = "GSE220113",
  
  caption = "Differential Expression Analysis using DESeq2",
  
  pCutoff = 0.05,
  
  FCcutoff = 1,
  
  pointSize = 3,
  
  labSize = 4,
  
  col = c(
    "grey70",
    "forestgreen",
    "royalblue",
    "red2"
  ),
  
  colAlpha = 1,
  
  drawConnectors = TRUE,
  
  widthConnectors = 0.5,
  
  legendPosition = "right"
)

#==============================================================
# 10. Save High Resolution Figure
#==============================================================

png(
  "Volcano_Plot.png",
  width = 1800,
  height = 1600,
  res = 300
)

EnhancedVolcano(
  
  res,
  
  lab = annotated_results$SYMBOL,
  
  x = "log2FoldChange",
  
  y = "padj",
  
  title = "Smoking vs Non-Smoking",
  
  subtitle = "GSE220113",
  
  caption = "Differential Expression Analysis using DESeq2",
  
  pCutoff = 0.05,
  
  FCcutoff = 1,
  
  pointSize = 3,
  
  labSize = 4,
  
  col = c(
    "grey70",
    "forestgreen",
    "royalblue",
    "red2"
  ),
  
  colAlpha = 1,
  
  drawConnectors = TRUE,
  
  widthConnectors = 0.5,
  
  legendPosition = "right"
)

dev.off()

select_genes <- c(
  "GPR15",
  "HOMER3",
  "GPR55",
  "ARHGAP6",
  "DNASE1L3",
  "RARA",
  "PIK3R6",
  "ITGAL",
  "CD226",
  "HFE"
)

EnhancedVolcano(
  res,
  lab = annotated_results$SYMBOL,
  selectLab = select_genes,
  
  x = "log2FoldChange",
  y = "padj",
  
  title = "Differential Gene Expression",
  subtitle = "Smoking vs Non-Smoking",
  caption = "Dataset: GSE220113",
  
  pCutoff = 0.05,
  FCcutoff = 1,
  
  pointSize = 3,
  labSize = 5,
  
  drawConnectors = TRUE,
  
  col = c(
    "grey70",
    "forestgreen",
    "royalblue",
    "red2"
  ),
  
  legendPosition = "right"
)
