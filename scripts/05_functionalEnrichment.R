#Functional Enrichment Analysis

#==============================================================
# 1. Clear Workspace
#==============================================================

rm(list = ls())

#==============================================================
# 2. Load Libraries
#==============================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

#==============================================================
# 3. Load Annotated DESeq2 Results
#==============================================================

sig_genes <- read.csv(
  "Significant_Genes_Annotated.csv",
  stringsAsFactors = FALSE
)

# Inspect data

head(sig_genes)

dim(sig_genes)

#==============================================================
# 4. Extract Entrez Gene IDs
#==============================================================

gene_list <- sig_genes$ENTREZID

length(gene_list)

head(gene_list)

#==============================================================
# 5. Gene Ontology (Biological Process)
#==============================================================

go_bp <- enrichGO(
  gene          = gene_list,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# View results

head(go_bp)

#==============================================================
# 6. Save GO Biological Process Results
#==============================================================

write.csv(
  as.data.frame(go_bp),
  "GO_Biological_Process.csv",
  row.names = FALSE
)

#==============================================================
# 7. GO Dot Plot
#==============================================================

dotplot(
  go_bp,
  showCategory = 20,
  title = "GO Biological Processes"
)

#==============================================================
# 8. GO Bar Plot
#==============================================================

barplot(
  go_bp,
  showCategory = 20,
  title = "Top GO Biological Processes"
)

#==============================================================
# 9. Gene Ontology (Molecular Function)
#==============================================================

go_mf <- enrichGO(
  gene          = gene_list,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

write.csv(
  as.data.frame(go_mf),
  "GO_Molecular_Function.csv",
  row.names = FALSE
)

#==============================================================
# 10. Gene Ontology (Cellular Component)
#==============================================================

go_cc <- enrichGO(
  gene          = gene_list,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

write.csv(
  as.data.frame(go_cc),
  "GO_Cellular_Component.csv",
  row.names = FALSE
)

#==============================================================
# 11. KEGG Pathway Enrichment
#==============================================================

kegg <- enrichKEGG(
  gene         = gene_list,
  organism     = "hsa",
  pvalueCutoff = 0.05
)

# View results

head(kegg)

# Save results

write.csv(
  as.data.frame(kegg),
  "KEGG_Pathways.csv",
  row.names = FALSE
)

#==============================================================
# 12. KEGG Dot Plot
#==============================================================

dotplot(
  kegg,
  showCategory = 20,
  title = "KEGG Pathway Enrichment"
)

#==============================================================
# 13. KEGG Bar Plot
#==============================================================

barplot(
  kegg,
  showCategory = 20,
  title = "Top KEGG Pathways"
)

#==============================================================
# 14. Enrichment Map
#==============================================================
library(enrichplot)

go_bp <- pairwise_termsim(go_bp)


emapplot(go_bp)

#==============================================================
# 15. Gene-Concept Network
#==============================================================
cnetplot(
  go_bp,
  showCategory = 5
)
