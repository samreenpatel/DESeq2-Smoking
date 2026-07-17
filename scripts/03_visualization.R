#Visualization

#==============================================================
# 1. Clear Workspace
#==============================================================

rm(list = ls())

#==============================================================
# 2. Load Libraries
#==============================================================

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(matrixStats)

#==============================================================
# 3. Load DESeq2 Objects
#==============================================================

load("DESeq2_results.RData")

#==============================================================
# 4. MA Plot
#==============================================================

png("MA_Plot.png", width = 1800, height = 1600, res = 300)

plotMA(
  res,
  ylim = c(-5,5),
  main = "MA Plot: Smoking vs Non-Smoking"
)

dev.off()

#==============================================================
# 5. Plot Counts for Top Significant Genes
#==============================================================

top6 <- rownames(sig_genes)[1:6]

png("Top6_Gene_Counts.png", width = 2200, height = 1800, res = 300)

par(mfrow = c(2,3))

for(gene in top6){
  
  plotCounts(
    dds,
    gene = gene,
    intgroup = "Smoking"
  )
  
}

mtext(
  "Top Differentially Expressed Genes",
  outer = TRUE,
  cex = 1.4
)

dev.off()

#==============================================================
# 6. Variance Stabilizing Transformation
#==============================================================

vsd <- vst(
  dds,
  blind = FALSE
)

#==============================================================
# 7. PCA Plot
#==============================================================

png("PCA_Plot.png", width = 1800, height = 1600, res = 300)

plotPCA(
  vsd,
  intgroup = "Smoking"
)

dev.off()

#==============================================================
# 8. Sample Distance Heatmap
#==============================================================

sampleDists <- dist(
  t(
    assay(vsd)
  )
)

sampleDistMatrix <- as.matrix(sampleDists)

png("Sample_Distance_Heatmap.png", width = 1800, height = 1800, res = 300)

pheatmap(
  sampleDistMatrix,
  display_numbers = FALSE,
  fontsize = 8,
  border_color = "grey",
  main = "Sample Distance Heatmap"
)

dev.off()

#==============================================================
# 9. Heatmap of Top 50 Significant Genes
#==============================================================

top50 <- head(
  order(res$padj),
  50
)

mat <- assay(vsd)[top50, ]

mat <- mat - rowMeans(mat)

png("Top50_Gene_Heatmap.png", width = 1800, height = 2200, res = 300)

pheatmap(
  mat,
  annotation_col = as.data.frame(colData(vsd)[,"Smoking",drop=FALSE]),
  show_rownames = FALSE,
  main = "Top 50 Differentially Expressed Genes"
)

dev.off()