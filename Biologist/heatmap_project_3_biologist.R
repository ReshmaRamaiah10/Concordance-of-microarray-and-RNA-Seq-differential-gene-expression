#load required libraries
library(pheatmap)

#read filtered gene list
limma <- read.csv("DNADamage_deseq_norm_counts1.csv")

# Plotting heatmap (this works)
lim_hm <-pheatmap(as.matrix(limma[-1]), color = colorRampPalette(c("red", "yellow", "green"))(100),
            cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = FALSE, show_colnames = TRUE,
            labels_row = limma$Genes, scale = 'row', fontsize_row = 1.8)

