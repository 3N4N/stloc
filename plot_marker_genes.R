require(ggforce)
require(patchwork)
require(ggpubr)

counts_raw <- read.delim("./datasets/skin_cancer/GSE144239_ST_P2_S1_counts.tsv", header = TRUE, row.names = 1)

coords_raw <- do.call(rbind, strsplit(rownames(counts_raw), "x"))
coords <- apply(coords_raw, 1:2, as.numeric)
colnames(coords) <- c("x","y")
rownames(coords) <- rownames(counts_raw)
counts <- t(counts_raw)


clusterData <- read.delim("./datasets/skin_cancer/reference_markers_for_NMF.tsv", header = TRUE)
clusterGenes <- clusterData[,8]
clusterGenes <- unique(clusterGenes)
commonGenes <- intersect(rownames(counts), clusterGenes)
clusterData <- as.data.frame(clusterData)

if (!file.exists("output/plot_expr")) {
  system("mkdir output/plot_expr")
}

plotGeneExpr = function(df_res) {
  plot_expr <- ggplot(df_res, aes(x = x, y = -y)) +
    geom_point(aes(colour = expr), size = 3) +
    theme_minimal() +
    theme(panel.grid = element_blank()) +
    theme(axis.text = element_blank()) +
    xlab("") +
    ylab("") +
    labs(colour = "") +
    theme(legend.position = "bottom") +
    theme(plot.title = element_text(hjust = 0.5, face = "italic")) +
    scale_color_viridis_c() +
    coord_fixed() +
    guides(colour = guide_colourbar(title.position = "top",
                                    title.hjust = 0.5)) +
    theme(legend.key.width = unit(0.5, "inches")) +
    theme(plot.title = element_text(size = 20)) +
    theme(axis.title = element_text(size = 15)) +
    theme(legend.title = element_text(size = 15)) +
    labs(colour = "Gene Expression") +
    NULL

  expr_leg = as_ggplot(get_legend(plot_expr))
  scater::multiplot(plot_expr + theme(legend.position = "none")
                    + theme(plot.margin = margin(10,0,-10,0)),
                    expr_leg,
                    layout = matrix(
                      c(1,1,1,1,1,1,
                        1,1,1,1,1,1,
                        1,1,1,1,1,1,
                        3,3,2,2,3,3), ncol = 6, byrow = TRUE))

}

for (gene in commonGenes) {
  print(gene)
  df_res <- data.frame(x = coords[,"x"],
                       y = coords[,"y"],
                       expr = counts[gene,])
  # print(df_res)
  pdf(paste0("output/plot_expr/", gene, ".pdf"),
      height = 6, width = 10, onefile = F)
  plotGeneExpr(df_res)
  dev.off()

  # break
}
