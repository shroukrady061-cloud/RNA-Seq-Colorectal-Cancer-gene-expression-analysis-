ggsave("figures/PCA_Plot.png", p_pca, width = 8, height = 6, dpi = 300)
pcaData <- plotPCA(
  vsd,
  intgroup = c("cell", "time"),
  returnData = TRUE
)

percentVar <- round(100 * attr(pcaData, "percentVar"))

p_pca <- ggplot(
  pcaData,
  aes(
    PC1,
    PC2,
    color = cell,
    shape = factor(time)
  )
) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA Plot") +
  theme_minimal()

p_pca
ggsave(
  "figures/PCA_Plot.png",
  p_pca,
  width = 8,
  height = 6,
  dpi = 300
)
p_volcano <- ggplot(
  res_df,
  aes(
    x = log2FoldChange,
    y = -log10(padj)
  )
) +
  geom_point(aes(color = significant), alpha = 0.6) +
  scale_color_manual(
    values = c(
      "Upregulated" = "red",
      "Downregulated" = "blue",
      "Not significant" = "grey"
    )
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  labs(
    title = "Volcano Plot",
    x = "log2 Fold Change",
    y = "-log10 Adjusted p-value"
  ) +
  theme_minimal()

p_volcano

ggsave(
  "figures/Volcano_Plot.png",
  p_volcano,
  width = 8,
  height = 6,
  dpi = 300
)
p_go_up <- dotplot(
  ego_up,
  showCategory = 15,
  title = "GO Biological Process - Upregulated Genes"
)

p_go_up
ggsave(
  "figures/GO_Upregulated.png",
  p_go_up,
  width = 8,
  height = 6,
  dpi = 300
)
p_go_down <- dotplot(
  ego_down,
  showCategory = 15,
  title = "GO Biological Process - Downregulated Genes"
)

p_go_down

ggsave(
  "figures/GO_Downregulated.png",
  p_go_down,
  width = 8,
  height = 6,
  dpi = 300
)
p_kegg_up <- dotplot(
  ekegg_up,
  showCategory = 15,
  title = "KEGG Pathways - Upregulated Genes"
)

p_kegg_up

ggsave(
  "figures/KEGG_Upregulated.png",
  p_kegg_up,
  width = 8,
  height = 6,
  dpi = 300
)
p_kegg_down <- dotplot(
  ekegg_down,
  showCategory = 15,
  title = "KEGG Pathways - Downregulated Genes"
)

p_kegg_down

ggsave(
  "figures/KEGG_Downregulated.png",
  p_kegg_down,
  width = 8,
  height = 6,
  dpi = 300
)
