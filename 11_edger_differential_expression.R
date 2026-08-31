library(edgeR)
library(dplyr)
library(ggplot2)

pseudobulk_counts <- readRDS(
  "data/processed/epithelial_pseudobulk_counts.rds"
)

pseudobulk_metadata <- read.csv(
  "results/epithelial_pseudobulk_metadata.csv",
  stringsAsFactors = FALSE
)

rownames(pseudobulk_metadata) <-
  pseudobulk_metadata$sample_id

pseudobulk_metadata <- pseudobulk_metadata[
  colnames(pseudobulk_counts),
  ,
  drop = FALSE
]

identical(
  rownames(pseudobulk_metadata),
  colnames(pseudobulk_counts)
)

pseudobulk_metadata$Patient <- factor(
  pseudobulk_metadata$Patient
)

pseudobulk_metadata$Class <- factor(
  pseudobulk_metadata$Class,
  levels = c("Normal", "Tumor")
)


dge <- DGEList(
  counts = pseudobulk_counts,
  samples = pseudobulk_metadata
)

design <- model.matrix(
  ~ Patient + Class,
  data = pseudobulk_metadata
)

colnames(design)

keep_genes <- filterByExpr(
  dge,
  design = design
)

table(keep_genes)

dge <- dge[
  keep_genes,
  ,
  keep.lib.sizes = FALSE
]

dim(dge)


dge <- calcNormFactors(
  dge,
  method = "TMM"
)


dge$samples[
  ,
  c(
    "sample_id",
    "Patient",
    "Class",
    "lib.size",
    "norm.factors"
  )
]


dge <- estimateDisp(
  dge,
  design,
  robust = TRUE
)

fit <- glmQLFit(
  dge,
  design,
  robust = TRUE
)

qlf <- glmQLFTest(
  fit,
  coef = "ClassTumor"
)


de_results <- topTags(
  qlf,
  n = Inf,
  sort.by = "PValue"
)$table

de_results$gene <- rownames(de_results)

de_results <- de_results %>%
  relocate(gene)

head(de_results, 20)


de_results <- de_results %>%
  mutate(
    significance = case_when(
      FDR < 0.05 & logFC >= 1 ~ "Up in tumor",
      FDR < 0.05 & logFC <= -1 ~ "Down in tumor",
      TRUE ~ "Not significant"
    )
  )

table(de_results$significance)


write.csv(
  de_results,
  file = "results/epithelial_tumor_vs_normal_edger.csv",
  row.names = FALSE
)

saveRDS(
  fit,
  file = "data/processed/epithelial_edger_fit.rds"
)

colnames(design)
table(keep_genes)
table(de_results$significance)


if (!requireNamespace("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel")
}

library(ggrepel)



mds_result <- plotMDS(
  dge,
  plot = FALSE
)

mds_data <- data.frame(
  Dimension_1 = mds_result$x,
  Dimension_2 = mds_result$y,
  pseudobulk_metadata
)

mds_plot <- ggplot(
  mds_data,
  aes(
    x = Dimension_1,
    y = Dimension_2
  )
) +
  geom_line(
    aes(group = Patient),
    colour = "grey70",
    linewidth = 0.6
  ) +
  geom_point(
    aes(colour = Class),
    size = 3
  ) +
  geom_text_repel(
    aes(label = sample_id),
    size = 3,
    max.overlaps = Inf
  ) +
  scale_colour_manual(
    values = c(
      "Normal" = "#F8766D",
      "Tumor" = "#00BFC4"
    )
  ) +
  theme_classic() +
  labs(
    title = "Epithelial pseudobulk expression",
    x = "Leading logFC dimension 1",
    y = "Leading logFC dimension 2"
  )

mds_plot


ggsave(
  filename = "figures/epithelial_pseudobulk_mds.png",
  plot = mds_plot,
  width = 8,
  height = 6,
  dpi = 300
)


volcano_labels <- bind_rows(
  de_results %>%
    filter(significance == "Up in tumor") %>%
    arrange(FDR) %>%
    slice_head(n = 5),
  de_results %>%
    filter(significance == "Down in tumor") %>%
    arrange(FDR) %>%
    slice_head(n = 5)
)


de_results <- de_results %>%
  mutate(
    minus_log10_FDR = -log10(
      pmax(FDR, 1e-300)
    )
  )


volcano_plot <- ggplot(
  de_results,
  aes(
    x = logFC,
    y = minus_log10_FDR,
    colour = significance
  )
) +
  geom_point(
    alpha = 0.6,
    size = 1.2
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed",
    colour = "grey40"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    colour = "grey40"
  ) +
  geom_text_repel(
    data = volcano_labels,
    aes(label = gene),
    size = 3,
    max.overlaps = Inf
  ) +
  scale_colour_manual(
    values = c(
      "Down in tumor" = "#377EB8",
      "Not significant" = "grey75",
      "Up in tumor" = "#E41A1C"
    )
  ) +
  theme_classic() +
  labs(
    title = "Tumor versus normal epithelial cells",
    subtitle = "Paired patient-level pseudobulk edgeR analysis",
    x = "log2 fold change",
    y = "-log10 FDR",
    colour = NULL
  )

volcano_plot



ggsave(
  filename = "figures/epithelial_volcano_plot.png",
  plot = volcano_plot,
  width = 8,
  height = 6,
  dpi = 300
)


if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}


top_heatmap_genes <- bind_rows(
  de_results %>%
    filter(significance == "Up in tumor") %>%
    arrange(FDR) %>%
    slice_head(n = 10),
  de_results %>%
    filter(significance == "Down in tumor") %>%
    arrange(FDR) %>%
    slice_head(n = 10)
) %>%
  pull(gene) %>%
  unique()


log_cpm <- cpm(
  dge,
  log = TRUE,
  prior.count = 2
)

heatmap_matrix <- log_cpm[
  top_heatmap_genes,
  ,
  drop = FALSE
]


heatmap_zscore <- t(
  scale(
    t(heatmap_matrix)
  )
)


heatmap_annotation <- data.frame(
  Patient = pseudobulk_metadata$Patient,
  Class = pseudobulk_metadata$Class
)

rownames(heatmap_annotation) <-
  rownames(pseudobulk_metadata)


pheatmap::pheatmap(
  heatmap_zscore,
  annotation_col = heatmap_annotation,
  cluster_cols = FALSE,
  border_color = NA,
  fontsize_row = 9,
  filename = "figures/epithelial_top_de_heatmap.png",
  width = 10,
  height = 8
)


file.exists(
  "figures/epithelial_top_de_heatmap.png"
)

pheatmap::pheatmap(
  heatmap_zscore,
  annotation_col = heatmap_annotation,
  cluster_cols = FALSE,
  border_color = NA,
  fontsize_row = 9,
  width = 10,
  height = 8
)



saveRDS(
  list(
    dge = dge,
    design = design,
    fit = fit,
    qlf = qlf,
    metadata = pseudobulk_metadata
  ),
  file = "data/processed/epithelial_edger_analysis.rds"
)

top20_results <- de_results %>%
  select(
    gene,
    logFC,
    logCPM,
    F,
    PValue,
    FDR,
    significance
  ) %>%
  head(20)

print(
  tibble::as_tibble(top20_results),
  n = 20
)






