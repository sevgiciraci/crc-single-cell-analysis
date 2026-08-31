library(dplyr)
library(ggplot2)
library(msigdbr)
library(fgsea)


de_results <- read.csv(
  "results/epithelial_tumor_vs_normal_edger.csv",
  stringsAsFactors = FALSE
)

dim(de_results)

names(de_results)


rank_table <- de_results %>%
  filter(
    !is.na(gene),
    !is.na(logFC),
    !is.na(F)
  ) %>%
  group_by(gene) %>%
  slice_max(
    order_by = abs(F),
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()

gene_ranks <- sign(rank_table$logFC) *
  sqrt(rank_table$F)

names(gene_ranks) <- rank_table$gene

gene_ranks <- sort(
  gene_ranks,
  decreasing = TRUE
)

head(gene_ranks)
tail(gene_ranks)
length(gene_ranks)


hallmark_data <- msigdbr(
  db_species = "HS",
  species = "Homo sapiens",
  collection = "H"
)

dim(hallmark_data)

hallmark_pathways <- split(
  hallmark_data$gene_symbol,
  hallmark_data$gs_name
)

length(hallmark_pathways)

set.seed(42)

fgsea_results <- fgseaMultilevel(
  pathways = hallmark_pathways,
  stats = gene_ranks,
  minSize = 15,
  maxSize = 500,
  eps = 0
)

fgsea_results <- fgsea_results %>%
  as.data.frame() %>%
  arrange(padj)

fgsea_results %>%
  select(
    pathway,
    NES,
    pval,
    padj,
    size
  ) %>%
  head(20) %>%
  tibble::as_tibble() %>%
  print(n = 20)

fgsea_export <- fgsea_results %>%
  mutate(
    leadingEdge = vapply(
      leadingEdge,
      paste,
      collapse = ";",
      FUN.VALUE = character(1)
    )
  )

write.csv(
  fgsea_export,
  file = "results/epithelial_hallmark_fgsea.csv",
  row.names = FALSE
)

file.exists(
  "results/epithelial_hallmark_fgsea.csv"
)

fgsea_results %>%
  select(pathway, NES, pval, padj, size) %>%
  head(20) %>%
  tibble::as_tibble() %>%
  print(n = 20)



pathway_plot_data <- fgsea_results %>%
  filter(padj < 0.05) %>%
  mutate(
    pathway_label = gsub(
      "^HALLMARK_",
      "",
      pathway
    ),
    pathway_label = gsub(
      "_",
      " ",
      pathway_label
    ),
    direction = ifelse(
      NES > 0,
      "Tumor enriched",
      "Normal enriched"
    )
  )

pathway_plot_data %>%
  select(
    pathway_label,
    NES,
    padj,
    direction
  ) %>%
  tibble::as_tibble() %>%
  print(n = Inf)



pathway_plot <- ggplot(
  pathway_plot_data,
  aes(
    x = reorder(pathway_label, NES),
    y = NES,
    fill = direction
  )
) +
  geom_col(
    width = 0.75
  ) +
  coord_flip() +
  geom_hline(
    yintercept = 0,
    color = "black",
    linewidth = 0.4
  ) +
  scale_fill_manual(
    values = c(
      "Tumor enriched" = "#00A9B7",
      "Normal enriched" = "#F8766D"
    )
  ) +
  labs(
    title = "Hallmark pathways in epithelial cells",
    subtitle = "Paired tumor versus normal pseudobulk analysis",
    x = NULL,
    y = "Normalized enrichment score (NES)",
    fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "top"
  )

pathway_plot



ggsave(
  filename = "figures/epithelial_hallmark_pathways.png",
  plot = pathway_plot,
  width = 9,
  height = 7,
  dpi = 300
)



