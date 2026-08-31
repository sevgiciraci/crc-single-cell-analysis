library(Seurat)
library(dplyr)
library(ggplot2)

crc <- readRDS(
  "data/processed/crc_annotated.rds"
)


Idents(crc) <- "seurat_clusters"

levels(Idents(crc))

cluster_markers <- FindAllMarkers(
  object = crc,
  assay = "RNA",
  only.pos = TRUE,
  test.use = "wilcox",
  min.pct = 0.25,
  logfc.threshold = 0.25,
  max.cells.per.ident = 1000,
  random.seed = 42,
  verbose = TRUE
)


dim(cluster_markers)
head(cluster_markers)

write.csv(
  cluster_markers,
  file = "results/all_cluster_markers.csv",
  row.names = FALSE
)

top10_markers <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 10,
    with_ties = FALSE
  ) %>%
  ungroup()


write.csv(
  top10_markers,
  file = "results/top10_markers_per_cluster.csv",
  row.names = FALSE
)


write.csv(
  top10_markers,
  file = "results/top10_markers_per_cluster.csv",
  row.names = FALSE
)



top10_markers %>%
  select(
    cluster,
    gene,
    avg_log2FC,
    pct.1,
    pct.2,
    p_val_adj
  ) %>%
  head(30) %>%
  print(n = 30)
