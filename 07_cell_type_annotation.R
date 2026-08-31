library(Seurat)
library(dplyr)
library(ggplot2)

crc <- readRDS(
  "data/processed/crc_harmony_umap_clustered.rds"
)


umap_author_types <- DimPlot(
  crc,
  reduction = "umap",
  group.by = "Author_cell_type",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.2
)

umap_author_types


ggsave(
  filename = "figures/umap_author_cell_types.png",
  plot = umap_author_types,
  width = 9,
  height = 6,
  dpi = 300
)


cluster_author_table <- table(
  Cluster = crc$seurat_clusters,
  Author_type = crc$Author_cell_type
)

cluster_author_table

cluster_author_percent <- round(
  prop.table(
    cluster_author_table,
    margin = 1
  ) * 100,
  digits = 1
)

cluster_author_percent

write.csv(
  as.data.frame.matrix(cluster_author_percent),
  file = "results/cluster_author_celltype_percent.csv"
)


canonical_markers <- c(
  "EPCAM", "KRT8", "KRT18",
  "CD3D", "CD3E", "TRAC",
  "NKG7", "GNLY", "KLRD1",
  "MS4A1", "CD79A", "CD37",
  "LST1", "LYZ", "FCER1G",
  "TPSAB1", "TPSB2", "KIT",
  "COL1A1", "COL1A2", "DCN",
  "PECAM1", "VWF", "EMCN"
)

marker_dotplot <- DotPlot(
  crc,
  features = canonical_markers,
  group.by = "seurat_clusters"
) +
  RotatedAxis()

marker_dotplot

ggsave(
  filename = "figures/canonical_marker_dotplot.png",
  plot = marker_dotplot,
  width = 14,
  height = 7,
  dpi = 300
)

cluster_to_celltype <- c(
  "0"  = "T/NK cells",
  "1"  = "T/NK cells",
  "2"  = "Epithelial cells",
  "3"  = "Myeloid cells",
  "4"  = "B cells",
  "5"  = "B cells",
  "6"  = "Endothelial cells",
  "7"  = "B cells",
  "8"  = "Fibroblasts",
  "9"  = "Epithelial cells",
  "10" = "Fibroblasts",
  "11" = "Epithelial cells",
  "12" = "Fibroblasts",
  "13" = "Fibroblasts",
  "14" = "Myeloid cells",
  "15" = "Epithelial cells",
  "16" = "T/NK cells",
  "17" = "Fibroblasts",
  "18" = "Fibroblasts",
  "19" = "Mast cells",
  "20" = "Myeloid cells"
)

crc$Broad_cell_type <- unname(
  cluster_to_celltype[
    as.character(crc$seurat_clusters)
  ]
)

sum(is.na(crc$Broad_cell_type))

table(crc$Broad_cell_type)

umap_broad_types <- DimPlot(
  crc,
  reduction = "umap",
  group.by = "Broad_cell_type",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.2
)

umap_broad_types

ggsave(
  filename = "figures/umap_broad_cell_types.png",
  plot = umap_broad_types,
  width = 9,
  height = 6,
  dpi = 300
)


annotation_decisions <- data.frame(
  Cluster = names(cluster_to_celltype),
  Broad_cell_type = unname(cluster_to_celltype)
)

write.csv(
  annotation_decisions,
  file = "results/cluster_annotation_decisions.csv",
  row.names = FALSE
)


saveRDS(
  crc,
  file = "data/processed/crc_annotated.rds",
  compress = FALSE
)


file.exists("data/processed/crc_annotated.rds")
