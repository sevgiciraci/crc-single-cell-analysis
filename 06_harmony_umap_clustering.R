library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)

crc <- readRDS(
  "data/processed/crc_normalized_pca.rds"
)

crc <- RunHarmony(
  object = crc,
  group.by.vars = "Patient",
  reduction.use = "pca",
  dims.use = 1:30,
  reduction.save = "harmony",
  verbose = TRUE
)

Reductions(crc)
dim(Embeddings(crc, reduction = "harmony"))

crc <- FindNeighbors(
  crc,
  reduction = "harmony",
  dims = 1:30,
  verbose = TRUE
)

crc <- FindClusters(
  crc,
  resolution = 0.4,
  random.seed = 42,
  verbose = TRUE
)

table(crc$seurat_clusters)

table(Idents(crc))

crc <- RunUMAP(
  crc,
  reduction = "harmony",
  dims = 1:30,
  n.neighbors = 30,
  min.dist = 0.3,
  seed.use = 42,
  verbose = TRUE
)

umap_clusters <- DimPlot(
  crc,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.2
) +
  NoLegend()

umap_clusters




