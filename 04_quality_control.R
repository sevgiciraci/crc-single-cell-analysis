library(Seurat)
library(ggplot2)
library(patchwork)

crc <- readRDS(
  "data/processed/crc_matched_raw_seurat.rds"
)

crc

dir.create(
  "figures",
  showWarnings = FALSE
)

crc[["percent.mt"]] <- PercentageFeatureSet(
  crc,
  pattern = "^MT-"
)

summary(crc$nFeature_RNA)
summary(crc$nCount_RNA)
summary(crc$percent.mt)

quantile(
  crc$nFeature_RNA,
  probs = c(0, 0.01, 0.05, 0.50, 0.95, 0.99, 1)
)

quantile(
  crc$nCount_RNA,
  probs = c(0, 0.01, 0.05, 0.50, 0.95, 0.99, 1)
)

quantile(
  crc$percent.mt,
  probs = c(0, 0.01, 0.05, 0.50, 0.95, 0.99, 1)
)

qc_violin <- VlnPlot(
  crc,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  group.by = "Class",
  pt.size = 0,
  ncol = 3,
  layer = "counts"
)

qc_violin

ggsave(
  filename = "figures/qc_violin_by_class.png",
  plot = qc_violin,
  width = 12,
  height = 4,
  dpi = 300
)

qc_scatter <- FeatureScatter(
  crc,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA",
  group.by = "Class"
)

qc_scatter


ggsave(
  filename = "figures/qc_counts_vs_features.png",
  plot = qc_scatter,
  width = 7,
  height = 5,
  dpi = 300
)

saveRDS(
  crc,
  file = "data/processed/crc_with_qc_metrics.rds",
  compress = FALSE
)

file.exists("data/processed/crc_with_qc_metrics.rds")



summary(crc$nFeature_RNA)
summary(crc$nCount_RNA)
summary(crc$percent.mt)



quantile(
  crc$nFeature_RNA,
  probs = c(0, 0.01, 0.05, 0.50, 0.95, 0.99, 1)
)

quantile(
  crc$nCount_RNA,
  probs = c(0, 0.01, 0.05, 0.50, 0.95, 0.99, 1)
)

quantile(
  crc$percent.mt,
  probs = c(0, 0.01, 0.05, 0.50, 0.95, 0.99, 1)
)

crc$qc_pass <- (
  crc$nFeature_RNA > 200 &
    crc$nFeature_RNA < 6000 &
    crc$nCount_RNA > 1000 &
    crc$percent.mt < 20
)

table(crc$qc_pass)

saveRDS(
  crc,
  file = "data/processed/crc_with_qc_metrics.rds",
  compress = FALSE
)