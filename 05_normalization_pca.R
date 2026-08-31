library(Seurat)
library(ggplot2)
library(patchwork)

crc <- readRDS(
  "data/processed/crc_with_qc_metrics.rds"
)

crc
table(crc$qc_pass)


crc <- NormalizeData(
  crc,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = TRUE
)

Layers(crc[["RNA"]])

crc <- FindVariableFeatures(
  crc,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = TRUE
)


length(VariableFeatures(crc))
head(VariableFeatures(crc), 20)

variable_plot <- VariableFeaturePlot(crc)

top10_genes <- head(
  VariableFeatures(crc),
  10
)

variable_plot <- LabelPoints(
  plot = variable_plot,
  points = top10_genes,
  repel = TRUE
)

variable_plot


ggsave(
  filename = "figures/highly_variable_genes.png",
  plot = variable_plot,
  width = 7,
  height = 5,
  dpi = 300
)

table(crc$qc_pass)
Layers(crc[["RNA"]])
length(VariableFeatures(crc))


crc <- ScaleData(
  crc,
  features = VariableFeatures(crc),
  verbose = TRUE
)


crc <- RunPCA(
  crc,
  features = VariableFeatures(crc),
  npcs = 50,
  verbose = TRUE
)

crc[["pca"]]


elbow_plot <- ElbowPlot(
  crc,
  ndims = 50
)

elbow_plot


ggsave(
  filename = "figures/pca_elbow_plot.png",
  plot = elbow_plot,
  width = 7,
  height = 5,
  dpi = 300
)


pca_sd <- Stdev(crc, reduction = "pca")

pca_variance <- pca_sd^2

pca_percent <- (
  pca_variance /
    sum(pca_variance)
) * 100

pca_summary <- data.frame(
  PC = seq_along(pca_percent),
  Percent_variance = pca_percent,
  Cumulative_variance = cumsum(pca_percent)
)

head(pca_summary, 20)



write.csv(
  pca_summary,
  file = "results_pca_variance.csv",
  row.names = FALSE
)


dir.create(
  "results",
  showWarnings = FALSE
)

write.csv(
  pca_summary,
  file = "results/pca_variance_explained.csv",
  row.names = FALSE
)


pca_class_plot <- DimPlot(
  crc,
  reduction = "pca",
  group.by = "Class",
  pt.size = 0.2
)

pca_class_plot

ggsave(
  filename = "figures/pca_by_class.png",
  plot = pca_class_plot,
  width = 7,
  height = 5,
  dpi = 300
)


pca_patient_plot <- DimPlot(
  crc,
  reduction = "pca",
  group.by = "Patient",
  pt.size = 0.2
)

pca_patient_plot


ggsave(
  filename = "figures/pca_by_patient.png",
  plot = pca_patient_plot,
  width = 8,
  height = 6,
  dpi = 300
)

loading_plot <- VizDimLoadings(
  crc,
  dims = 1:2,
  reduction = "pca",
  nfeatures = 15
)

loading_plot

ggsave(
  filename = "figures/pca_loadings_pc1_pc2.png",
  plot = loading_plot,
  width = 11,
  height = 5,
  dpi = 300
)


saveRDS(
  crc,
  file = "data/processed/crc_normalized_pca.rds",
  compress = FALSE
)

file.exists("data/processed/crc_normalized_pca.rds")


head(pca_summary, 20)












