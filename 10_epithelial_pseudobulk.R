library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

crc <- readRDS(
  "data/processed/crc_annotated.rds"
)

epithelial_cell_counts <- crc[[]] %>%
  filter(
    Broad_cell_type == "Epithelial cells"
  ) %>%
  count(
    Patient,
    Class,
    name = "n_cells"
  ) %>%
  pivot_wider(
    names_from = Class,
    values_from = n_cells,
    values_fill = 0
  ) %>%
  arrange(Patient)

epithelial_cell_counts


eligible_patients <- epithelial_cell_counts %>%
  filter(
    Normal >= 20,
    Tumor >= 20
  ) %>%
  pull(Patient)

eligible_patients
length(eligible_patients)


cells_to_keep <- colnames(crc)[
  crc$Broad_cell_type == "Epithelial cells" &
    crc$Patient %in% eligible_patients
]

length(cells_to_keep)


epithelial <- subset(
  crc,
  cells = cells_to_keep
)

epithelial
table(
  epithelial$Patient,
  epithelial$Class
)

epithelial$sample_id <- paste(
  epithelial$Patient,
  epithelial$Class,
  sep = "_"
)

table(epithelial$sample_id)


pseudobulk_result <- AggregateExpression(
  epithelial,
  assays = "RNA",
  group.by = "sample_id",
  return.seurat = FALSE,
  verbose = FALSE
)

pseudobulk_counts <- pseudobulk_result$RNA

dim(pseudobulk_counts)
colnames(pseudobulk_counts)


pseudobulk_metadata <- epithelial[[]] %>%
  distinct(
    sample_id,
    Patient,
    Class
  ) %>%
  mutate(
    sample_id = gsub(
      "_",
      "-",
      sample_id,
      fixed = TRUE
    )
  ) %>%
  as.data.frame()

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

anyNA(pseudobulk_metadata$Patient)

saveRDS(
  pseudobulk_counts,
  file = "data/processed/epithelial_pseudobulk_counts.rds"
)

write.csv(
  pseudobulk_metadata,
  file = "results/epithelial_pseudobulk_metadata.csv",
  row.names = FALSE
)


