library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

crc <- readRDS(
  "data/processed/crc_annotated.rds"
)

cell_metadata <- crc[[]]

head(cell_metadata)


sample_information <- cell_metadata %>%
  distinct(
    Sample,
    Patient,
    Class
  )

all_cell_types <- sort(
  unique(cell_metadata$Broad_cell_type)
)

cell_composition <- cell_metadata %>%
  count(
    Sample,
    Broad_cell_type,
    name = "n_cells"
  ) %>%
  complete(
    Sample,
    Broad_cell_type = all_cell_types,
    fill = list(n_cells = 0)
  ) %>%
  left_join(
    sample_information,
    by = "Sample"
  ) %>%
  group_by(Sample) %>%
  mutate(
    total_cells = sum(n_cells),
    proportion = n_cells / total_cells
  ) %>%
  ungroup()


head(cell_composition)

write.csv(
  cell_composition,
  file = "results/cell_type_composition_by_sample.csv",
  row.names = FALSE
)

composition_barplot <- ggplot(
  cell_composition,
  aes(
    x = Sample,
    y = proportion,
    fill = Broad_cell_type
  )
) +
  geom_col(width = 0.85) +
  scale_y_continuous(
    labels = percent_format()
  ) +
  scale_fill_brewer(
    palette = "Set2"
  ) +
  labs(
    x = "Sample",
    y = "Cell-type proportion",
    fill = "Cell type"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1
    )
  )

composition_barplot



ggsave(
  filename = "figures/cell_composition_stacked_barplot.png",
  plot = composition_barplot,
  width = 12,
  height = 6,
  dpi = 300
)


cell_composition$Class <- factor(
  cell_composition$Class,
  levels = c("Normal", "Tumor")
)

composition_paired_plot <- ggplot(
  cell_composition,
  aes(
    x = Class,
    y = proportion,
    group = Patient
  )
) +
  geom_line(
    alpha = 0.5,
    colour = "grey50"
  ) +
  geom_point(
    aes(colour = Class),
    size = 2
  ) +
  facet_wrap(
    ~ Broad_cell_type,
    scales = "free_y",
    ncol = 3
  ) +
  scale_y_continuous(
    labels = percent_format()
  ) +
  scale_colour_manual(
    values = c(
      "Normal" = "#F8766D",
      "Tumor" = "#00BFC4"
    )
  ) +
  labs(
    x = NULL,
    y = "Cell-type proportion"
  ) +
  theme_classic() +
  theme(
    legend.position = "none"
  )

composition_paired_plot


ggsave(
  filename = "figures/cell_composition_paired.png",
  plot = composition_paired_plot,
  width = 10,
  height = 9,
  dpi = 300
)



composition_wide <- cell_composition %>%
  select(
    Patient,
    Class,
    Broad_cell_type,
    proportion
  ) %>%
  pivot_wider(
    names_from = Class,
    values_from = proportion
  )



composition_statistics <- composition_wide %>%
  group_by(Broad_cell_type) %>%
  summarise(
    median_normal = median(Normal),
    median_tumor = median(Tumor),
    median_paired_difference = median(Tumor - Normal),
    p_value = wilcox.test(
      Tumor,
      Normal,
      paired = TRUE,
      exact = FALSE
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    FDR = p.adjust(
      p_value,
      method = "BH"
    )
  ) %>%
  arrange(FDR)

composition_statistics


write.csv(
  composition_statistics,
  file = "results/cell_composition_statistics.csv",
  row.names = FALSE
)


