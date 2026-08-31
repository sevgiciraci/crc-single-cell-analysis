annotation_file <- "data/raw/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt"

annotation <- read.delim(
  annotation_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE

)

dim(annotation)
names(annotation)
head(annotation)

object.size(annotation) / 1024^2

patient_class_table <- table(
  annotation$Patient,
  annotation$Class
)

patient_class_table

matched_patients <- rownames(patient_class_table)[
  patient_class_table[, "Normal"] > 0 &
    patient_class_table[, "Tumor"] > 0
]

matched_patients
length(matched_patients)

selected_annotation <- annotation[
  annotation$Patient %in% matched_patients,
]

dim(selected_annotation)
table(selected_annotation$Class)

table(
  selected_annotation$Cell_type,
  selected_annotation$Class
)


selected_cell_ids <- selected_annotation$Index

writeLines(
  selected_cell_ids,
  "data/selected_cell_ids.txt"
)

file.exists("data/selected_cell_ids.txt")

length(readLines("data/selected_cell_ids.txt"))




