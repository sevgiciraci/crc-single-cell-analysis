library(Seurat)
library(Matrix)

matrix_file <- "data/processed/matrix.mtx"
features_file <- "data/processed/features.tsv"
barcodes_file <- "data/processed/barcodes.tsv"
annotation_file <- "data/raw/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt"

file.exists(
  c(
    matrix_file,
    features_file,
    barcodes_file,
    annotation_file
  )
)

features <- readLines(features_file)
barcodes <- readLines(barcodes_file)

length(features)
length(barcodes)

counts_triplet <- readMM(matrix_file)

dim(counts_triplet)

class(counts_triplet)

sum(duplicated(features))

features <- make.unique(features)

counts <- as(counts_triplet, "CsparseMatrix")

class(counts)

rm(counts_triplet)
gc()


rownames(counts) <- features
colnames(counts) <- barcodes

dim(counts)
head(rownames(counts))
head(colnames(counts))


annotation <- read.delim(
  annotation_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

metadata <- annotation[
  match(barcodes, annotation$Index),
  ,
  drop = FALSE
]


sum(is.na(metadata$Index))
identical(metadata$Index, barcodes)

sum(is.na(metadata$Index))
identical(metadata$Index, barcodes)

rownames(metadata) <- metadata$Index
metadata$Index <- NULL

names(metadata)[names(metadata) == "Cell_type"] <- "Author_cell_type"
names(metadata)[names(metadata) == "Cell_subtype"] <- "Author_cell_subtype"

head(metadata)

crc <- CreateSeuratObject(
  counts = counts,
  project = "CRC_GSE132465",
  meta.data = metadata,
  min.cells = 3,
  min.features = 200
)



