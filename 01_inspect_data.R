count_file <- "data/raw/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt"

file.info(count_file)$size / 1024^3

connection <- file(count_file, open = "r")

header_line <- readLines(connection, n = 1)
first_data_line <- readLines(connection, n = 1)

close(connection)

nchar(header_line)
nchar(first_data_line)

substr(header_line, 1, 300)
substr(first_data_line, 1, 300)

rm(header_line, first_data_line, connection)
gc()

