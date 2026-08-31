from pathlib import Path
import numpy as np


# ---------------------------------------------------------
# File paths
# ---------------------------------------------------------

project_dir = Path.cwd()

count_file = (
    project_dir
    / "data"
    / "raw"
    / "GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt"
)

selected_cells_file = (
    project_dir
    / "data"
    / "selected_cell_ids.txt"
)

output_dir = project_dir / "data" / "processed"
output_dir.mkdir(parents=True, exist_ok=True)

matrix_file = output_dir / "matrix.mtx"
features_file = output_dir / "features.tsv"
barcodes_file = output_dir / "barcodes.tsv"


# ---------------------------------------------------------
# Check that the required files exist
# ---------------------------------------------------------

if not count_file.exists():
    raise FileNotFoundError(f"Count file not found: {count_file}")

if not selected_cells_file.exists():
    raise FileNotFoundError(
        f"Selected-cell file not found: {selected_cells_file}"
    )


# ---------------------------------------------------------
# Read the selected cell names
# ---------------------------------------------------------

with selected_cells_file.open("r") as handle:
    requested_cells = {
        line.strip()
        for line in handle
        if line.strip()
    }

print(f"Requested cells: {len(requested_cells):,}")


# ---------------------------------------------------------
# Read the header and locate the selected cells
# ---------------------------------------------------------

with count_file.open("r") as input_handle:
    header = input_handle.readline().rstrip("\r\n").split("\t")

cell_names = header[1:]

selected_indices = []
selected_cell_names = []

for index, cell_name in enumerate(cell_names):
    if cell_name in requested_cells:
        selected_indices.append(index)
        selected_cell_names.append(cell_name)

selected_indices = np.asarray(selected_indices, dtype=np.int64)

missing_cells = requested_cells.difference(selected_cell_names)

print(f"Cells found in count matrix: {len(selected_cell_names):,}")
print(f"Cells missing from count matrix: {len(missing_cells):,}")

if missing_cells:
    example_missing = sorted(missing_cells)[:10]
    raise ValueError(
        "Some selected cells were not found. "
        f"Examples: {example_missing}"
    )


# ---------------------------------------------------------
# Write cell names in the same order as the matrix columns
# ---------------------------------------------------------

with barcodes_file.open("w") as output_handle:
    for cell_name in selected_cell_names:
        output_handle.write(f"{cell_name}\n")


# ---------------------------------------------------------
# Stream through the large matrix
# ---------------------------------------------------------

gene_count = 0
nonzero_count = 0

with count_file.open("r") as input_handle, \
        features_file.open("w") as feature_handle, \
        matrix_file.open("w+b") as matrix_handle:

    # Skip the header, which has already been examined.
    input_handle.readline()

    matrix_handle.write(
        b"%%MatrixMarket matrix coordinate integer general\n"
    )
    matrix_handle.write(
        b"% Created from GSE132465 raw UMI count matrix\n"
    )

    # Matrix Market requires dimensions before the values.
    # We reserve a fixed-width line and fill it after processing.
    dimension_position = matrix_handle.tell()
    matrix_handle.write(b" " * 99 + b"\n")

    for line in input_handle:
        line = line.rstrip("\r\n")

        if not line:
            continue

        first_tab = line.find("\t")

        if first_tab == -1:
            raise ValueError(
                f"Invalid row encountered after gene {gene_count}"
            )

        gene_name = line[:first_tab]

        all_values = np.fromstring(
            line[first_tab + 1:],
            dtype=np.int32,
            sep="\t"
        )

        if len(all_values) != len(cell_names):
            raise ValueError(
                f"Unexpected number of values for gene {gene_name}: "
                f"{len(all_values)} instead of {len(cell_names)}"
            )

        gene_count += 1
        feature_handle.write(f"{gene_name}\n")

        selected_values = all_values[selected_indices]
        nonzero_columns = np.flatnonzero(selected_values)

        if nonzero_columns.size > 0:
            entries = "".join(
                f"{gene_count} {column + 1} "
                f"{int(selected_values[column])}\n"
                for column in nonzero_columns
            )

            matrix_handle.write(entries.encode("ascii"))
            nonzero_count += nonzero_columns.size

        if gene_count % 500 == 0:
            print(
                f"Processed {gene_count:,} genes; "
                f"non-zero entries: {nonzero_count:,}"
            )

    dimension_line = (
        f"{gene_count} "
        f"{len(selected_cell_names)} "
        f"{nonzero_count}"
    ).encode("ascii")

    if len(dimension_line) > 99:
        raise ValueError("Matrix dimension line is unexpectedly long.")

    matrix_handle.seek(dimension_position)
    matrix_handle.write(
        dimension_line.ljust(99, b" ") + b"\n"
    )


print("")
print("Conversion completed successfully.")
print(f"Genes: {gene_count:,}")
print(f"Cells: {len(selected_cell_names):,}")
print(f"Non-zero entries: {nonzero_count:,}")
print(f"Output directory: {output_dir}")
