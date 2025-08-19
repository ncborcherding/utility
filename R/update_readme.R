# R/update_readme.R

# This script generates a summary of the cohort and updates the README.md file.

library(Seurat)

# 1. Load the final integrated Seurat object
# This assumes the script is run from the root of the project directory
seurat_object_path <- "./data/integratedSeuratObject_scVI_clustered.rds"
if (!file.exists(seurat_object_path)) {
  stop("Final Seurat object not found. Please run the integration pipeline first.")
}
object.integrated <- readRDS(seurat_object_path)

# 2. Calculate summary statistics
n_cells <- ncol(object.integrated)
n_tissues <- length(unique(object.integrated$tissue))

# Infer patient IDs by splitting 'orig.ident'
patient_ids <- sapply(strsplit(object.integrated$orig.ident, "\\."), `[`, 1)
n_patients <- length(unique(patient_ids))

# Count the number of unique tumor samples
# Assuming 'tissue' metadata column contains 'Tumor' for tumor samples
tumor_samples <- unique(object.integrated$orig.ident[grepl("Tumor", object.integrated$tissue, ignore.case = TRUE)])
n_tumors <- length(tumor_samples)


# 3. Create the Markdown summary table
# Using kableExtra for a nicely formatted table, but a simple paste should also work.
# For simplicity and fewer dependencies, I'll use simple string pasting.
summary_table_md <- paste(
  "| Metric          | Count     |",
  "|-----------------|-----------|",
  paste("| Total Cells     |", format(n_cells, big.mark = ",")),
  paste("| Unique Patients |", n_patients),
  paste("| Unique Tissues  |", n_tissues),
  paste("| Tumor Samples   |", n_tumors),
  sep = "\n"
)

# Add a title for the new section
summary_section <- paste(
  "### Cohort Summary",
  "",
  "This table provides a high-level overview of the integrated dataset.",
  "",
  summary_table_md,
  "",
  sep = "\n"
)


# 4. Read the README.md and insert the table
readme_path <- "./README.md"
readme_content <- readLines(readme_path)

# Find the line where the detailed cohort information starts
insertion_point <- which(readme_content == "#### Cohort Information")

if (length(insertion_point) == 0) {
  stop("Could not find the '#### Cohort Information' section in README.md")
}

# Check if our summary section already exists to prevent duplicates
if (any(grepl("### Cohort Summary", readme_content))) {
    print("Cohort summary section already exists. Skipping update.")
} else {
    # Split the README content and insert the new section
    before <- readme_content[1:(insertion_point - 1)]
    after <- readme_content[insertion_point:length(readme_content)]

    # Add a blank line for spacing
    new_readme_content <- c(before, "", summary_section, after)

    # 5. Write the updated content back to README.md
    writeLines(new_readme_content, readme_path)
    print("README.md has been updated with the cohort summary.")
}
