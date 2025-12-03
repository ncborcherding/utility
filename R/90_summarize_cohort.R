# R/90_summarize_cohort.R

suppressPackageStartupMessages({
  library(Seurat)
  library(yaml)
  library(dplyr)
  library(data.table)
  library(knitr)
})

source("R/00_utils.R")

summarise_cohort <- function(config_path = "config.yaml") {
  log_message("--- Summarizing Cohort and Updating README ---")
  
  config <- yaml::read_yaml(config_path)
  
  # Use the Annotated object (has finalized labels)
  input_path <- file.path(config$paths$results_dir, "FINAL_integrated_object.rds")
  
  if (!file.exists(input_path)) {
    # Fallback to integrated if annotated doesn't exist yet
    input_path <- file.path(config$paths$results_dir, "FINAL_integrated_object.rds")
  }
  if (!file.exists(input_path)) stop("No final object found to summarize.")
  
  obj <- readRDS(input_path)
  
  # --- 1. Gather Statistics ---
  meta <- obj@meta.data
  
  total_cells <- nrow(meta)
  n_sequencing_runs <- length(unique(meta$orig.ident))
  
  # Parse Tissue and Patients if available
  # Check for 'tissue' column, default to NA
  tissues <- unique(sub("\\..*", "", meta$orig.ident))
  n_tissues <- length(tissues)
  
  # Infer patients from orig.ident (assuming format like PatientID_RunID or similar)
  # Adjust regex based on your specific naming convention
  # Assuming: ProjectID_PatientID_...
  patients <- unique(sub("^([^\\.]+\\.[^\\.]+\\.[^\\.]+)\\..*$", "\\1", meta$orig.ident)) 
  n_patients <- length(patients)
  
  # VDJ Statistics
  has_vdj <- "CTaa" %in% colnames(meta)
  n_with_tcr <- if(has_vdj) sum(!is.na(meta$CTaa)) else 0
  
  # --- 2. Generate Summary Table (Markdown) ---
  summary_df <- data.frame(
    Metric = c("Total Cells", "Sequencing Runs", "Unique Tissues", "Unique Patients", "Cells with TCR"),
    Count = c(
      format(total_cells, big.mark = ","),
      format(n_sequencing_runs, big.mark = ","),
      format(n_tissues, big.mark = ","),
      format(n_patients, big.mark = ","),
      format(n_with_tcr, big.mark = ",")
    )
  )
  
  table_md <- knitr::kable(summary_df, format = "markdown", align = "l")
  
  summary_section <- paste(
    "### Cohort Summary",
    "",
    paste0("Last Updated: ", Sys.Date()),
    "",
    paste(table_md, collapse = "\n"),
    "",
    sep = "\n"
  )
  
  # --- 3. Update README.md ---
  readme_path <- "README.md"
  if (file.exists(readme_path)) {
    content <- readLines(readme_path)
    
    # Define delimiters for the summary section to make future updates scalable
    start_mark <- "<!-- COHORT_SUMMARY_START -->"
    end_mark <- "<!-- COHORT_SUMMARY_END -->"
    
    if (any(grepl(start_mark, content))) {
      # Replace existing block
      start_idx <- grep(start_mark, content)[1]
      end_idx <- grep(end_mark, content)[1]
      
      new_content <- c(
        content[1:(start_idx)],
        summary_section,
        content[end_idx:length(content)]
      )
    } else {
      # Append after "Cohort Information"
      insert_idx <- grep("#### Cohort Information", content)[1]
      if (is.na(insert_idx)) insert_idx <- length(content)
      
      new_content <- c(
        content[1:insert_idx],
        "",
        start_mark,
        summary_section,
        end_mark,
        "",
        content[(insert_idx+1):length(content)]
      )
    }
    
    writeLines(new_content, readme_path)
    log_message("README.md updated.")
  } else {
    log_message("README.md not found, skipping update.")
  }
  
  # --- 4. Export Detailed Tables ---
  summary_dir <- file.path(config$paths$results_dir, "10_summary_tables")
  dir.create(summary_dir, showWarnings = FALSE)
  
  # Cells per tissue
  if("tissue" %in% colnames(meta)) {
    tissue_stats <- meta %>% group_by(tissue) %>% tally()
    fwrite(tissue_stats, file.path(summary_dir, "cells_per_tissue.csv"))
  }
  
  # Cells per cluster
  cluster_col <- "consensus_celltype"
  if(cluster_col %in% colnames(meta)) {
    cluster_stats <- meta %>% group_by(!!sym(cluster_col)) %>% tally()
    fwrite(cluster_stats, file.path(summary_dir, "cells_per_cluster.csv"))
  }
  
  log_message("Summary tables exported to ", summary_dir)
}