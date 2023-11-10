PrepSCTFindMarkers <- function(object, assay = "SCT", verbose = TRUE) {
  med.tmp <- SCTResults(object = object[[assay]], slot = "cell.attributes")
  observed_median_umis <-  median(med.tmp[,"umi"])
  model.list <- slot(object = object[[assay]], name = "SCTModel.list")
  median_umi.status <- lapply(X = model.list,
                              FUN = function(x) { return(tryCatch(
                                expr = slot(object = x, name = 'median_umi'),
                                error = function(...) {return(NULL)})
                              )})
  if (any(is.null(x = unlist(x = median_umi.status)))){
    # For old SCT objects  median_umi is set to median umi as calculated from obserbed UMIs
    slot(object = object[[assay]], name = "SCTModel.list") <- lapply(X = model.list,
                                                                     FUN = UpdateSlots)
    SCTResults(object = object[[assay]], slot = "median_umi") <- observed_median_umis
    
  }
  model_median_umis <- SCTResults(object = object[[assay]], slot = "median_umi")
  min_median_umi <- min(unlist(x = observed_median_umis))
  if (verbose) {
    message(paste0("Found ",
                   length(x = levels(x = object[[assay]])),
                   " SCT models.",
                   " Recorrecting SCT counts using minimum median counts: ",
                   min_median_umi))
  }
  umi.assay <- unique(
    x = unlist(
      x = SCTResults(object = object[[assay]], slot = "umi.assay")
    )
  )
  if (length(x = umi.assay) > 1) {
    stop("Multiple UMI assays are used for SCTransform: ",
         paste(umi.assay, collapse = ", ")
    )
  }
  raw_umi <- GetAssayData(object = object, assay = umi.assay, slot = "counts")
  corrected_counts <- Matrix(
    nrow = nrow(x = raw_umi),
    ncol = ncol(x = raw_umi),
    data = 0,
    dimnames = dimnames(x = raw_umi),
    sparse = TRUE
  )
  cell_attr <- SCTResults(object = object[[assay]], slot = "cell.attributes")
  model_pars_fit <- SCTResults(object = object[[assay]], slot = "feature.attributes")[, c("theta", "(Intercept)", "log_umi")]
  arguments <- SCTResults(object = object[[assay]], slot = "arguments")
  model_str <- SCTResults(object = object[[assay]], slot = "model")
  set_median_umi <- rep(min_median_umi, length(levels(x = object[[assay]])))
  names(set_median_umi) <- levels(x = object[[assay]])
  set_median_umi <- as.list(set_median_umi)
  # correct counts
    model_genes <- rownames(x = model_pars_fit)
    x <- list(
      model_str = model_str,
      arguments = arguments,
      model_pars_fit = as.matrix(x = model_pars_fit),
      cell_attr = cell_attr
    )
    cells <- rownames(x = cell_attr)
    umi <- raw_umi[model_genes, cells]
    
    umi_corrected <- correct_counts(
      x = x,
      umi = umi,
      verbosity = 0,
      scale_factor = min_median_umi
    )
    corrected_counts[rownames(umi_corrected), colnames(umi_corrected)] <- umi_corrected
    
  corrected_data <- log1p(corrected_counts)
  suppressWarnings({object <- SetAssayData(object = object,
                                           assay = assay,
                                           slot = "counts",
                                           new.data = corrected_counts)})
  suppressWarnings({object <- SetAssayData(object = object,
                                           assay = assay,
                                           slot = "data",
                                           new.data = corrected_data)})
  SCTResults(object = object[[assay]], slot = "median_umi") <- set_median_umi
  
  return(object)
}

