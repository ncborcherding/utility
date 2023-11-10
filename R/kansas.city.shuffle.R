library(scran)
library(Matrix)
#Reorganize single-cell objects lists to generate a sparse matrix 
#list with the same genes organized by the top N most highly 
#variable genes. This is to minimize the memory foot print 
#before integration.

kansas.city.shuffle <- function(sce.list,
                                number.var.features = 2000,
                                assay = NULL, 
                                format = NULL) {
  
  if(is.null(format)) {
    if(inherits(x=sce.list[[1]], what ="Seurat")) {
      format <- "Seurat"
      if(is.null(assay)) {
        assay <- "RNA"
      }
    } else {
      format <- "SCE"
      if(is.null(assay)) {
        assay <- "logcounts"
      }
    }
  }
  #Removing other assay data 
  if(format == "Seurat") {
    list.tmp <- lapply(sce.list, function(x) {
      tmp <- as.SingleCellExperiment(x, assay = assay)
      tmp <- assay(x, "logcounts")
      tmp
    })
  } else {
    list.tmp <- lapply(sce.list, function(x) {
      tmp <- assay(x, assay)
      tmp
    })
  }
  for(i in seq_along(list.tmp)) {
    tmp.dim <- dim(list.tmp[[i]])[2]
    if(i == 1) {
      num.cells <- tmp.dim
    } else {
      num.cells <- c(num.cells, tmp.dim)
    }
  }
  #Combining the data
  combined.experiments <- merge.sparse(list.tmp)
  #Selecting the variable genes and subsetting the massive combined data
  var.features <- modelGeneVar(combined.experiments)
  genes <- getTopHVGs(var.features, n=number.var.features)
  combined.experiments <- combined.experiments[genes,]
  #remaking a list to be used with fastMNN
  new.sce.list <- lapply(num.cells, function(x) {
    matrix.subset <- combined.experiments[,1:x]
    matrix.subset
  })
  return(new.sce.list)
}

merge.sparse <- function(list) {
  RowMergeMatricesList <- getFromNamespace("RowMergeMatricesList", "SeuratObject")
  all.mat <- list
  all.colnames <- all.rownames <- vector(
    mode = 'list',
    length = length(x = all.mat)
  )
  for (i in seq_along(along.with = all.mat)) {
    if (is.data.frame(x = all.mat[[1]])) {
      all.mat[[i]] <- as.matrix(x = all.mat[[i]])
    }
    all.rownames[[i]] <- rownames(x = all.mat[[i]])
    all.colnames[[i]] <- colnames(x = all.mat[[i]])
  }
  use.cbind <- all(duplicated(x = all.rownames)[2:length(x = all.rownames)])
  if (isTRUE(x = use.cbind)) {
    new.mat <- do.call(what = cbind, args = all.mat)
  } else {
    all.mat <- lapply(X = all.mat, FUN = as, Class = "RsparseMatrix")
    all.names <- unique(x = unlist(x = all.rownames))
    new.mat <- RowMergeMatricesList(
      mat_list = all.mat,
      mat_rownames = all.rownames,
      all_rownames = all.names
    )
    rownames(x = new.mat) <- make.unique(names = all.names)
  }
  colnames(x = new.mat) <- make.unique(names = unlist(x = all.colnames))
  return(new.mat)
}

