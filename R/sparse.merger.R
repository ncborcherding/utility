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

sparse.expander <- function(sparse.matrix, 
                            sc) {
  sparse.matrix  <- sparse.matrix [,intersect(colnames(sparse.matrix ), colnames(sc))]
  SC.barcodes <- colnames(sc)[which(colnames(sc) %!in% colnames(sparse.matrix))]
  emp.matrix <- Matrix(data = 0, 
                       nrow = nrow(sparse.matrix), 
                       ncol = length(SC.barcodes),
                       dimnames = list(rownames(sparse.matrix), SC.barcodes),
                       sparse = TRUE)
  sparse.matrix <- merge.sparse(list(sparse.matrix, emp.matrix))
  return(sparse.matrix)
}