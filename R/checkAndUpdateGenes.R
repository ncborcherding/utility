library(HGNChelper)
"%!in%" <- Negate("%in%")

gene.symbols <- HGNChelper::hgnc.table

checkAndUpdateGenes <- function(genes, gene.symbols) {
  # Initialize a vector to hold the updated gene symbols
  updated.genes <- vector("character", length = length(genes))
  
  # Initialize a vector to hold genes not found
  not.found.genes <- character(0)
  
  for (i in seq_along(genes)) {
    gene <- genes[i]
    
    # Check if the gene is in the Symbol column
    if (gene %in% gene.symbols$Symbol & gene %!in% gene.symbols$Approved.Symbol) {
      # Update to Approved.Symbol
      updated.genes[i] <- gene.symbols$Approved.Symbol[gene.symbols$Symbol == gene][1]
    } else if (gene %in% gene.symbols$Approved.Symbol) {
      # Gene is already an approved symbol
      updated.genes[i] <- gene
    } else {
      # Gene not found in either column, add to not.found.genes
      not.found.genes <- c(not.found.genes, gene)
      updated.genes[i] <- NA  # Mark as NA or keep original, depending on your preference
    }
  }
  
  # Return a list containing the updated genes and the not found genes
  list(updated_genes = updated.genes, not_found_genes = not.found.genes)
}
