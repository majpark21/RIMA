#' Correlation of Paired Expression
#'
#' Calculate the correlation of gene expression values between pairs of matched neighbourhoods
#' @param milos List of 2 Milo objects with a filled nhoodExpression slot.
#' @param dt_match A data.table with 2 columns containing the nhood-nhood matching. The first (resp. second) column contains the names of nhoods in the 1st (resp. 2nd) Milo object.
#' @param genes A vector of genes for which to calculate the CoPE score. The genes must be present in both Milo's rownames.
#' @param method Correlation method. Must be one of c('pearson', 'kendall', 'spearman').
#'
#' @returns A data table with 2 columns. Each row contains the gene name and its associated cope score.
#' @export
#'
#' @examples
calculate_cope <- function(milos, dt_match, genes, method="spearman"){
  if (!method %in% c("pearson", "kendall", "spearman")) {
    stop("'method' must be one of: 'pearson', 'kendall', 'spearman'.")
  }
  # Subset to genes of interest and reorder the nhoods according to matches, transpose to put genes in columns
  nhoods1 <- nhoodExpression(milos[[1]])
  nhoods1 <- t(nhoods1[genes, dt_match[[1]], drop=FALSE])

  nhoods2 <- nhoodExpression(milos[[2]])
  nhoods2 <- t(nhoods2[genes, dt_match[[2]], drop=FALSE])

  out <- cor(nhoods1, nhoods2, method = method)

  # Only keep same gene's correlation
  out <- diag(out)

  out <- data.table(gene=names(out), cope=out)

  return(out)
}
