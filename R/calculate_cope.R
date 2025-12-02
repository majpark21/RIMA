#' Calculate Correlation of Paired Expression (CoPE)
#'
#' Calculates the correlation of gene expression values between pairs of matched neighbourhoods
#' to identify conserved genes across datasets.
#'
#' @param milos List of 2 Milo objects with a filled nhoodExpression slot.
#' @param dt_match A data.table with 2 columns containing the nhood-nhood matching.
#'   The first column contains neighbourhood names from the 1st Milo object,
#'   the second column contains neighbourhood names from the 2nd Milo object.
#' @param genes A character vector of gene names for which to calculate the CoPE score.
#'   Genes must be present in both Milo objects' rownames.
#'   If \code{NULL}, uses all available genes.
#' @param method Character string specifying the correlation method.
#'   Must be one of \code{"pearson"}, \code{"kendall"}, or \code{"spearman"}.
#'   Default is \code{"spearman"}.
#'
#' @returns A data.table with 2 columns:
#'   \item{gene}{Gene name.}
#'   \item{cope}{Correlation of paired expression score.}
#'
#' @export
#'
#' @examples
#' # Not run: calculate_cope(milos, dt_match, genes = c("gene1", "gene2"))
calculate_cope <- function(milos, dt_match, genes=NULL, method="spearman"){
  if (!method %in% c("pearson", "kendall", "spearman")) {
    stop("'method' must be one of: 'pearson', 'kendall', 'spearman'.")
  }
  if(is.null(genes)){
    genes <- rownames(milos[[1]])
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
