#' Scramble cell identity in a Milo object, and calculate scrambled nhoodExpression
#'
#' @param milo A Milo object with filled nhoods slot.
#' @param n Number of permutations.
#' @param replace Is the resampling to be done with replacement? If FALSE, this is a simple permutation.
#' @param weight_sample A vector of same length as number of cells that defines how likely a cell is to be picked. If specified, sampling is always done with replacement.
#' @param assay Assay to calculate neighbourhoods' expression.
#'
#' @returns A list of length n, where each entry is a matrix of dimensions [n_cell, n_nhood], containing the nhoodExpression of scrambled neighbourhoods.
#' @export
#'
#' @examples
#' @noRd
.scramble_nhoodExpression <- function(milo, n=5, replace=FALSE, weight_sample=NULL, assay="logcounts"){
    # Returns number of cells in each nhood
    .nhoodSize <- function(milo){
      out <- Matrix::colSums(nhoods(milo))
      return(out)
    }

  # Permute (or resample) cell gene counts, but keep original nhood structure (number nhoods, size nhoods...)
  mat <- SummarizedExperiment::assay(milo, assay)
  mat_nhoods <- nhoods(milo)
  all_cells <- colnames(milo)
  if(!is.null(weight_sample)){replace <- TRUE}
  nhood_sizes <- .nhoodSize(milo)

  l_scrambled_nhoodExpression <- list()
  for(ii in 1:n){
    permut_cells <- sample(
      all_cells,
      size = length(all_cells),
      replace = replace,
      prob = weight_sample
    )
    permut_mat <- mat[, permut_cells]
    colnames(permut_mat) <- all_cells

    # Average gene expression in permuted nhoods, use consistent matrix format with NhoodExpression
    l_scrambled_nhoodExpression[[ii]] <- permut_mat %*% mat_nhoods
    l_scrambled_nhoodExpression[[ii]] <- l_scrambled_nhoodExpression[[ii]] / rep(nhood_sizes, each=nrow(l_scrambled_nhoodExpression[[ii]]))
    l_scrambled_nhoodExpression[[ii]] <- as.matrix(l_scrambled_nhoodExpression[[ii]])
  }
  return(l_scrambled_nhoodExpression)
}
