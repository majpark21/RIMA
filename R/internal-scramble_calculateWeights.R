#' Calculate Cell Sampling Weights for Scrambling
#'
#' Internal function that calculates weights for cell resampling during scrambling,
#' based on the inverse frequency of a specified label to balance sampling across groups.
#'
#' @param milo A Milo object with filled nhoods slot.
#' @param col_scramble_label Character string specifying a column in the Milo's \code{colData}.
#'   Label frequencies are estimated from this column. If \"false\", returns \code{NULL} weights
#'   for unweighted uniform sampling. Default is \"false\".
#'
#' @returns A numeric vector of the same length as the number of cells in \code{milo}.
#'   Each value represents the sampling weight for that cell based on the inverse frequency
#'   of its label. Weights are normalized to sum to 1.
#'   Returns \code{NULL} if \code{col_scramble_label == \"false\"}.
#'
#' @details
#' Weights are calculated as the inverse of the label frequency, divided by the sum of all
#' inverse frequencies to normalize. This upsamples cells from rare labels and downsamples
#' cells from common labels, creating a balanced representation during scrambling.
#'
#' @keywords internal
#' @noRd
.scramble_calculateWeights <- function(milo, col_scramble_label="false"){
  if(col_scramble_label=="false"){
    weight_cells <- NULL
  } else {
    weight_cells <- 1/table(SummarizedExperiment::colData(milo)[[col_scramble_label]])
    weight_cells <- weight_cells / sum(weight_cells)
    weight_cells <- weight_cells[milo[[col_scramble_label]]]
  }
  return(weight_cells)
}
