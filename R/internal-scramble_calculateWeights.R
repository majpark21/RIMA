#' Get the weights to randomly permute cell identity and create scrambled neighbourhoods
#'
#' The weight of sampling each cell is based on the inverse frequency of its associated label.
#' @param milo A Milo object with filled nhoods slot.
#' @param col_scramble_label The name of a column in the Milo's colData. Label frequencies are estimated from this column. If "false", returns NULL weights.
#'
#' @returns A numerical vector of same length as the number of cells in milo. NULL if col_scramble_label=="false".
#'
#' @examples
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
