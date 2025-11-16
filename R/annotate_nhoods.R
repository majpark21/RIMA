#' Annotate neighbourhoods
#'
#' Do a majority voting annotations of neighbourhoods, i.e. the most label among
#' the neighbourhood's cells wins.
#' @param milo A Milo object with a filled nhoods slot.
#' @param cols_annot The columns from colData to use for annotation.
#'
#' @details This is a wrapper around miloR::annotateNhoods
#'
#' @returns A data.table with the neighbourhoods' annotations. The two first
#'   columns contain a nhood label and the name of the nhood. Then a pair of
#'   column for each annotation indicating the winning label, and the fraction
#'   of cells in the neighbourhood with the winning label.
#' @export
#'
#' @examples
annotate_nhoods <- function(milo, cols_annot = c("celltype")) {
  nhoods_mat <- nhoods(milo)
  nhood_annot <- data.table(Nhood = 1:ncol(nhoods_mat),
                            Nhood_center = colnames(nhoods_mat))
  l_annot <- lapply(cols_annot, function(x)
    miloR::annotateNhoods(milo, nhood_annot, coldata_col = x))
  merge_by <- function(x, y) {
    merge(x, y, by = c("Nhood", "Nhood_center"))
  }
  dt_annot <- Reduce(merge_by, l_annot)
  return(dt_annot)
}
