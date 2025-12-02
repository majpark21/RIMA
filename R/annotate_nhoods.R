#' Annotate Neighbourhoods with Majority Voting
#'
#' Performs majority voting annotation of neighbourhoods, where the most common label
#' among a neighbourhood's cells wins. This is a wrapper around \code{miloR::annotateNhoods}.
#'
#' @param milo A Milo object with a filled nhoods slot.
#' @param cols_annot Character vector of column names from \code{colData} to use for annotation.
#'
#' @details
#' For each neighbourhood and each annotation column, the function identifies the most
#' frequent label among cells in that neighbourhood and records the fraction of cells
#' with the winning label.
#'
#' @returns A data.table with the following columns:
#'   \item{Nhood}{Neighbourhood identifier (numeric).}
#'   \item{Nhood_center}{Centre cell name of the neighbourhood.}
#'   For each annotation column, two columns are added:
#'   \item{<col>}{The winning annotation label for that neighbourhood.}
#'   \item{<col>.fraction}{The fraction of cells with the winning label.}
#'
#' @export
#'
#' @examples
#' # Not run: annotate_nhoods(milo, cols_annot = c("celltype", "stage"))
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
