#' Match Neighbourhoods Using Weighted Bipartite Graph Matching
#'
#' Identifies one-to-one matching between neighbourhood pairs from two datasets
#' using a weighted bipartite graph matching algorithm that maximizes the sum of similarities.
#'
#' @param dt_significant_nhoodnhood A data.table with at least 3 columns as returned by
#'   \code{calculate_nhoodnhood_significance()}. Must contain columns named "sim" and two neighbourhood
#'   identifier columns. The first two columns contain the neighbourhood pairs to be matched.
#' @param filter_insignificant Logical. If \code{TRUE}, filters rows in \code{dt_significant_nhoodnhood}
#'   where the column specified by \code{col_significance} is not \code{TRUE}. Default is \code{TRUE}.
#' @param col_significance Character string specifying the name of a logical column in
#'   \code{dt_significant_nhoodnhood} used for filtering rows (when \code{filter_insignificant = TRUE}).
#'   Default is \code{"is_significant"}.
#' @param return_sim Logical. If \code{TRUE}, adds a 3rd column to the output containing
#'   the similarity value of each matched nhood-nhood pair. Default is \code{TRUE}.
#'
#' @returns A data.table with neighbourhood pairs. The first two columns contain the matched
#'   neighbourhood names. If \code{return_sim = TRUE}, a third column \code{sim} contains
#'   the similarity scores.
#'
#' @export
#'
#' @examples
#' # Not run: dt_match <- match_nhoods(dt_sims_significant)
match_nhoods <- function(dt_significant_nhoodnhood,
                         filter_insignificant = TRUE,
                         col_significance = "is_significant",
                         return_sim = TRUE) {
  if (filter_insignificant) {
    dt_significant_nhoodnhood <- dt_significant_nhoodnhood[get(col_significance) == TRUE]
  }
  wbgm <- .weightedBipartiteGraphMatching(dt_significant_nhoodnhood,
                                          col_sim = "sim",
                                          return_sim = return_sim)
  return(wbgm)
}
