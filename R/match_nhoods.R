#' Matches neighbourhoods with Weighted Bipartite Graph Matching
#'
#' @param dt_significant_nhoodnhood A data.table with at least 3 columns as returned by 'calculate_nhoodnhood_significance'. Must contain a column named "sim". The 2 first columns must be contain the nhoods to be matched.
#' @param filter_insignificant Logical, whether to filter the rows in 'dt_significant_nhoodnhood', where the column 'col_significance' is not TRUE.
#' @param col_significance If filter_insignificant is TRUE; name of a logical column in dt_significant_nhoodnhood used to filter the rows where its value is TRUE.
#' @param return_sim Logical. If TRUE, add a 3rd column to the output that contains the similarity value of the nhood-nhood pair.
#'
#' @returns A data.table with 2 columns containing the nhood-nhood matching. An additional column is present if return_sim is TRUE.
#' @export
#'
#' @examples
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
