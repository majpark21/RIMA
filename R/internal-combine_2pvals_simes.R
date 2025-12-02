#' Combine Two P-Values Using Simes Method
#'
#' Internal function that combines two p-values using a special case of Simes' method
#' optimized for speed and ease of use within data.table operations.
#'
#' @param pval1,pval2 Numeric variables containing the two p-values to be combined.
#'
#' @returns A numeric value containing the combined p-value using Simes' method.
#'  Values are capped at 1.0.
#'
#' @details
#' Simes' method for combining two p-values is applied as follows:\n
#' 1. Sort the two p-values.\n
#' 2. Calculate: combined = min(1, 2 * min(p1/1, p2/2))\n
#'
#' @keywords internal
#' @noRd
.combine_2pvals_simes <- function(pval1, pval2){
  pvals <- c(pval1, pval2)
  ordered_pvals <- sort(pvals)
  adj_pvals <- 2 *(ordered_pvals/(1:2))
  return(min(1, adj_pvals))
}
