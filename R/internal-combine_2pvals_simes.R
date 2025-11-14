#' Simes pvalue combination special case with 2 values
#'
#' Use non-generalizable Simes for speed up and ease of use in data.table
#' @param pval1,pval2 The two pvalues to be combined.
#'
#' @returns A single, combined pvalue.
.combine_2pvals_simes <- function(pval1, pval2){
  pvals <- c(pval1, pval2)
  ordered_pvals <- sort(pvals)
  adj_pvals <- 2 *(ordered_pvals/(1:2))
  return(min(1, adj_pvals))
}
