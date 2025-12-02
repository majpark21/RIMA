#' Calculate P-Values for Neighbourhood-Neighbourhood Similarities
#'
#' Internal function that calculates the statistical significance of neighbourhood-neighbourhood
#' similarity values by comparing them to similarities between original and scrambled neighbourhoods.
#'
#' @param dt_sim_true,dt_sim_scrambled data.table with 3 columns as returned by \code{calculate_similarities()}.
#'   The first 2 columns indicate neighbourhood pairs, the third column (\code{col_sim}) indicates similarity values.
#'   \code{dt_sim_true} contains similarities between original neighbourhoods,
#'   \code{dt_sim_scrambled} contains similarities between original and scrambled neighbourhoods.
#' @param alpha Numeric. Significance level for calling edges significant. Float between 0 and 1.
#'   Default is 0.05.
#' @param adjust Character string or \code{NULL}. Method for p-value adjustment.
#'   Must be a valid method in \code{p.adjust()} or \code{NULL} for no adjustment.
#'   Default is "holm".
#' @param col_sim Character string. Name of the column holding similarity values in both data.tables.
#'   Default is "sim".
#' @param col_group Character string or \code{NULL}. Name of a column in \code{dt_sim_true}
#'   specifying groups for null population estimation. If set to one of the neighbourhood columns,
#'   p-values are estimated separately for each neighbourhood (comparing scrambled pairs
#'   involving that specific neighbourhood). If \code{NULL}, uses all scrambled similarities.
#'   Default is \code{NULL}.
#'
#' @details
#' The function:
#' 1. Estimates a null distribution using the empirical CDF of scrambled similarities.
#' 2. Assigns empirical p-values to true similarities by comparing to the null distribution.
#' 3. Applies multiple testing correction if specified.
#' P-values represent the proportion of scrambled similarities that are >= true similarity.
#'
#' @returns A data.table based on \code{dt_sim_true} with additional columns:
#'   \item{pval}{Empirical p-value for each neighbourhood pair.}
#'   \item{pval_adjusted}{Adjusted p-value after multiple testing correction.}
#'
#' @keywords internal
#' @noRd
.get_nhoodnhood_pval <- function(dt_sim_true,
                                 dt_sim_scrambled,
                                 alpha = 0.05,
                                 adjust = "holm",
                                 col_sim = "sim",
                                 col_group = NULL) {
  # Basic function applied to each group to compare true similarities to scrambled similarities. From this returns pvalues, or the significance cutoff of similarity.
  #
  # First estimates null population by calculating ECDF of scrambled similarities. True similarities are compared to this null distribution and assigned a pvalue. The pvalue represent the proportion of scrambled similarities that are greater than the true similarity.
  # @param v_sim_true Numerical vector holding the similarities between original neighbourhoods.
  # @param v_sim_scrambled Numerical vector the similarities between original and scrambled neighbourhoods.
  # @param alpha Level of significance to call significant edges. Float between 0 and 1. Only used if return_cutoff is TRUE.
  # @param adjust Method to correct pvalues for multiple testing. If NULL does not do pvalue correction. Must be a valid method in 'p.adjust'.
  # @param return_cutoff If TRUE, instead of calculating pvalues, return similarity cutoff of significance.
  #
  # @returns A list of 2 vectors of same length as v_sim_true. The first represents the pvalues associated with the true similarities, the second represents the pvalues after adjustment.
  .compare_true_scrambled <- function(v_sim_true,
                                      v_sim_scrambled,
                                      alpha = 0.05,
                                      adjust = "holm",
                                      return_cutoff = FALSE) {
    # Special case speed-up: no need to calculate costly ecdf
    if ((return_cutoff) & (is.null(adjust))) {
      cutoff <- quantile(v_sim_scrambled, probs = 1 - alpha)
      return(cutoff)
    }
    # All other cases need to calculate ecdf
    fn_ecdf <- ecdf(v_sim_scrambled)  # Estimate p-values function using scrambled distribution
    pvals <- 1 - fn_ecdf(v_sim_true)  # Assign p-values to authentic distribution
    if (!is.null(adjust)) {
      pvals_adj <- p.adjust(pvals, method = adjust)
    } else {
      pvals_adj <- pvals
    }
    if (return_cutoff) {
      signif_sim <- which(pvals <= alpha)
      cutoff <- min(v_sim_true[signif_sim])
      return(cutoff)
    } else {
      return(list(pvals = pvals, pvals_adj = pvals_adj))
    }
  }

  cols_original <- colnames(dt_sim_true)
  dt_sim_true <- copy(dt_sim_true)  # Avoid modifying original object
  dt_sim_scrambled <- copy(dt_sim_scrambled)

  col_pval <- "pval"
  col_pval_adj <- "pval_adjusted"
  col_nhoods <- setdiff(colnames(dt_sim_true), col_sim)
  # Make sure nhoods are not integer because confuses grouping
  for (col in col_nhoods) {
    set(dt_sim_true, j = col, value = as.character(dt_sim_true[[col]]))
    set(dt_sim_scrambled, j = col, value = as.character(dt_sim_scrambled[[col]]))
  }

  if (is.null(col_group)) {
    tmp <- .compare_true_scrambled(
      v_sim_true = dt_sim_true[[col_sim]],
      v_sim_scrambled = dt_sim_scrambled[[col_sim]],
      alpha = alpha,
      adjust = adjust,
      return_cutoff = FALSE
    )
    dt_sim_true[[col_pval]] <- tmp$pvals
    dt_sim_true[[col_pval_adj]] <- tmp$pvals_adj
  } else {
    setkeyv(dt_sim_true, col_group)
    setkeyv(dt_sim_scrambled, col_group)

    groups <- unique(dt_sim_true[[col_group]])
    pb <- txtProgressBar(min = 0,
                         max = length(groups),
                         style = 3)
    for (idx_groups in seq(length(groups))) {
      curr_group <- groups[idx_groups]
      tmp <- .compare_true_scrambled(
        v_sim_true = dt_sim_true[curr_group, get(col_sim)],
        v_sim_scrambled = dt_sim_scrambled[curr_group, get(col_sim)],
        alpha = alpha,
        adjust = adjust,
        return_cutoff = FALSE
      )
      dt_sim_true[curr_group, (col_pval) := tmp$pvals]
      dt_sim_true[curr_group, (col_pval_adj) := tmp$pvals_adj]
      setTxtProgressBar(pb, idx_groups)
    }
    close(pb)
  }

  return(dt_sim_true)
}
