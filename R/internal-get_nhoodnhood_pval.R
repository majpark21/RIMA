#' Calculate the significance of nhood-nhood similarity values, by comparing them to similarities between original and scrambled nhoods
#'
#' First estimates null population by calculating ECDF of scrambled similarities. True similarities are compared to this null distribution and assigned a pvalue. The pvalue represent the proportion of scrambled similarities that are greater than the true similarity.
#' @param dt_sim_true,dt_sim_scrambled data.table with 3 columns, as returned by calculate_similarities(). The first 2 columns indicate the pair of neighbourhoods, the third indicates the similarity value. dt_sim_true holds the similarities between original neighbourhoods, dt_sim_scrambled holds the similarities between original and scrambled neighbourhoods.
#' @param alpha  Level of significance to call significant edges. Float between 0 and 1.
#' @param adjust Method to correct pvalues for multiple testing. If NULL does not do pvalue correction. Must be a valid method in 'p.adjust'.
#' @param col_sim Name of the column holding the similarity values in dt_sim_true and dt_sim_scrambled.
#' @param col_group Name of a column in dt_sim_true. Specifies the groups that are used to define the null population in pvalue calculation. If named after one of the two first columns containing neighbourhood names, the pvalues will be estimated for each neighbourhood individually. Specifically, this means that all true nhood1-nhoodX similarity values will be compared to the similarity values of all scrambled pairs involving nhood1 (i.e. nhood1-scrambled1, nhood1-scrambled2...). If NULL, the scrambled similarities across every pair (not just those containing nhood1) will be used to estimate the NULL population (not recommended).
#'
#' @returns A data.table based on 'dt_sim_true', with 2 extra columns representing the pvalue and the adjusted pvalue associated with each similarity value.
#'
#' @examples
.get_nhoodnhood_pval <- function(dt_sim_true,
                                 dt_sim_scrambled,
                                 alpha = 0.05,
                                 adjust = "fdr",
                                 col_sim = "sim",
                                 col_group = NULL) {
  #' Basic function applied to each group to compare true similarities to scrambled similarities. From this returns pvalues, or the significance cutoff of similarity.
  #'
  #' First estimates null population by calculating ECDF of scrambled similarities. True similarities are compared to this null distribution and assigned a pvalue. The pvalue represent the proportion of scrambled similarities that are greater than the true similarity.
  #' @param v_sim_true Numerical vector holding the similarities between original neighbourhoods.
  #' @param v_sim_scrambled Numerical vector the similarities between original and scrambled neighbourhoods.
  #' @param alpha Level of significance to call significant edges. Float between 0 and 1. Only used if return_cutoff is TRUE.
  #' @param adjust Method to correct pvalues for multiple testing. If NULL does not do pvalue correction. Must be a valid method in 'p.adjust'.
  #' @param return_cutoff If TRUE, instead of calculating pvalues, return similarity cutoff of significance.
  #'
  #' @returns A list of 2 vectors of same length as v_sim_true. The first represents the pvalues associated with the true similarities, the second represents the pvalues after adjustment.
  .compare_true_scrambled <- function(v_sim_true,
                                      v_sim_scrambled,
                                      alpha = 0.05,
                                      adjust = "fdr",
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
    dt_sim_true[[col]] <- as.character(dt_sim_true[[col]])
    dt_sim_scrambled[[col]] <- as.character(dt_sim_scrambled[[col]])
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
