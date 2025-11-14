#' Pipeline to give pvalue to nhood-nhood edges
#'
#' This function takes in both Milos objects and the similarities between  does the scrambling of nhoods, and calculates the pvalues against scrambled similarities.
#' @param milos A list of 2 Milo objects, the 1st one is kept intact, the 2nd one is scrambled.
#' @param dt_sims A data.table with 3 columns. The first 2 columns indicate the pair of neighbourhoods, the third indicates the similarity value. This contains the similarity values between intact neighbourhoods.
#' @param col_scramble_label The name of a column in the 2nd Milo's colData. Label frequencies are estimated from this column and used in the scrambling process to weight the sampling probability of each cell.. Cells with most common samples are less likely to be sampled.
#' @param n_scrambles Number of rounds of cell resampling.
#' @param assay Assay from milos from which to calculate the neighbourhood expression.
#' @param sim_method Similarity metric. Must be one of c('pearson', 'kendall', 'spearman').
#' @param adjust Method to correct pvalues for multiple testing. If NULL does not do pvalue correction. Must be a valid method in 'p.adjust'.
#' @param alpha_adjust Level of significance to call significant edges. Float between 0 and 1.
#' @param direction "Direction of the matching. Must be one of c('lr', 'rl', 'b')
#' If 'lr' ('lr'; short for left-right), nhoods in 1st Milo of milos (i.e. the left object) are kept intact and nhoods in 2nd Milo of milos (i.e. the right object) are scrambled to derive the significance of nhood-nhood similarities.
#' If 'rl' (right-left), this role of the 1st and 2nd Milo are inverted.
#' If 'b' (bidirectional), the nhood-nhood significance is calculated in both directions and combined with Simes' method.
#' The 'b' option is the most conservative and slower (because it runs both lr and rl) but accounts for the difference of heterogeneity between the 2 milos and is most reliable.",
#'
#' @returns A data.table based on 'dt_sim_true', with 2 extra columns representing the pvalue and the adjusted pvalue associated with each similarity value.
#' @export
#'
#' @examples
calculate_nhoodnhood_significance <- function(milos,
                                              dt_sims,
                                              n_scrambles,
                                              col_scramble_label = "false",
                                              assay = "logcounts",
                                              sim_method = "spearman",
                                              adjust = "holm",
                                              alpha_adjust = 0.05,
                                              direction = c('b', 'lr', 'rl')) {
  direction <- match.arg(direction)
  # Arguments for .calculate_nhoodnhood_significance_oneway() common to all calls, no matter the direction
  ls_args_common <- list(
    col_scramble_label = col_scramble_label,
    n_scrambles = n_scrambles,
    assay = assay,
    sim_method =  sim_method,
    adjust = adjust,
    alpha_adjust = alpha_adjust
  )
  name_milo1 <- colnames(dt_sims)[1]
  name_milo2 <- colnames(dt_sims)[2]

  # Left -> Right: 1st Milo intact, 2nd Milo scrambled ============================================
  if (direction == 'lr') {
    cols_sim <- c(name_milo1, name_milo2, "sim")
    ls_args <- list(
      milos = list(milos[[1]], milos[[2]]),
      name_intact = name_milo1,
      name_scramble = name_milo2,
      dt_sims_true = dt_sims[, ..cols_sim]
    )
    ls_args <- append(ls_args, ls_args_common)
    dt_sims_withSignif <- do.call(.calculate_nhoodnhood_significance_oneway, ls_args)
  }
  # Right -> Left: 2nd Milo intact, 1st Milo scrambled ============================================
  else if (direction == 'rl') {
    cols_sim <- c(name_milo2, name_milo1, "sim")
    ls_args <- list(
      milos = list(milos[[2]], milos[[1]]),
      name_intact = name_milo2,
      name_scramble = name_milo1,
      dt_sims_true = dt_sims[, ..cols_sim]
    )
    ls_args <- append(ls_args, ls_args_common)
    dt_sims_withSignif <- do.call(.calculate_nhoodnhood_significance_oneway, ls_args)
  }
  # Both directions: does lr, then rl, and combine the pvalues from both ==========================
  else if (direction == 'b') {
    ls_pvals <- list()

    # left -> right
    cols_sim <- c(name_milo1, name_milo2, "sim")
    ls_args <- list(
      milos = list(milos[[1]], milos[[2]]),
      name_intact = name_milo1,
      name_scramble = name_milo2,
      dt_sims_true = dt_sims[, ..cols_sim]
    )
    ls_args <- append(ls_args, ls_args_common)
    ls_pvals[[1]] <- do.call(.calculate_nhoodnhood_significance_oneway, ls_args)

    # right -> left
    cols_sim <- c(name_milo2, name_milo1, "sim")
    ls_args <- list(
      milos = list(milos[[2]], milos[[1]]),
      name_intact = name_milo2,
      name_scramble = name_milo1,
      dt_sims_true = dt_sims[, ..cols_sim]
    )
    ls_args <- append(ls_args, ls_args_common)
    ls_pvals[[2]] <- do.call(.calculate_nhoodnhood_significance_oneway, ls_args)

    # Return p-values from both directions; combine with Simes method and adjust the combined Pvalues
    cat("Combining left-right and right-left pvalues...")
    dt_sims_withSignif <- merge(
      ls_pvals[[1]],
      ls_pvals[[2]],
      by = c(name_milo1, name_milo2),
      suffixes = c("_lr", "_rl")
    )
    dt_sims_withSignif[, pval_combined := .combine_2pvals_simes(pval_lr, pval_rl), by =
                         1:nrow(dt_sims_withSignif)]
    dt_sims_withSignif[, pval_adjusted_combined := p.adjust(pval_combined, method =
                                                              adjust)]
    dt_sims_withSignif[, is_significant := pval_adjusted_combined <= alpha_adjust]

    # Consistent column names and order with unidirectional
    dt_sims_withSignif[, sim := sim_lr]  # duplicated sim_lr and sim_rl if sim_method is symmetrical
    setcolorder(dt_sims_withSignif, c(name_milo1, name_milo2, "sim"))
  }

  return(dt_sims_withSignif)
}
