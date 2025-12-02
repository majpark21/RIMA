#' Calculate Neighbourhood-Neighbourhood Significance (One Direction)
#'
#' Internal function that performs significance testing in one direction.
#' This function does the scrambling of neighbourhoods and calculates p-values against
#' scrambled similarities. It runs in only one direction (one Milo intact, the other scrambled).
#'
#' @param milos A list of 2 Milo objects. The 1st is kept intact, the 2nd is scrambled.
#' @param name_intact,name_scramble Character strings specifying the column names for the
#'   neighbourhoods from the 1st (intact) and 2nd (scrambled) Milo objects.
#' @param col_scramble_label Character string specifying a column in the 2nd Milo's \\code{colData}.
#'   Label frequencies are estimated from this column to weight cell sampling during scrambling.
#'   Cells with common labels are less likely to be sampled. Use \"false\" for unweighted sampling.
#' @param n_scrambles Integer. Number of rounds of cell resampling.
#' @param assay Character string specifying the assay to use for calculating neighbourhood expression.
#' @param sim_method Character string specifying the similarity metric.
#'   Must be one of \"pearson\", \"kendall\", or \"spearman\".
#' @param dt_sims_true A data.table with 3 columns as returned by \\code{calculate_similarities()}.
#'   The first 2 columns indicate neighbourhood pairs, the third column (named \"sim\")
#'   contains similarity values between intact neighbourhoods.
#' @param adjust Character string specifying the p-value adjustment method.
#'   Must be a valid method in \\code{p.adjust()}. If \\code{NULL}, no adjustment is applied.
#' @param alpha_adjust Numeric. Significance level for calling edges significant.
#'
#' @details
#' This function:
#' 1. Scrambles neighbourhoods in the 2nd Milo, optionally using label weights.
#' 2. Calculates similarities between intact nhoods from 1st Milo and scrambled nhoods from 2nd.
#' 3. Calculates and adjusts p-values by comparing true to scrambled similarities.
#'
#' @returns A data.table with columns from \\code{dt_sims_true} plus:
#'   \\item{pval}{Empirical p-value for each neighbourhood pair.}
#'   \\item{pval_adjusted}{Adjusted p-value after multiple testing correction.}
#'   \\item{is_significant}{Logical indicating significance at \\code{alpha_adjust}.}
#'
#' @keywords internal
#' @noRd
.calculate_nhoodnhood_significance_oneway <- function(milos,
                                                      dt_sims_true,
                                                      n_scrambles,
                                                      name_intact = "milo_intact",
                                                      name_scramble = "milo_scramble",
                                                      col_scramble_label = "false",
                                                      assay = "logcounts",
                                                      sim_method = "spearman",
                                                      adjust="holm",
                                                      alpha_adjust=0.05) {
  # Generate null (scrambled) distribution of similarities, upsample rare labels through resampling weights
  l_sims_scrambled <- list()

  # Weight is inversely proportional to the number of cells of the corresponding celltype
  weight_cells <- .scramble_calculateWeights(milo = milos[[2]], col_scramble_label = col_scramble_label)
  # Generate scrambled nhood average expressions
  cat("Generate scrambled neighbourhoods...")
  l_scrambled_assay_nhoods <- .scramble_nhoodExpression(
    milos[[2]],
    # already subset here for speed up
    replace = TRUE,
    n = n_scrambles,
    weight_sample = weight_cells,
    assay = assay
  )

  # Calculate similarities against scrambled nhood average expression
  cat("Calculate similarities against scrambled nhood average expression...\n")
  l_matExpression <- list(nhoodExpression(milos[[1]]))
  for (ii in 1:n_scrambles) {
    l_matExpression[[2]] <- l_scrambled_assay_nhoods[[ii]]
    dt_sims_scrambled <- calculate_similarities(milos = l_matExpression, method = sim_method)

    setnames(dt_sims_scrambled, c(name_intact, name_scramble, "sim"))
    l_sims_scrambled[[ii]] <- dt_sims_scrambled
  }
  dt_sims_scrambled <- rbindlist(l_sims_scrambled, idcol = "iteration")

  # Filter insignificant nhood-nhood edges
  cat("Get pvalues nhood-nhood edges...\n")
  dt_sims_withSignif <- .get_nhoodnhood_pval(
    dt_sim_true = dt_sims_true,
    dt_sim_scrambled = dt_sims_scrambled,
    alpha = alpha_adjust,
    adjust = adjust,
    col_sim = "sim",
    col_group = name_intact
  )

  return(dt_sims_withSignif)
}
