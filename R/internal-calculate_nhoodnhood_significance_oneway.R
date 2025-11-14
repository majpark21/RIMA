#' Pipeline to give pvalue to nhood-nhood edges, oneway, not exposed to user
#'
#' This function does the scrambling of nhoods, and calculates the pvalues against scrambled similarities. Runs only in one direction, i.e. one of the Milos is kept intact, the other is scrambled.
#' @param milos A list of 2 Milo objects, the 1st one is kept intact, the 2nd one is scrambled.
#' @param name_intact,name_scramble Name given to the column with neighbourhoods from 1st (resp. 2nd) Milo object.
#' @param col_scramble_label The name of a column in the 2nd Milo's colData. Label frequencies are estimated from this column and used in the scrambling process to weight the sampling probability of each cell.. Cells with most common samples are less likely to be sampled.
#' @param n_scrambles Number of rounds of cell resampling.
#' @param assay Assay from milos from which to calculate the neighbourhood expression.
#' @param sim_method Similarity metric. Must be one of c('pearson', 'kendall', 'spearman').
#' @param dt_sims_true A data.table with 3 columns. The first 2 columns indicate the pair of neighbourhoods, the third indicates the similarity value. This contains the similarity values between intact neighbourhoods.
#' @param adjust Method to correct pvalues for multiple testing. If NULL does not do pvalue correction. Must be a valid method in 'p.adjust'.
#' @param alpha_adjust Level of significance to call significant edges. Float between 0 and 1.
#'
#' @details
#' This function does:
#' 1. Scramble neighbourhoods in the 2nd Milo, optionally taking labels to determine the resampling weights.
#' 2. Calculate the similarities between true nhoods from the 1st Milo and the scrambled nhoods from the 2nd Milo
#' 3. Calculate and adjust pvalues for each nhood-nhood edges by comparing dt_sims_true to the scrambled ones.
#'
#' @returns A data.table with 3 columns
#'
#' @examples
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
  weight_cells <- .scramble_calculateWeights(milo = milos[[2]], col_scramble = col_scramble_label)
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
