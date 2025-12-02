#' Plot Paired Expression Across Matched Neighbourhoods
#'
#' Creates scatter plots comparing gene expression values between pairs of matched neighbourhoods
#' from two datasets. Each gene is displayed in a separate facet.
#'
#' @param milos List of 2 Milo objects with a filled nhoodExpression slot.
#' @param dt_match A data.table with 2 columns containing the nhood-nhood matching.
#'   The first column contains neighbourhood names from the 1st Milo object,
#'   the second column contains neighbourhood names from the 2nd Milo object.
#' @param genes A character vector of gene names to plot. Genes must be present in both
#'   Milo objects' rownames.
#'
#' @returns A \code{ggplot} object with scatter plots of paired expression values.
#'   Each gene is displayed in a separate facet.
#'
#' @export
#'
#' @examples
#' # Not run: plot_paired_expression(milos, dt_match, genes = c("gene1", "gene2", "gene3"))
plot_paired_expression <- function(milos, dt_match, genes){
  dt_plot <- RIMA::get_paired_expression(milos, dt_match, genes)
  # Plot genes in the same order as they were passed
  dt_plot[, gene := factor(gene, levels=genes)]
  p <- ggplot(dt_plot, aes(x=expression_1, y=expression_2)) +
    geom_point() +
    facet_wrap("gene", scales="free") +
    theme_bw()

  return(p)
}
