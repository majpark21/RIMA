#' Plot paired expression
#'
#' @param milos List of 2 Milo objects with a filled nhoodExpression slot.
#' @param dt_match A data.table with 2 columns containing the nhood-nhood matching. The first (resp. second) column contains the names of nhoods in the 1st (resp. 2nd) Milo object.
#' @param genes A vector of genes names. The genes must be present in both Milo's rownames.
#'
#' @returns A ggplot object, where each gene is plotted as a scatterplot in its own facet.
#' @export
#'
#' @examples
plot_paired_expression <- function(milos, dt_match, genes){
  dt_plot <- RIMA::get_paired_expression(milos, dt_match, genes)
  p <- ggplot(dt_plot, aes(x=expression_1, y=expression_2)) +
    geom_point() +
    facet_wrap("gene") +
    theme_bw()

  return(p)
}
