#' Heatmap: Visualize Edges Trimmed by Significance Testing
#'
#' Displays a heatmap showing the proportion of neighbourhood-neighbourhood edges retained
#' after significance filtering. Neighbourhoods are annotated and aggregated, showing how
#' many edges linking annotation groups passed significance thresholds.
#'
#' @param milos A list of 2 Milo objects.
#' @param dt_sims_withSignif A data.table with at least 3 columns. Each row represents a pair
#'   of neighbourhoods across the two Milos (first two columns contain neighbourhood identifiers).
#'   Must include a logical column named \code{"is_significant"}.
#' @param cols_label Character vector of length 2 specifying column names from the Milo's
#'   \code{colData}. Neighbourhoods will be annotated according to these columns and the
#'   ratio of retained edges will be calculated between annotated groups.
#' @param cluster_rows Logical. If \code{TRUE}, reorders rows using hierarchical clustering.
#'   Default is \code{TRUE}.
#' @param cluster_cols Logical. If \code{TRUE}, reorders columns using hierarchical clustering.
#'   Default is \code{TRUE}.
#'
#' @returns A \code{ggplot} object displaying a heatmap with the fraction of retained edges.
#'
#' @export
#'
#' @examples
#' # Not run: plot_trimming_map(milos, dt_sims_withSignif)
plot_trimming_map <- function(milos,
                              dt_sims_withSignif,
                              cols_label = c("celltype", "celltype"),
                              cluster_rows = TRUE,
                              cluster_cols = TRUE) {
  # Check arguments, set defaults
  if (length(cols_label) != 2) {
    stop("cols_label must either be of length 2.")
  }

  name1 <- colnames(dt_sims_withSignif)[1]
  name2 <- colnames(dt_sims_withSignif)[2]
  dt_n_edges <- .calculate_trimmed_ratios(milos, dt_sims_withSignif, cols_label)

  cols <- c(paste0(cols_label[1], "_", name1),
            paste0(cols_label[2], "_", name2),
            "ratio_retained")
  # Determine row order (hclust)
  if(cluster_rows){
    tmp <- copy(dt_n_edges[, ..cols])
    cast_formula <- as.formula(paste0(paste0(cols_label[1], "_", name1), " ~ ", paste0(cols_label[2], "_", name2)))
    tmp <- as.matrix(dcast(tmp, cast_formula, value.var = "ratio_retained"))
    rn <- tmp[, paste0(cols_label[1], "_", name1)]
    tmp <- tmp[, colnames(tmp) != paste0(cols_label[1], "_", name1)]
    tmp <- apply(tmp, 2, as.numeric)
    rownames(tmp) <- rn
    tmp[which(is.na(tmp))] <- 0
    ord_rows <- hclust(dist(tmp, method = "euclidean"), method = "ward.D")
    col <- paste0(cols_label[1], "_", name1)
    dt_n_edges[, (col) := factor(get(col), levels = rownames(tmp)[ord_rows$order])]
  }
  # Determine col order (hclust)
  if(cluster_cols){
    tmp <- copy(dt_n_edges[, ..cols])
    cast_formula <- as.formula(paste0(paste0(cols_label[2], "_", name2), " ~ ", paste0(cols_label[1], "_", name1)))
    tmp <- as.matrix(dcast(tmp, cast_formula, value.var = "ratio_retained"))
    rn <- tmp[, paste0(cols_label[2], "_", name2)]
    tmp <- tmp[, colnames(tmp) != paste0(cols_label[2], "_", name2)]
    tmp <- apply(tmp, 2, as.numeric)
    rownames(tmp) <- rn
    tmp[which(is.na(tmp))] <- 0
    ord_cols <- hclust(dist(tmp, method = "euclidean"), method = "ward.D")
    col <- paste0(cols_label[2], "_", name2)
    dt_n_edges[, (col) := factor(get(col), levels = rownames(tmp)[ord_cols$order])]
  }

  p <- ggplot(dt_n_edges, aes(x = .data[[paste0(cols_label[2], "_", name2)]], y = .data[[paste0(cols_label[1], "_", name1)]])) +
    geom_tile(aes(fill = ratio_retained), color="black") +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    viridis::scale_fill_viridis(option = "E") +
    labs(fill = "Fraction retained edges") +
    theme_minimal()

  return(p)
}
