#' Heatmap: how many edges linking were trimmed by RIMA?
#'
#' This function shows how many nhood-nhood edges were trimmed by RIMA. The nhoods are first annotated to calculate the number of (un)significant edges linking 2 neighbourhoods' annotation.
#' @param milos A list of 2 Milo objects.
#' @param dt_sims_withSignif A data.table with at least 3 columns. Each row represents a pair of nhoods across the Milos (column 1 and 2). A logical column named 'is_significant'
#' @param cols_label Vector of length 2, names of columns in the Milo's colData. Neighbourhoods will be annotated accroding to these columns and the number of (trimmed) edges will be calculated between the annotations.
#' @param cluster_rows,cluster_cols Logical, whether to reorder the rows and columns with hierarchical clustering.
#'
#' @returns A ggplot object
#' @export
#'
#' @examples
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
