#' Heatmap: Visualize Neighbourhood Matches by Annotation
#'
#' Displays a heatmap showing which neighbourhoods were matched together, aggregated by
#' their annotations. Neighbourhoods are annotated with \code{annotate_nhoods(), and
#' the number of matches between annotation labels is counted and displayed.
#'
#' @param milos A list of 2 Milo objects.
#' @param dt_match A data.table with at least 3 columns. Each row represents a matched pair
#'   of neighbourhoods across the two Milos (first two columns contain neighbourhood identifiers).
#' @param cols_label Character vector of length 2 specifying column names from the Milo's
#'   \code{colData}. Neighbourhoods will be annotated according to these columns and the
#'   number of matches will be calculated between the annotated groups.
#' @param cluster_rows Logical. If \code{TRUE}, reorders rows using hierarchical clustering.
#'   Default is \code{TRUE}.
#' @param cluster_cols Logical. If \code{TRUE}, reorders columns using hierarchical clustering.
#'   Default is \code{TRUE}.
#' @param transform Character string or \code{NULL}. Specifies a transformation to apply to
#'   match counts for visualization. Valid values include \code{"log10"}, \code{"log1p"}, \code{"sqrt"}.
#'   See \code{?ggplot2::scale_fill_continuous} for more options. Default is \code{NULL}.
#'
#' @returns A \code{ggplot} object displaying a heatmap of neighbourhood match counts.
#'
#' @export
#'
#' @examples
#' # Not run: plot_matches_map(milos, dt_match, cols_label = c("celltype", "celltype"))
plot_matches_map <- function(milos,
                              dt_match,
                              cols_label = c("celltype", "celltype"),
                              cluster_rows = TRUE,
                              cluster_cols = TRUE,
                         transform=NULL) {
  # Check arguments, set defaults
  if (length(cols_label) != 2) {
    stop("cols_label must either be of length 2.")
  }

  name1 <- colnames(dt_match)[1]
  name2 <- colnames(dt_match)[2]
  annot1 <- annotate_nhoods(milos[[1]], cols_label[1])
  annot2 <- annotate_nhoods(milos[[2]], cols_label[2])

  # Annotate the matches
  cols <- c("Nhood_center", cols_label[1])
  dt_match_annot <- merge(dt_match,
                    annot1[, ..cols],
                    by.x = name1,
                    by.y = "Nhood_center")
  cols <- c("Nhood_center", cols_label[2])
  dt_match_annot <- merge(
    dt_match_annot,
    annot2[, ..cols],
    by.x = name2,
    by.y = "Nhood_center",
    suffixes = c(paste0("_", name1), paste0("_", name2))
  )
  # Use consistent naming, mimic suffixes effect if cols are different
  if(cols_label[1] != cols_label[2]){
    setnames(
      dt_match_annot,
      cols_label,
      c(paste0(cols_label[1], "_", name1), paste0(cols_label[2], "_", name2))
    )
  }


  # Count the number of matches linking 2 labels
  dt_n_matches <- as.data.table(table(
    dt_match_annot[, get(paste0(cols_label[1], "_", name1))],
    dt_match_annot[, get(paste0(cols_label[2], "_", name2))]),
    keep.rownames = T)
  setnames(dt_n_matches, c(
    paste0(cols_label[1], "_", name1),
    paste0(cols_label[2], "_", name2),
    "n_matches"
  ))

  cols <- c(paste0(cols_label[1], "_", name1),
            paste0(cols_label[2], "_", name2),
            "n_matches")
  # Determine row order (hclust)
  if(cluster_rows){
    tmp <- copy(dt_n_matches[, ..cols])
    cast_formula <- as.formula(paste0(paste0(cols_label[1], "_", name1), " ~ ", paste0(cols_label[2], "_", name2)))
    tmp <- as.matrix(dcast(tmp, cast_formula, value.var = "n_matches"))
    rn <- tmp[, paste0(cols_label[1], "_", name1)]
    tmp <- tmp[, colnames(tmp) != paste0(cols_label[1], "_", name1)]
    tmp <- apply(tmp, 2, as.numeric)
    rownames(tmp) <- rn
    tmp[which(is.na(tmp))] <- 0
    ord_rows <- hclust(dist(tmp, method = "euclidean"), method = "ward.D")
    col <- paste0(cols_label[1], "_", name1)
    dt_n_matches[, (col) := factor(get(col), levels = rownames(tmp)[ord_rows$order])]
  }
  # Determine col order (hclust)
  if(cluster_cols){
    tmp <- copy(dt_n_matches[, ..cols])
    cast_formula <- as.formula(paste0(paste0(cols_label[2], "_", name2), " ~ ", paste0(cols_label[1], "_", name1)))
    tmp <- as.matrix(dcast(tmp, cast_formula, value.var = "n_matches"))
    rn <- tmp[, paste0(cols_label[2], "_", name2)]
    tmp <- tmp[, colnames(tmp) != paste0(cols_label[2], "_", name2)]
    tmp <- apply(tmp, 2, as.numeric)
    rownames(tmp) <- rn
    tmp[which(is.na(tmp))] <- 0
    ord_cols <- hclust(dist(tmp, method = "euclidean"), method = "ward.D")
    col <- paste0(cols_label[2], "_", name2)
    dt_n_matches[, (col) := factor(get(col), levels = rownames(tmp)[ord_cols$order])]
  }

  p <- ggplot(dt_n_matches, aes(x = .data[[paste0(cols_label[2], "_", name2)]], y = .data[[paste0(cols_label[1], "_", name1)]])) +
    geom_tile(aes(fill = n_matches), color="black") +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    labs(fill = "Number matches") +
    theme_minimal()
  if(!is.null(transform)){
    p <- p + viridis::scale_fill_viridis(option = "E", transform=transform)
  } else {
    p <- p + viridis::scale_fill_viridis(option = "E")
  }

  return(p)
}
