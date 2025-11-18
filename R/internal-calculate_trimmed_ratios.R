#' Calculate the number of nhood-nhood edges that were trimmed after comparing true similarities to scrambled ones.
#'
#' @param milos A list of 2 Milo objects.
#' @param dt_sims_withSignif A data.table with at least 3 columns. Each row represents a pair of nhoods across the Milos (column 1 and 2). A logical column named 'is_significant'
#' @param cols_label Vector of length 2, names of columns in the Milo's colData. Neighbourhoods will be annotated accroding to these columns and the number of (trimmed) edges will be calculated between the annotations.
#'
#' @returns A data.table with the number of (trimmed) edges.
#'
.calculate_trimmed_ratios <- function(milos,
                                      dt_sims_withSignif,
                                      cols_label = c("celltype", "celltype")) {
  name1 <- colnames(dt_sims_withSignif)[1]
  name2 <- colnames(dt_sims_withSignif)[2]
  annot1 <- annotate_nhoods(milos[[1]], cols_label[1])
  annot2 <- annotate_nhoods(milos[[2]], cols_label[2])

  # Annotate the matches
  cols <- c("Nhood_center", cols_label[1])
  dt_edges <- merge(dt_sims_withSignif,
                    annot1[, ..cols],
                    by.x = name1,
                    by.y = "Nhood_center")
  cols <- c("Nhood_center", cols_label[2])
  dt_edges <- merge(
    dt_edges,
    annot2[, ..cols],
    by.x = name2,
    by.y = "Nhood_center",
    suffixes = c(paste0("_", name1), paste0("_", name2))
  )
  # Use consistent naming
  if(cols_label[1] != cols_label[2]){
    setnames(
      dt_edges,
      cols_label,
      c(paste0(cols_label[1], "_", name1), paste0(cols_label[2], "_", name2))
    )
  }

  # Count the number of significant edges linking 2 labels
  dt_n_signifEdges <- as.data.table(table(
    dt_edges[is_significant == T, get(paste0(cols_label[1], "_", name1))],
    dt_edges[is_significant == T, get(paste0(cols_label[2], "_", name2))]),
    keep.rownames = T)
  setnames(dt_n_signifEdges, c(
    paste0(cols_label[1], "_", name1),
    paste0(cols_label[2], "_", name2),
    "n_signif_edges"
  ))

  # Count the number of all edges linking 2 labels
  dt_n_allEdges <- as.data.table(table(
    dt_edges[, get(paste0(cols_label[1], "_", name1))],
    dt_edges[, get(paste0(cols_label[2], "_", name2))]),
    keep.rownames = T)
  setnames(dt_n_allEdges, c(
    paste0(cols_label[1], "_", name1),
    paste0(cols_label[2], "_", name2),
    "n_all_edges"
  ))

  dt_n_edges <- merge(
    dt_n_allEdges,
    dt_n_signifEdges,
    by = c(
      paste0(cols_label[1], "_", name1),
      paste0(cols_label[2], "_", name2)
    ),
    all.x = TRUE
  )
  dt_n_edges[, ratio_retained := n_signif_edges / n_all_edges]

  return(dt_n_edges)
}
