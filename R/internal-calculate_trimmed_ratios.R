#' Calculate Proportions of Trimmed Edges
#'
#' Internal function that calculates the number of neighbourhood-neighbourhood edges trimmed
#' after significance filtering, aggregating by neighbourhood annotations.
#'
#' @param milos A list of 2 Milo objects.
#' @param dt_sims_withSignif A data.table with at least 3 columns. Each row represents a pair
#'   of neighbourhoods across the two Milos (first two columns contain neighbourhood identifiers).
#'   Must include a logical column named \"is_significant\".
#' @param cols_label Character vector of length 2 specifying column names from the Milo's
#'   \\code{colData}. Neighbourhoods will be annotated according to these columns.
#'
#' @returns A data.table with columns representing:\n
#'   The 1st and 2nd columns: neighbourhood annotation pairs.\n
#'   \\item{n_all_edges}{Total number of edges linking the annotation pair.}\n
#'   \\item{n_signif_edges}{Number of significant edges (passed threshold).}\n
#'   \\item{ratio_retained}{Fraction of edges retained (n_signif / n_all).}
#'
#' @keywords internal
#' @noRd
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
  # Use consistent naming, mimic suffixes effect if cols are different
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
