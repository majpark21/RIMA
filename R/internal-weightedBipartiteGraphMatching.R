#' Weighted Bipartite Graph Matching
#'
#' @param dt_significant_nhoodnhood A data.table with 3 columns. The first 2 columns indicate the pair of neighbourhoods, the third indicates the similarity value. This should only contain the significant nhood-nhood edges to be used for matching.
#' @param col_sim Column name of dt_significant_nhoodnhood containing the similarity values.
#' @param col_nhoods1 Column name of dt_significant_nhoodnhood containing the names of the 1st (resp. 2nd) sets of nhoods to be matched. If NULL, automatically set to the 1st (resp. 2nd) column name of dt_significant_nhoodnhood
#' @param return_sim Logical. If TRUE, add a 3rd column to the output that contains the similarity value of the nhood-nhood pair.
#'
#' @returns A data.table with 2 columns where each row represent a pair of matched neighbourhoods.
#'
.weightedBipartiteGraphMatching <- function(dt_significant_nhoodnhood,
                                            col_sim = "sim",
                                            col_nhoods1 = NULL,
                                            col_nhoods2 = NULL,
                                            return_sim = TRUE) {
  if (is.null(col_nhoods1)) {
    col_nhoods1 <- colnames(dt_significant_nhoodnhood)[1]
  }
  if (is.null(col_nhoods2)) {
    col_nhoods2 <- colnames(dt_significant_nhoodnhood)[2]
  }
  sel_cols <- c(col_nhoods1, col_nhoods2, col_sim)
  dt_significant_nhoodnhood <- dt_significant_nhoodnhood[, ..sel_cols]

  setnames(dt_significant_nhoodnhood, col_sim, "weight")
  original_cols <- setdiff(colnames(dt_significant_nhoodnhood), c("weight"))

  # Ensure names of nodes are unique for each atlas
  for (ii in seq_along(original_cols)) {
    dt_significant_nhoodnhood[[original_cols[ii]]] <- paste(paste0(".GROUP", ii), dt_significant_nhoodnhood[[original_cols[ii]]], sep =
                                                              "__")
  }

  # vertex type, one atlas is coded as up, the other as down
  dt_types <- list(
    up = unique(dt_significant_nhoodnhood[[1]]),
    down = unique(dt_significant_nhoodnhood[[2]])
  )
  dt_types <- rbindlist(lapply(dt_types, as.data.table, keep.rownames =
                                 T), idcol = "type")
  setnames(dt_types, "V1", "name")
  setcolorder(dt_types, "name")
  dt_types[, type := type == "up"]  # convert to Bool for matching

  g <- igraph::graph_from_data_frame(dt_significant_nhoodnhood,
                                     directed = F,
                                     vertices = dt_types)
  # NULL weights and types will read directly from graph's attributes
  match <- igraph::max_bipartite_match(g, weights = NULL, types = NULL)

  out <- as.data.table(match$matching, keep.rownames = T)
  # Remove cells for which there is no match left
  out <- out[(!is.na(V1)) & (!is.na(V2))]
  # Each match is listed twice, because it is the same in both direction. Keep only one
  out <- out[grepl(".GROUP1", out$V1), ]
  out$V1 <- sub(".GROUP1__", "", out$V1)
  out$V2 <- sub(".GROUP2__", "", out$V2)
  setnames(out, original_cols)


  if (return_sim) {
    # Restore original nhood labels
    dt_significant_nhoodnhood[[col_nhoods1]] <- sub(".GROUP1__", "", dt_significant_nhoodnhood[[col_nhoods1]])
    dt_significant_nhoodnhood[[col_nhoods2]] <- sub(".GROUP2__", "", dt_significant_nhoodnhood[[col_nhoods2]])

    out <- merge(out,
                 dt_significant_nhoodnhood,
                 by = c(col_nhoods1, col_nhoods2))
    setnames(out, "weight", col_sim)
  }

  return(out)
}
