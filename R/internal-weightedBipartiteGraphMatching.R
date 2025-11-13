#' Weighted Bipartite Graph Matching
#'
#' @param dt_signif_sims A data.table with 3 columns. The first 2 columns indicate the pair of neighbourhoods, the third indicates the similarity value. This should only contain the significant nhood-nhood edges to be used for matching.
#' @param col_sim Column name of dt_signif_sims containing the similarity values.
#'
#' @returns A data.table where each row represent a pair of matched neighbourhoods.
#'
.weightedBipartiteGraphMatching <- function(dt_signif_sims, col_sim="sim"){
  setnames(dt_signif_sims, col_sim, "weight")
  original_cols <- setdiff(colnames(dt_mnn), "weight")

  # Ensure names of nodes are unique for each atlas
  for(ii in seq_along(original_cols)){
    dt_signif_sims[[original_cols[ii]]] <- paste(paste0(".GROUP", ii), dt_signif_sims[[original_cols[ii]]], sep="__")
  }

  # vertex type, one atlas is coded as up, the other as down
  dt_types <- list(up=unique(dt_signif_sims[[1]]), down=unique(dt_signif_sims[[2]]))
  dt_types <- rbindlist(lapply(dt_types, as.data.table, keep.rownames=T), idcol = "type")
  setnames(dt_types, "V1", "name")
  setcolorder(dt_types, "name")
  dt_types[, type := type=="up"]  # convert to Bool for matching

  g <- graph_from_data_frame(dt_signif_sims, directed = F, vertices = dt_types)
  # NULL weights and types will read directly from graph's attributes
  match <- max_bipartite_match(g, weights = NULL, types = NULL)

  out <- as.data.table(match$matching, keep.rownames=T)
  # Remove cells for which there is no match left
  out <- out[(!is.na(V1)) & (!is.na(V2))]
  # Each match is listed twice, because it is the same in both direction. Keep only one
  out <- out[str_detect(V1, ".GROUP1")]
  out[, V1 := str_replace(V1, ".GROUP1__", "")]
  out[, V2 := str_replace(V2, ".GROUP2__", "")]
  setnames(out, original_cols)

  return(out)
}
