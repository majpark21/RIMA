#' Extract neighbourhood coordinates from reduced dimension for visualisation
#'
#' @param milo A Milo object.
#' @param dimred The name of the reduced dimension.
#' @param use_indexCell Logical. If TRUE, use the coordinates of the neighbourhood index cells to represent the neighbourhoods. If FALSE, calculate the average coordinates of the individual cells (slower).
#' @param sel_dims Numeric of length 2. Which 2 coordinates to use.
#'
#' @returns A data.table with 3 columns, indicating the nhood and its 2 coordinates in the reduced dimension.
#' @export
#'
#' @examples
extract_centroid_nhoods <- function(milo, dimred="UMAP", use_indexCell=T, sel_dims=c(1,2)){
  .nhoodCells <- function(milo){
    l_nhood_cells <- apply(nhoods(milo), 2, function(x){names(x)[x==1]})
    return(l_nhood_cells)
  }

  mat <- SingleCellExperiment::reducedDim(milo, dimred)[, sel_dims]
  if(use_indexCell){
    l_nhood_centers <- lapply(milo@nhoodIndex, function(x) colnames(milo)[x])
    names(l_nhood_centers) <- unlist(milo@nhoodIndex)
    l_nhood_centers <- lapply(l_nhood_centers, function(x) mat[x,])
  } else {
    l_nhood_cells <- .nhoodCells(milo)
    l_nhood_coords <- lapply(l_nhood_cells, function(x) mat[x,])
    l_nhood_avg <- lapply(l_nhood_coords, function(x) matrix(colMeans(x), ncol = 2))

    l_nhood_centers <- lapply(names(l_nhood_avg), function(x) rdist::cdist(l_nhood_coords[[x]], l_nhood_avg[[x]], metric = "euclidean"))
    names(l_nhood_centers) <- names(l_nhood_avg)
    l_nhood_centers <- lapply(l_nhood_centers, which.min)
    l_nhood_centers <- lapply(names(l_nhood_centers), function(x) l_nhood_coords[[x]][l_nhood_centers[[x]],])
    names(l_nhood_centers) <- names(l_nhood_avg)
  }

  l_nhood_centers <- lapply(l_nhood_centers, function(x) data.table(matrix(x, ncol=2)))
  out <- rbindlist(l_nhood_centers, idcol = "nhood")
  out[, V1 := as.numeric(V1)]
  out[, V2 := as.numeric(V2)]
  return(out)
}
