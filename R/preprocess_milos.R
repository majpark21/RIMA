#' Preprocess a pair of Milo objects for RIMA
#'
#' @param milo1,milo2 Milo objects or path to RDS files with Milo objects. Their nhoods slot must be populated.
#' Rownames (i.e. features/gene names) must be specified to ensure that the right features
#' are compared; feature mapping between the two objects can be provided, see 'dt_features'.
#' @param calculate_similarities
#' @param dt_features
#'
#' @returns
#' @export
#'
#' @examples
preprocess_milos <- function(milo1,
                             milo2,
                             assay = "logcounts",
                             calculate_similarities = TRUE,
                             dt_features = NULL) {
  require(miloR)
  require(SingleCellExperiment)
  # Check inputs -------------
  cat("Reading milo1...\n")
  if (is.character(milo1)) {
    milo1 <- readRDS(milo1)
  }
  if (!class(milo1) %in% "Milo") {
    stop(sprintf("milo1 must be a MiloR object, found: %s", class(milo1)))
  }
  if (all(dim(nhoods(milo1)) == c(1, 1))) {
    stop("Neighbourhoods of milo1 are not defined.")
  }

  cat("Reading milo2...\n")
  if (is.character(milo2)) {
    milo2 <- readRDS(milo2)
  }
  if (!class(milo2) %in% "Milo") {
    stop(sprintf("milo2 must be a MiloR object, found: %s"),
         class(milo2))
  }
  if (all(dim(nhoods(milo1)) == c(1, 1))) {
    stop("Neighbourhoods of milo2 are not defined.")
  }

  if (!is.null(dt_features)) {
    if (ncol(dt_features) != 2) {
      warning(
        sprintf(
          "dt_features has %s columns, only 2 needed. The 1st column (%s), and 2nd column (%s) will be used to map the genes of milo1 and milo2 respectively.",
          ncol(dt_features),
          colnames(dt_features)[1],
          colnames(dt_features)[2]
        )
      )
    }
  }


  # Subset to/reorder genes of interest -------------
  cat("Subset to/reorder genes of interest...\n")
  if (is.null(dt_features)) {
    if (any(sort(unique(rownames(milo1))) != sort(unique(rownames(milo2))))) {
      stop(
        "milo1 and milo2 must have the same rownames, otherwise dt_features must be specified."
      )
    }
    # Reorder features to same order
    milo2 <- milo2[rownames(milo1), ]
  } else {
    milo1 <- milo1[dt_features[[1]], ]
    milo2 <- milo2[dt_features[[2]], ]
    # Rename features in milo2 to their milo1's counterpart
    rownames(milo2) <- dt_features[[1]]
  }


  # Calculate nhoods similarities -----------------
  if (calculate_similarities) {
    cat("Calculate nhoods similarities...\n")
    if (all(dim(milo1@nhoodExpression) == c(1, 1))) {
      milo1 <- miloR::calcNhoodExpression(milo1, assay, subset.row = NULL)
    }
    if (all(dim(milo2@nhoodExpression) == c(1, 1))) {
      milo2 <- miloR::calcNhoodExpression(milo2, assay, subset.row = NULL)
    }
  }

  return(list(milo1, milo2))
}
