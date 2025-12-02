#' Preprocess Milo Objects for RIMA Analysis
#'
#' Prepares a pair of Milo objects for RIMA analysis by ensuring consistent feature sets,
#' calculating neighbourhood expression, and optionally mapping features across datasets.
#'
#' @param milo1,milo2 Milo objects or paths to RDS files containing Milo objects.
#'   Their \code{nhoods} slot must be populated.
#'   Rownames (feature/gene names) must be specified to ensure correct feature comparison.
#'   Feature mapping between objects can be provided via \code{dt_features}.
#' @param assay Character string specifying which assay to use for calculating
#'   neighbourhood expression. Default is \code{"logcounts"}.
#' @param calculate_expression Logical. If \code{TRUE}, calculates neighbourhood expression
#'   from the specified assay. Required for downstream RIMA functions. Default is \code{TRUE}.
#' @param dt_features A data.frame or data.table with 2 columns mapping features across the
#'   two Milos when features are not fully overlapping.
#'   The 1st column contains feature names in \code{milo1} (resp. 2nd column for \code{milo2}).
#'   Useful for cross-species comparisons using orthologs, or for analyzing a feature subset.
#'   If \code{NULL}, rownames must be identical between both Milos.
#'
#' @returns A list of 2 Milo objects with the same rownames and calculated neighbourhood expression
#'   (if \code{calculate_expression = TRUE}).
#'
#' @export
#'
#' @examples
#' # Not run: milos <- preprocess_milos(milo1, milo2)
preprocess_milos <- function(milo1,
                             milo2,
                             assay = "logcounts",
                             calculate_expression = TRUE,
                             dt_features = NULL) {
  # require(miloR)
  # require(SingleCellExperiment)
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


  # Calculate nhoods expression -----------------
  if (calculate_expression) {
    cat("Calculate nhoods expression...\n")
    if (all(dim(milo1@nhoodExpression) == c(1, 1))) {
      milo1 <- calcNhoodExpression(milo1, assay, subset.row = NULL)
    }
    if (all(dim(milo2@nhoodExpression) == c(1, 1))) {
      milo2 <- calcNhoodExpression(milo2, assay, subset.row = NULL)
    }
  }

  if (is.null(colnames(nhoodExpression(milo1)))) {
    colnames(nhoodExpression(milo1)) <- paste0("milo1_", 1:ncol(nhoodExpression(milo1)))
  }
  if (is.null(colnames(nhoodExpression(milo2)))) {
    colnames(nhoodExpression(milo2)) <- paste0("milo2_", 1:ncol(nhoodExpression(milo2)))
  }

  return(list(milo1, milo2))
}
