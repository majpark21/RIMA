#' Preprocess a pair of Milo objects for RIMA
#'
#' @param milo1,milo2 Milo objects or path to RDS files with Milo objects. Their nhoods slot must be populated.
#' Rownames (i.e. features/gene names) must be specified to ensure that the right features
#' are compared; feature mapping between the two objects can be provided, see 'dt_features'.
#' @param assay Assay to calculate neighbourhoods' expression.
#' @param calculate_expression Whether to calculate neighbourhood's expression. Needed for downstream steps.
#' @param dt_features DataFrame with 2 columns that maps features across both milos in case the features are not fully overlapping.
#' The 1st (resp. 2nd) column refers to rownames in milo1 (resp. milo2).
#' This is useful for example to use orthologs to compare atlases of different species.
#' This can also be used to perform the matching on a subset of features.
#' If not specified, rownames must fully match between both milos.
#'
#' @returns A list of 2 Milo objects, with the same row names.
#'
#' @export
#'
#' @examples
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
