#' Calculate Neighbourhood-Neighbourhood Similarities
#'
#' Calculates pairwise similarity scores between neighbourhoods across two datasets
#' based on their expression profiles.
#'
#' @param milos List of 2 Milo objects with a filled nhoodExpression slot,
#'   or alternatively a list of 2 matrices representing neighbourhood expression.
#' @param method Character string specifying the similarity metric.
#'   Must be one of \code{"pearson"}, \code{"kendall"}, or \code{"spearman"}.
#'   Default is \code{"spearman"}.
#'
#' @returns A data.table with 3 columns:
#'   \item{nhoods1}{Neighbourhood name from the 1st Milo object.}
#'   \item{nhoods2}{Neighbourhood name from the 2nd Milo object.}
#'   \item{sim}{Similarity value between the pair of neighbourhoods.}
#'
#' @export
#'
#' @examples
#' # Not run: calculate_similarities(milos, method = "spearman")
calculate_similarities <- function(milos, method = "spearman") {
  # require(data.table)
  if (length(milos) != 2) {
    stop("'milos' must be a list of length 2.")
  }
  if (!method %in% c("pearson", "kendall", "spearman")) {
    stop("'method' must be one of: 'pearson', 'kendall', 'spearman'.")
  }

  milo1 <- milos[[1]]
  milo2 <- milos[[2]]
  if (is.matrix(milo1)) {
    nhoods1 <- milo1
  } else if (inherits(milo1, "Milo")) {
    nhoods1 <- nhoodExpression(milo1)
  } else {
    stop("'milo1' must be of class Milo or a matrix.")
  }
  if (is.matrix(milo2)) {
    nhoods2 <- milo2
  } else if (inherits(milo2, "Milo")) {
    nhoods2 <- nhoodExpression(milo2)
  } else {
    stop("'milo2' must be of class Milo or a matrix.")
  }

  out <- cor(nhoods1, nhoods2, method = method)

  out <- data.table::melt(
    data.table::as.data.table(out, keep.rownames = TRUE),
    id.vars = "rn",
    variable.name = "col",
    value.name = "sim"
  )
  data.table::setnames(out, c("rn", "col"), c("nhoods1", "nhoods2"))

  # Make sure no factors
  out[, nhoods1 := as.character(nhoods1)]
  out[, nhoods2 := as.character(nhoods2)]

  return(out)
}
