#' Calculate neighbourhood-neighbourhood similarities
#'
#' @param milos List of 2 Milo objects with a filled nhoodExpression slot. Alternatively, can pass a list of 2 matrices representing neighbourhood expression.
#' @param method Similarity metric. Must be one of c('pearson', 'kendall', 'spearman').
#'
#' @returns A data.table with 3 columns. The first 2 columns indicate the pair of neighbourhoods, the thrid indicate the similarity metric.
#' @export
#'
#' @examples
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
  } else if (class(milo1) == "Milo") {
    nhoods1 <- nhoodExpression(milo1)
  } else {
    stop("'milo1' must be of class Milo or a matrix.")
  }
  if (is.matrix(milo2)) {
    nhoods2 <- milo2
  } else if (class(milo2) == "Milo") {
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
