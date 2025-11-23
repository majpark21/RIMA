#' Plot paired expression
#'
#' @param milos List of 2 Milo objects with a filled nhoodExpression slot.
#' @param dt_match A data.table with 2 columns containing the nhood-nhood matching. The first (resp. second) column contains the names of nhoods in the 1st (resp. 2nd) Milo object.
#' @param genes A vector of gene names. The genes must be present in both Milo's rownames. If NULL, uses all available genes.
#'
#' @returns A data.table
#' @export
#'
#' @examples
get_paired_expression <- function(milos, dt_match, genes=NULL){
  if(is.null(genes)){
    genes <- rownames(milos[[1]])
  }
  nhoods1 <- nhoodExpression(milos[[1]])
  nhoods1 <- as.data.table(t(nhoods1[genes, dt_match[[1]], drop=FALSE]))
  nhoods1[, TMP__PAIRDX := 1:.N]

  nhoods2 <- nhoodExpression(milos[[2]])
  nhoods2 <- as.data.table(t(nhoods2[genes, dt_match[[2]], drop=FALSE]))
  nhoods2[, TMP__PAIRDX := 1:.N]

  out <- list()
  for(gene in genes){
    cols <- c("TMP__PAIRDX", gene)
    out[[gene]] <- merge(nhoods1[, ..cols], nhoods2[, ..cols], by="TMP__PAIRDX")
    setnames(out[[gene]], c("TMP__PAIRDX", "expression_1", "expression_2"))
  }
  out <- rbindlist(out, idcol="gene")
  out[, TMP__PAIRDX := NULL]

  return(out)
}
