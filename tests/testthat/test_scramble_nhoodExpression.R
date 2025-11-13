test_that(
  "Scramble_nhoodExpression returns a list of n matrices of dimensions [n_genes, n_nhoods].",
  {
    set.seed(7)
    milo <- make_test_milo()
    milo <- calcNhoodExpression(milo)
    n_nhoods <- dim(nhoods(milo))[2]
    n_genes <- nrow(milo)

    n_scramble <- 5
    cell_weights <- .scramble_calculateWeights(milo, col_scramble = "Sample")

    out <- .scramble_nhoodExpression(
      milo,
      n = n_scramble,
      replace = FALSE,
      weight_sample = cell_weights,
      assay = "logcounts"
    )

    expect_is(out, "list")
    expect_length(out, n_scramble)
    expect_true(all(sapply(out, function(x) dim(x) == c(n_genes, n_nhoods))))
  }
)
