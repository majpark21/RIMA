test_that(
  "Weights for scrambling is inversely proportional to the frequency of labels.",
  {
    milo1 <- make_test_milo()

    weights <- .getWeightScrambling(milo1, col_scramble = "Sample")
    expect_length(weights, ncol(milo1))
    expect_equal(sum(unique(weights)), 1, tolerance = 0.01)
  }
)
