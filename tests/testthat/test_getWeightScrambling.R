test_that(
  "Weights for scrambling is inversely proportional to the frequency of labels.",
  {
    milo <- make_test_milo()

    weights <- .scramble_calculateWeights(milo, col_scramble = "Sample")
    expect_length(weights, ncol(milo))
    expect_equal(sum(unique(weights)), 1, tolerance = 0.01)
  }
)
