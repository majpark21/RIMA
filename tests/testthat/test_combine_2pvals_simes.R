test_that(
  "test a few cases for combine_2pvals_simes, make sure that combination behaves symmetrically.",
  {
    expect_equal(.combine_2pvals_simes(0.8, 0.2), 0.4, tolerance = 0.01)
    expect_equal(.combine_2pvals_simes(0.2, 0.8), 0.4, tolerance = 0.01)
    expect_equal(.combine_2pvals_simes(0, 0.99), 0)
    expect_equal(.combine_2pvals_simes(0.99, 0), 0)
    expect_equal(.combine_2pvals_simes(1, 0.001), 0.002, tolerance = 0.01)
    expect_equal(.combine_2pvals_simes(0.001, 1), 0.002, tolerance = 0.01)
  }
)
