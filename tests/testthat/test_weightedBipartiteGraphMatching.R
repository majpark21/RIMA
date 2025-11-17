test_that(
  "test that .weightedBipartiteGraphMatching returns data.table with 2 columns where each row represent a pair of matched neighbourhoods.",
  {
    milo1 <- make_test_milo()
    milo2 <- make_test_milo()

    milos <- preprocess_milos(milo1, milo2)
    dt_sims <- calculate_similarities(milos)
    dt_sims_withSignif <- calculate_nhoodnhood_significance(
      milos,
      dt_sims,
      n_scrambles = 5,
      col_scramble_label = "Sample",
      alpha = 0.05,
      direction = "b"
    )

    dt_significant_nhoodnhood <- dt_sims_withSignif[is_significant==TRUE]

    out_nosim <- .weightedBipartiteGraphMatching(dt_significant_nhoodnhood,
                                                col_sim = "sim",
                                                col_nhoods1 = NULL,
                                                col_nhoods2 = NULL,
                                                return_sim = FALSE)

    out_wsim <- .weightedBipartiteGraphMatching(dt_significant_nhoodnhood,
                                                 col_sim = "sim",
                                                 col_nhoods1 = "nhoods1",
                                                 col_nhoods2 = "nhoods2",
                                                 return_sim = TRUE)
    l_outs <- list(out_nosim, out_wsim)

    expect_all_true(unlist(lapply(l_outs, inherits, what="data.table")))
    expect_equal(ncol(out_nosim), 2)
    expect_equal(ncol(out_wsim), 3)

    # Match the same dataset on itself, so must be an exact match with similarity 1
    expect_equal(nrow(out_nosim), ncol(nhoods(milo1)))
    expect_equal(nrow(out_wsim), ncol(nhoods(milo1)))
    expect_equal(out_nosim[["nhoods1"]], out_nosim[["nhoods2"]])
    expect_equal(out_wsim[["nhoods1"]], out_wsim[["nhoods2"]])
    expect_all_equal(out_wsim[["sim"]], 1)
  }
)
