test_that(
  "test tjat match_nhoods returns a one-to-one matching in a 2-column data.table",
  {
    set.seed(7)
    n_scrambles <- 5
    milo1 <- make_test_milo()
    milo2 <- make_test_milo()
    colnames(miloR::nhoods(milo2)) <- paste0("milo2_", colnames(miloR::nhoods(milo2)))  # ensure no confusion since both milos have same nhood names

    milos <- preprocess_milos(milo1, milo2)
    dt_sims_true <- calculate_similarities(milos, method = "spearman")

    dt_nhoodnhood_sig <- calculate_nhoodnhood_significance(
      milos,
      dt_sims_true,
      col_scramble_label = "false",
      n_scrambles = n_scrambles,
      assay = "logcounts",
      sim_method = "spearman",
      adjust = "holm",
      alpha_adjust = 0.05,
      direction = 'lr'
    )
    dt_nhoodnhood_sig <- dt_nhoodnhood_sig[is_significant==TRUE]

    out <- match_nhoods(dt_nhoodnhood_sig)

    expect_length(colnames(out), 2)
    expect_equal(uniqueN(out[[1]]), nrow(out))
    expect_equal(uniqueN(out[[2]]), nrow(out))
    expect_is(out[[1]], "character")
    expect_is(out[[2]], "character")
  }
)
