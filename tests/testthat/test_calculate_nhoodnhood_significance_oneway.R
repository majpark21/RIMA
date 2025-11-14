test_that(
  ".calculate_nhoodnhood_significance_oneway returns a copy 'dt_sim_true', with 2 extra columns representing the pvalue and the adjusted pvalue associated with each similarity value.",
  {
    set.seed(7)
    n_scrambles <- 5
    milo1 <- make_test_milo()
    milo2 <- make_test_milo()

    milos <- preprocess_milos(milo1, milo2)
    dt_sims_true <- calculate_similarities(milos, method = "spearman")

    out <- .calculate_nhoodnhood_significance_oneway(
      milos = milos,
      dt_sims_true = dt_sims_true,
      n_scrambles = n_scrambles,
      name_intact = "nhoods1",
      name_scramble = "nhoods2",
      col_scramble_label = "Sample",
      assay = "logcounts",
      sim_method = "spearman",
      adjust = "holm",
      alpha_adjust = 0.05
    )

    # out has its key set to nhoods1, need to set same key for comparing in tests
    setkeyv(dt_sims_true, "nhoods1")

    old_columns <- colnames(dt_sims_true)
    new_columns <- setdiff(colnames(out), colnames(dt_sims_true))
    expect_true(ncol(out)==ncol(dt_sims_true)+2)
    expect_equal(new_columns, c("pval", "pval_adjusted"))
    expect_equal(dt_sims_true[, ..old_columns], out[, ..old_columns])
    expect_is(out$pval, "numeric")
    expect_is(out$pval_adjusted, "numeric")
  }
)
