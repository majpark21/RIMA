test_that(
  "calculate_nhoodnhood_significance returns a data.table based on 'dt_sim_true', with 2 extra columns representing the pvalue and the adjusted pvalue associated with each similarity value. Test the 3 directions of matching",
  {
    set.seed(7)
    n_scrambles <- 5
    milo1 <- make_test_milo()
    milo2 <- make_test_milo()
    colnames(miloR::nhoods(milo2)) <- paste0("milo2_", colnames(miloR::nhoods(milo2)))  # ensure no confusion since both milos have same nhood names

    milos <- preprocess_milos(milo1, milo2)
    dt_sims_true <- calculate_similarities(milos, method = "spearman")

    out_lr <- calculate_nhoodnhood_significance(
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

    out_rl <- calculate_nhoodnhood_significance(
      milos,
      dt_sims_true,
      col_scramble_label = "false",
      n_scrambles = n_scrambles,
      assay = "logcounts",
      sim_method = "spearman",
      adjust = "holm",
      alpha_adjust = 0.05,
      direction = 'rl'
    )

    out_b <- calculate_nhoodnhood_significance(
      milos,
      dt_sims_true,
      col_scramble_label = "false",
      n_scrambles = n_scrambles,
      assay = "logcounts",
      sim_method = "spearman",
      adjust = "holm",
      alpha_adjust = 0.05
      # direction = c('b', 'lr', 'rl')  # test if default to b
    )

    test_output <- function(out, dt_sims_true){
      old_columns <- colnames(dt_sims_true)
      new_columns <- setdiff(colnames(out), colnames(dt_sims_true))
      expect_true(ncol(out)==ncol(dt_sims_true)+3)
      expect_equal(new_columns, c("pval", "pval_adjusted", "is_significant"))
      expect_equal(dt_sims_true[, ..old_columns], out[, ..old_columns])
      expect_is(out$pval, "numeric")
      expect_is(out$pval_adjusted, "numeric")
    }

    # out has its key set to the intact milo, need to set same key for comparing in tests
    dt_sims_true_lr <- copy(dt_sims_true)
    setkeyv(dt_sims_true_lr, "nhoods1")
    test_output(out_lr, dt_sims_true_lr)
    dt_sims_true_rl <- copy(dt_sims_true)
    setcolorder(dt_sims_true_rl, "nhoods2", "nhoods1")
    setkeyv(dt_sims_true_rl, "nhoods2")
    test_output(out_rl, dt_sims_true_rl)

    # b is a concatenation of both lr and rl + combined pvals
    dt_sims_true_b <- copy(dt_sims_true)
    setkeyv(dt_sims_true_b, c("nhoods1", "nhoods2"))
    old_columns <- colnames(dt_sims_true_b)
    new_columns <- setdiff(colnames(out_b), colnames(dt_sims_true_b))

    expect_true(ncol(out_b)==ncol(dt_sims_true_b)+9)
    expect_true(all(new_columns %in% c(
      paste0("pval_", c("lr", "rl", "combined")),
      paste0("pval_adjusted_", c("lr", "rl", "combined")),
      paste0("sim_", c("lr", "rl")),
      "is_significant"
    )))
    expect_equal(dt_sims_true_b[, ..old_columns], out_b[, ..old_columns])
    expect_is(out_b$pval_combined, "numeric")
    expect_is(out_b$pval_adjusted_combined, "numeric")
  }
)
