test_that(
  ".get_nhoodnhood_pval returns a copy 'dt_sim_true', with 2 extra columns representing the pvalue and the adjusted pvalue associated with each similarity value.",
  {
    set.seed(7)
    n_scrambles <- 5
    milo1 <- make_test_milo()
    milo2 <- make_test_milo()

    milos <- preprocess_milos(milo1, milo2)
    dt_sims_true <- calculate_similarities(milos, method="spearman")

    cell_weights <- .scramble_calculateWeights(milos[[2]], col_scramble = "Sample")
    l_scrambled_assay_nhoods <- .scramble_nhoodExpression(
      milos[[2]],
      n = n_scrambles,
      replace = FALSE,
      weight_sample = cell_weights,
      assay = "logcounts"
    )

    # Calculate similarities against scrambled nhood average expression
    l_sims_scrambled <- list()
    l_matExpression <- list(nhoodExpression(milos[[1]]))
    for (ii in 1:n_scrambles) {
      l_matExpression[[2]] <- l_scrambled_assay_nhoods[[ii]]
      dt_sims_scrambled <- calculate_similarities(milos = l_matExpression, method = "spearman")

      setnames(dt_sims_scrambled, c("nhoods1", "nhoods2", "sim"))
      l_sims_scrambled[[ii]] <- dt_sims_scrambled
    }
    dt_sims_scrambled <- rbindlist(l_sims_scrambled, idcol = "iteration")

    out <- .get_nhoodnhood_pval(
      dt_sims_true,
      dt_sims_scrambled,
      alpha = 0.05,
      adjust = "holm",
      col_sim = "sim",
      col_group = "nhoods1"
    )

    expect_is(out, "data.table")
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
