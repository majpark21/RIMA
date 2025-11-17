test_that(
  "annotate nhoods returns a data.table with the neighbourhoods' annotations.",
  {
    milo <- make_test_milo()
    cols_annot <- c("Sample", "Replicate")
    dt_annot <- annotate_nhoods(milo, cols_annot = cols_annot)

    fixed_cols <- c("Nhood", "Nhood_center")
    var_cols <- c(cols_annot, paste0(cols_annot, "_fraction"))

    expect_is(dt_annot, "data.table")
    expect_equal(nrow(dt_annot), ncol(nhoods(milo)))
    expect_equal(ncol(dt_annot), length(cols_annot)*2 + length(fixed_cols))

    expect_true(all(colnames(nhoods(milo)) %in% dt_annot[["Nhood_center"]]))
    expect_true(all(var_cols %in% colnames(dt_annot)))
  }
)
