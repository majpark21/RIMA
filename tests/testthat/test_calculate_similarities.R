test_that(
  "calculate similarities returns a data.table with 3 columns. The first 2 columns indicate the pair of neighbourhoods, the thrid indicate the similarity metric.",
  {
    milo1 <- make_test_milo()
    milo2 <- make_test_milo()

    milos <- preprocess_milos(milo1, milo2)
    dt_sims <- calculate_similarities(milos)

    expect_is(dt_sims, "data.table")
    expect_equal(colnames(dt_sims), c("nhoods1", "nhoods2", "sim"))
    expect_type(dt_sims$nhoods1, "character")
    expect_type(dt_sims$nhoods2, "character")
    expect_is(dt_sims$sim, "numeric")
  }
)
