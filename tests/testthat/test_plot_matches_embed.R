test_that(
  "plot_matches_embed returns a ggplot object with points and lines representing the matches.",
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
    dt_match <- RIMA::match_nhoods(dt_sims_withSignif[is_significant == TRUE])

    out <- RIMA::plot_matches_embed(milos, dt_match, dt_palette = NULL, cols_color = NULL, dimred="PCA")
    expect_is(out, "ggplot")
  }
)
