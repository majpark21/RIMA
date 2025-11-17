test_that(
  "extract_centroid_nhoods returns a data.table with the neighbourhoods' coordinates in reduced dimension.",
  {
    milo <- make_test_milo()

    dt_centroids_avg <- extract_centroid_nhoods(milo, dimred = "PCA", use_indexCell = F)
    dt_centroids_idxCell <- extract_centroid_nhoods(milo, dimred = "PCA", use_indexCell = T)
    l_dts <- list(dt_centroids_avg, dt_centroids_idxCell)

    expect_all_true(unlist(lapply(l_dts, inherits, what="data.table")))
    expect_all_true(unlist(lapply(l_dts, function(x) ncol(x)==3)))
    expect_all_true(unlist(lapply(l_dts, function(x) nrow(x)==ncol(nhoods(milo)))))

    expect_is(dt_centroids_avg[["nhood"]], "character")
    expect_is(dt_centroids_avg[["V1"]], "numeric")
    expect_is(dt_centroids_avg[["V2"]], "numeric")
    expect_is(dt_centroids_idxCell[["nhood"]], "character")
    expect_is(dt_centroids_idxCell[["V1"]], "numeric")
    expect_is(dt_centroids_idxCell[["V2"]], "numeric")
  }
)
