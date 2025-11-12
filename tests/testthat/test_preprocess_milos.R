test_that(
  "preprocess_milos returns a list of 2 Milo objects with identical rownames and defined nhoodExpression.",
  {
    milo1 <- make_test_milo()
    milo2 <- make_test_milo()

    milos <- preprocess_milos(milo1, milo2)
    expect_is(milos, "list")
    expect_length(milos, 2)
    expect_true(all(sapply(milos, function(x)
      inherits(x, "Milo"))))
    expect_true(all(sapply(milos, function(x)
      dim(miloR::nhoodExpression(x)) != c(1, 1))))
  }
)
