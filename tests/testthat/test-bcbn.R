test_that("B-CBN works", {
  test = bcbn()
  expect_equal(length(test), 100000)
})
