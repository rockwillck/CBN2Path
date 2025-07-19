test_that("CT-CBN works", {
  posetBC = read_poset(get_examples()[1])
  patternBC = read_pattern(get_examples()[1])
  bc = Spock$new(poset=posetBC$sets, numMutations=posetBC$mutations, genotypeMatrix=patternBC)
  expect_equal(ctcbn_single(bc)[[1]]$poset$mutations, 10)
  expect_equal(ctcbn_single(bc, epsilon = 0.7)[[1]]$summary[[14]], 1.350237879)
  
  posetHIV = read_poset(get_examples()[5])
  patternHIV = read_pattern(get_examples()[5])
  hiv = Spock$new(poset=posetHIV$sets, numMutations=posetHIV$mutations, genotypeMatrix=patternHIV)
  expect_equal(ctcbn_single(hiv)[[1]]$lambda[3,1], 0.189185)
  expect_equal(ctcbn_single(hiv)[[1]]$summary[[4]], -213.0950000)
  
  posetProstate = read_poset(get_examples()[5])
  patternProstate = read_pattern(get_examples()[5])
  prostate = Spock$new(poset=posetProstate$sets, numMutations=posetProstate$mutations, genotypeMatrix=patternProstate)
  expect_equal(ctcbn_single(prostate)[[1]]$poset$mutations, 9)
  expect_equal(ctcbn_single(prostate)[[1]][[1]]$sets[1,1], NA)
})
