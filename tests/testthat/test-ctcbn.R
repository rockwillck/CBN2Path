test_that("CT-CBN works", {
  posetBC = readPoset(getExamples()[1])
  patternBC = readPattern(getExamples()[1])
  bc = Spock$new(poset=posetBC$sets, numMutations=posetBC$mutations, genotypeMatrix=patternBC)
  expect_equal(ctcbnSingle(bc)[[1]]$poset$mutations, 10)
  expect_equal(ctcbnSingle(bc, epsilon = 0.7)[[1]]$summary[[14]], 1.350237879)
  
  posetHIV = readPoset(getExamples()[5])
  patternHIV = readPattern(getExamples()[5])
  hiv = Spock$new(poset=posetHIV$sets, numMutations=posetHIV$mutations, genotypeMatrix=patternHIV)
  expect_equal(ctcbnSingle(hiv)[[1]]$lambda[3,1], 0.189185)
  expect_equal(ctcbnSingle(hiv)[[1]]$summary[[4]], -213.0950000)
  
  posetProstate = readPoset(getExamples()[5])
  patternProstate = readPattern(getExamples()[5])
  prostate = Spock$new(poset=posetProstate$sets, numMutations=posetProstate$mutations, genotypeMatrix=patternProstate)
  expect_equal(ctcbnSingle(prostate)[[1]]$poset$mutations, 9)
  expect_equal(ctcbnSingle(prostate)[[1]][[1]]$sets[1,1], NA)
})
