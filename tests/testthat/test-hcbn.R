test_that("H-CBN works", {
  posetBC = readPoset(getExamples()[1])
  patternBC = readPattern(getExamples()[1])
  bc = Spock$new(poset=posetBC$sets, numMutations=posetBC$mutations, genotypeMatrix=patternBC)
  expect_equal(hcbn(bc)$poset$sets[3,2], 8)

  lambdaBC = readLambda(getExamples()[1])
  bcWLambda = Spock$new(poset=posetBC$sets, numMutations=posetBC$mutations, genotypeMatrix=patternBC, lambda = lambdaBC)
  expect_error(hcbn(bc, epsilon = 0.7))
  expect_equal(hcbn(bcWLambda, epsilon = 0.7)$summary[[10]], 1.424189)
})
