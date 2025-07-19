test_that("H-CBN works", {
  posetBC = read_poset(get_examples()[1])
  patternBC = read_pattern(get_examples()[1])
  bc = Spock$new(poset=posetBC$sets, numMutations=posetBC$mutations, genotypeMatrix=patternBC)
  expect_equal(hcbn(bc)[[1]]$poset$sets[3,2], 8)
  
  lambdaBC = read_lambda(get_examples()[1])
  bcWLambda = Spock$new(poset=posetBC$sets, numMutations=posetBC$mutations, genotypeMatrix=patternBC, lambda = lambdaBC)
  expect_error(hcbn(bc, epsilon = 0.7))
  expect_equal(hcbn(bcWLambda, epsilon = 0.7)[[1]]$summary[[10]], 1.424189)
})
