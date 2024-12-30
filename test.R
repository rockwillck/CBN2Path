library(CBN2Path)
for (i in 1:15) {
  testPoset = Spock$new(
    poset = read_poset(paste("~/Downloads/Ground_Truth/Posets/poset",i,sep=""))$sets,
    numMutations = read_poset(paste("~/Downloads/Ground_Truth/Posets/poset",i,sep=""))$mutations,
    genotypeMatrix = read_pattern("~/Downloads/Ground_Truth/genotype_115_12_35")
  )
  print(ctcbn(testPoset)$row)
}