library(CBN2Path)
pos = list()
for (i in 1:5) {
  pos[[i]] = read_poset(paste("~/Downloads/Ground_Truth/Posets/poset",i,sep=""))$sets
}
testPoset = Spock$new(
  poset = pos,
  numMutations = read_poset("~/Downloads/Ground_Truth/Posets/poset1")$mutations,
  genotypeMatrix = read_pattern("~/Downloads/Ground_Truth/genotype_115_12_35")
)
testPoset_single = Spock$new(
  poset = read_poset("~/Downloads/Ground_Truth/Posets/poset1")$sets,
  numMutations = read_poset("~/Downloads/Ground_Truth/Posets/poset1")$mutations,
  genotypeMatrix = read_pattern("~/Downloads/Ground_Truth/genotype_115_12_35")
)
results = ctcbn(testPoset)
for (r in results) {
  print(r$row)
}