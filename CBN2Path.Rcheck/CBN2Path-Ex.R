pkgname <- "CBN2Path"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('CBN2Path')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("Base2IndVec")
### * Base2IndVec

flush(stderr()); flush(stdout())

### Name: Base2IndVec
### Title: Base2IndVec
### Aliases: Base2IndVec

### ** Examples

vec<-c(0,1,0,1)
Base2IndVec(vec)



cleanEx()
nameEx("Base2Indexing")
### * Base2Indexing

flush(stderr()); flush(stdout())

### Name: Base2Indexing
### Title: Base2Indexing
### Aliases: Base2Indexing

### ** Examples

set.seed(100)
mat<-matrix(sample(c(0,1),16,replace=TRUE),4,4)
Index<-Base2Indexing(mat)



cleanEx()
nameEx("EdgeMarginalized")
### * EdgeMarginalized

flush(stderr()); flush(stdout())

### Name: EdgeMarginalized
### Title: EdgeMarginalized
### Aliases: EdgeMarginalized

### ** Examples

DAG<-matrix(c(2,2,4,1,3,3),3,2)
LAMBDA<-c(1,4,3,2.5,2)
x<-4
PathP<-PathProb_CBN(DAG, LAMBDA, x)
EdgeProb<-EdgeMarginalized(PathP,x)



cleanEx()
nameEx("GenotypeMatrix_Mutator")
### * GenotypeMatrix_Mutator

flush(stderr()); flush(stdout())

### Name: GenotypeMatrix_Mutator
### Title: GenotypeMatrix_Mutator
### Aliases: GenotypeMatrix_Mutator

### ** Examples

set.seed(100)
gMat<-matrix(sample(c(0,1),800,replace = TRUE),200,4)
gMat_mut<-GenotypeMatrix_Mutator(gMat,0.2,0.2)



cleanEx()
nameEx("Genotype_Feasibility")
### * Genotype_Feasibility

flush(stderr()); flush(stdout())

### Name: Genotype_Feasibility
### Title: Genotype_Feasibility
### Aliases: Genotype_Feasibility

### ** Examples

Geno4<-generate_matrix_genotypes(4)
DAG<-matrix(c(4,4,4,1,2,3),3,2)
x<-4
GenoF4<-Genotype_Feasibility(Geno4,DAG,x)



cleanEx()
nameEx("Jensen_Shannon_Divergence")
### * Jensen_Shannon_Divergence

flush(stderr()); flush(stdout())

### Name: Jensen_Shannon_Divergence
### Title: Jensen_Shannon_Divergence
### Aliases: Jensen_Shannon_Divergence

### ** Examples

set.seed(100)
gMat<-matrix(sample(c(0,1),800,replace = TRUE),200,4)
PathCT<-PathProb_Quartet_CTCBN(gMat)
PathH<-PathProb_Quartet_HCBN(gMat)
JSD<-Jensen_Shannon_Divergence(PathCT,PathH)



cleanEx()
nameEx("PathProb_CBN")
### * PathProb_CBN

flush(stderr()); flush(stdout())

### Name: PathProb_CBN
### Title: PathProb_CBN: quantifies pathway probabilities using the output
###   of CT-CBN or H-CBN
### Aliases: PathProb_CBN

### ** Examples

DAG<-matrix(c(2,2,4,1,3,3),3,2)
LAMBDA<-c(1,4,3,2.5,2)
x<-4
PathP<-PathProb_CBN(DAG, LAMBDA, x)



cleanEx()
nameEx("PathProb_Quartet_BCBN")
### * PathProb_Quartet_BCBN

flush(stderr()); flush(stdout())

### Name: PathProb_Quartet_BCBN
### Title: PathProb_Quartet_BCBN
### Aliases: PathProb_Quartet_BCBN

### ** Examples

set.seed(100)
gMat<-matrix(sample(c(0,1),800,replace = TRUE),200,4)
PathB<-PathProb_Quartet_BCBN(gMat)



cleanEx()
nameEx("PathProb_Quartet_CTCBN")
### * PathProb_Quartet_CTCBN

flush(stderr()); flush(stdout())

### Name: PathProb_Quartet_CTCBN
### Title: PathProb_Quartet_CTCBN
### Aliases: PathProb_Quartet_CTCBN

### ** Examples

set.seed(100)
gMat<-matrix(sample(c(0,1),800,replace = TRUE),200,4)
PathCT<-PathProb_Quartet_CTCBN(gMat)



cleanEx()
nameEx("PathProb_Quartet_HCBN")
### * PathProb_Quartet_HCBN

flush(stderr()); flush(stdout())

### Name: PathProb_Quartet_HCBN
### Title: PathProb_Quartet_HCBN
### Aliases: PathProb_Quartet_HCBN

### ** Examples

set.seed(100)
gMat<-matrix(sample(c(0,1),800,replace = TRUE),200,4)
PathH<-PathProb_Quartet_HCBN(gMat)



cleanEx()
nameEx("PathProb_Quartet_RCBN")
### * PathProb_Quartet_RCBN

flush(stderr()); flush(stdout())

### Name: PathProb_Quartet_RCBN
### Title: PathProb_Quartet_RCBN
### Aliases: PathProb_Quartet_RCBN

### ** Examples

set.seed(100)
gMat<-matrix(sample(c(0,1),800,replace = TRUE),200,4)
PathR<-PathProb_Quartet_RCBN(gMat)



cleanEx()
nameEx("PathProb_SSWM")
### * PathProb_SSWM

flush(stderr()); flush(stdout())

### Name: PathProb_SSWM
### Title: PathProb_SSWM
### Aliases: PathProb_SSWM

### ** Examples

F<-c(0,0.1,0.2,0.1,0.2,0.4,0.3,0.2,0.2,0.1,0,0.6,0.4,0.3,0.2,1)
x<-4
PathP<-PathProb_SSWM(F,x)



cleanEx()
nameEx("Path_Edge_Mapper")
### * Path_Edge_Mapper

flush(stderr()); flush(stdout())

### Name: Path_Edge_Mapper
### Title: Path_Edge_Mapper
### Aliases: Path_Edge_Mapper

### ** Examples

PEmap<-Path_Edge_Mapper(4)



cleanEx()
nameEx("Path_Normalization")
### * Path_Normalization

flush(stderr()); flush(stdout())

### Name: Path_Normalization
### Title: Path_Normalization
### Aliases: Path_Normalization

### ** Examples

DAG<-matrix(c(2,2,4,1,3,3),3,2)
LAMBDA<-c(1,4,3,2.5,2)
x<-4
PathP<-PathProb_CBN(DAG, LAMBDA, x)
PathN<-Path_Normalization(PathP, x)



cleanEx()
nameEx("Pathway_Compatibility_Quartet")
### * Pathway_Compatibility_Quartet

flush(stderr()); flush(stdout())

### Name: Pathway_Compatibility_Quartet
### Title: Pathway_Compatibility_Quartet
### Aliases: Pathway_Compatibility_Quartet

### ** Examples

set.seed(100)
gMat<-matrix(sample(c(0,1),800,replace = TRUE),200,4)
Pathway_Compatibility_Quartet(gMat)



cleanEx()
nameEx("Pathway_Feasibility")
### * Pathway_Feasibility

flush(stderr()); flush(stdout())

### Name: Pathway_Feasibility
### Title: Pathway_Feasibility
### Aliases: Pathway_Feasibility

### ** Examples

DAG<-matrix(c(4,4,4,1,2,3),3,2)
x<-4
PathF<-Pathway_Feasibility(DAG, x)



cleanEx()
nameEx("Pathway_Genotype_Compatiblility")
### * Pathway_Genotype_Compatiblility

flush(stderr()); flush(stdout())

### Name: Pathway_Genotype_Compatiblility
### Title: Pathway_Genotype_Compatiblility
### Aliases: Pathway_Genotype_Compatiblility

### ** Examples

Geno1<-c(1,0,1,0)
Geno2<-c(1,1,0,0)
Path<-c(1,2,3,4)
Pathway_Genotype_Compatiblility(Path,Geno1)
Pathway_Genotype_Compatiblility(Path,Geno2)



cleanEx()
nameEx("Pathway_Weighting_RCBN")
### * Pathway_Weighting_RCBN

flush(stderr()); flush(stdout())

### Name: Pathway_Weighting_RCBN
### Title: Pathway_Weighting_RCBN
### Aliases: Pathway_Weighting_RCBN

### ** Examples

DAG<-matrix(c(2,2,4,1,3,3),3,2)
LAMBDA<-c(1,4,3,2.5,2)
x<-4
PathP<-PathProb_CBN(DAG, LAMBDA, x)
EdgeProb<-EdgeMarginalized(PathP,x)
PEmap<-Path_Edge_Mapper(4)
W2<-Pathway_Weighting_RCBN(EdgeProb,PEmap)



cleanEx()
nameEx("Poset_Weighting_RCBN")
### * Poset_Weighting_RCBN

flush(stderr()); flush(stdout())

### Name: Poset_Weighting_RCBN
### Title: Poset_Weighting_RCBN
### Aliases: Poset_Weighting_RCBN

### ** Examples

set.seed(100)
LogLik<-runif(219)
W1<-Poset_Weighting_RCBN(LogLik)



cleanEx()
nameEx("Predictability")
### * Predictability

flush(stderr()); flush(stdout())

### Name: Predictability
### Title: Predictability
### Aliases: Predictability

### ** Examples

set.seed(100)
gMat<-matrix(sample(c(0,1),800,replace = TRUE),200,4)
PathCT<-PathProb_Quartet_CTCBN(gMat)
PathH<-PathProb_Quartet_HCBN(gMat)
PredC<-Predictability(PathCT,4)
PredH<-Predictability(PathH,4)



cleanEx()
nameEx("bcbn")
### * bcbn

flush(stderr()); flush(stdout())

### Name: bcbn
### Title: B-CBN
### Aliases: bcbn

### ** Examples

if (!require("rBCBN")) {
    install.packages(CBN2Path::getBCBNinstall(), repos = NULL, type = "source")
}
bcbn()



cleanEx()
nameEx("combinations")
### * combinations

flush(stderr()); flush(stdout())

### Name: combinations
### Title: combinations
### Aliases: combinations

### ** Examples

COMB<-combinations(10,4)



cleanEx()
nameEx("ctcbn")
### * ctcbn

flush(stderr()); flush(stdout())

### Name: ctcbn
### Title: CT-CBN
### Aliases: ctcbn

### ** Examples

example_path <- get_examples()[1]
bc <- Spock$new(
    poset = read_poset(example_path)$sets,
    numMutations = read_poset(example_path)$mutations,
    genotypeMatrix = read_pattern(example_path)
)
ctcbn(bc)



cleanEx()
nameEx("ctcbn_single")
### * ctcbn_single

flush(stderr()); flush(stdout())

### Name: ctcbn_single
### Title: CT-CBN Single Batch
### Aliases: ctcbn_single

### ** Examples

example_path <- get_examples()[1]
bc <- Spock$new(
    poset = read_poset(example_path)$sets,
    numMutations = read_poset(example_path)$mutations,
    genotypeMatrix = read_pattern(example_path)
)
ctcbn_single(bc)



cleanEx()
nameEx("generate_data")
### * generate_data

flush(stderr()); flush(stdout())

### Name: generate_data
### Title: Generate Data
### Aliases: generate_data

### ** Examples

poset <- matrix(0, 10, 10)

poset[1, 2] <- 1
poset[2, 3] <- 1
poset[3, 4] <- 1
poset[5, 4] <- 1
poset[6, 7] <- 1
poset[8, 9] <- 1
poset[8, 10] <- 1
poset[6, 9] <- 1

tr <- transitive_closure(poset)
theta <- c(0.8, 0.7, 0.6, 0.7, 0.4, 0.25, 0.6, 0.75, 0.5, 0.2)
eps <- 0.1
N <- 400

generate_data(tr, theta, eps, N)



cleanEx()
nameEx("generate_matrix_genotypes")
### * generate_matrix_genotypes

flush(stderr()); flush(stdout())

### Name: generate_matrix_genotypes
### Title: generate_matrix_genotypes
### Aliases: generate_matrix_genotypes

### ** Examples

Geno4<-generate_matrix_genotypes(4)



cleanEx()
nameEx("getBCBNinstall")
### * getBCBNinstall

flush(stderr()); flush(stdout())

### Name: getBCBNinstall
### Title: Get BCBN .tgz path
### Aliases: getBCBNinstall

### ** Examples

getBCBNinstall()



cleanEx()
nameEx("get_examples")
### * get_examples

flush(stderr()); flush(stdout())

### Name: get_examples
### Title: Get paths to examples
### Aliases: get_examples

### ** Examples

get_examples()



cleanEx()
nameEx("hcbn")
### * hcbn

flush(stderr()); flush(stdout())

### Name: hcbn
### Title: H-CBN
### Aliases: hcbn

### ** Examples

example_path <- get_examples()[1]
bc <- Spock$new(
    poset = read_poset(example_path)$sets,
    numMutations = read_poset(example_path)$mutations,
    genotypeMatrix = read_pattern(example_path)
)
hcbn(bc)
hcbn(c(bc, bc, bc))



cleanEx()
nameEx("hcbn_single")
### * hcbn_single

flush(stderr()); flush(stdout())

### Name: hcbn_single
### Title: H-CBN Single Batch
### Aliases: hcbn_single

### ** Examples

example_path <- get_examples()[1]
bc <- Spock$new(
    poset = read_poset(example_path)$sets,
    numMutations = read_poset(example_path)$mutations,
    genotypeMatrix = read_pattern(example_path)
)
hcbn_single(bc)



cleanEx()
nameEx("permutations")
### * permutations

flush(stderr()); flush(stdout())

### Name: permutations
### Title: permutations
### Aliases: permutations

### ** Examples

PERM<-permutations(4,4)



cleanEx()
nameEx("read_lambda")
### * read_lambda

flush(stderr()); flush(stdout())

### Name: read_lambda
### Title: Read a .lambda file
### Aliases: read_lambda

### ** Examples

bcPath <- get_examples()[1]
read_lambda(bcPath)



cleanEx()
nameEx("read_pattern")
### * read_pattern

flush(stderr()); flush(stdout())

### Name: read_pattern
### Title: Read a .pat file
### Aliases: read_pattern

### ** Examples

bcPath <- get_examples()[1]
read_pattern(bcPath)



cleanEx()
nameEx("read_poset")
### * read_poset

flush(stderr()); flush(stdout())

### Name: read_poset
### Title: Read a .poset file
### Aliases: read_poset

### ** Examples

bcPath <- get_examples()[1]
read_poset(bcPath)



cleanEx()
nameEx("read_time")
### * read_time

flush(stderr()); flush(stdout())

### Name: read_time
### Title: Read a .time file
### Aliases: read_time

### ** Examples

bcPath <- get_examples()[1]
read_pattern(bcPath)



cleanEx()
nameEx("transitive_closure")
### * transitive_closure

flush(stderr()); flush(stdout())

### Name: transitive_closure
### Title: Transitive Closure
### Aliases: transitive_closure

### ** Examples

poset <- matrix(0, 10, 10)

poset[1, 2] <- 1
poset[2, 3] <- 1
poset[3, 4] <- 1
poset[5, 4] <- 1
poset[6, 7] <- 1
poset[8, 9] <- 1
poset[8, 10] <- 1
poset[6, 9] <- 1
transitive_closure(poset)



cleanEx()
nameEx("visualize_cbn_model")
### * visualize_cbn_model

flush(stderr()); flush(stdout())

### Name: visualize_cbn_model
### Title: Visualize CBN Model
### Aliases: visualize_cbn_model

### ** Examples

poset <- read_poset(get_examples()[1])
visualize_cbn_model(poset)



cleanEx()
nameEx("visualize_fitness_landscape")
### * visualize_fitness_landscape

flush(stderr()); flush(stdout())

### Name: visualize_fitness_landscape
### Title: Visualize Fitness Landscape
### Aliases: visualize_fitness_landscape

### ** Examples

Genotypes <- c(
    "0000",
    "1000",
    "0100",
    "0010",
    "0001",
    "1100",
    "1010",
    "1001",
    "0110",
    "0101",
    "0011",
    "1110",
    "1101",
    "1011",
    "0111",
    "1111"
)
#
COLintensity <- c(0, rep(0.25, 4), rep(0.5, 6), rep(0.75, 4), 1)
visualize_fitness_landscape(COLintensity)



cleanEx()
nameEx("visualize_probabilities")
### * visualize_probabilities

flush(stderr()); flush(stdout())

### Name: visualize_probabilities
### Title: Visualize Pathway Probabilities
### Aliases: visualize_probabilities

### ** Examples

visualize_probabilities(c(0.05, 0.03, 0.12, 0.04, 0.02, 0, 0.05, 0.04, 0.05, 0.06, 0.04, 0.02, 0.03, 0.02, 0.05, 0.03, 0.01, 0.09, 0.06, 0.04, 0, 0.08, 0.05, 0.02))

visualize_probabilities(c(0.05, 0.03, 0.12, 0.04, 0.02, 0, 0.05, 0.04, 0.05, 0.06, 0.04, 0.02, 0.03, 0.02, 0.05, 0.03, 0.01, 0.09, 0.06, 0.04, 0, 0.08, 0.05, 0.02), geneNames = c("AAAA", "BBBB", "CCCC", "DDDD"))

mat <- matrix(c(0.1, 0.3, 0, 0.2, 0.4, 0, 0.2, 0.2, 0.1, 0, 0.2, 0.3), ncol = 2)
visualize_probabilities(mat, columnTitles = TRUE)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
