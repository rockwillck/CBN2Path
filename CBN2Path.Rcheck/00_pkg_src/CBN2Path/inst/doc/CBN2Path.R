## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
# if (!requireNamespace("BiocManager"))
#     install.packages("BiocManager")
# BiocManager::install("CBN2Path")

## -----------------------------------------------------------------------------
library(CBN2Path)

## -----------------------------------------------------------------------------
# The poset
DAG<-matrix(c(3,3,4,4,1,2,1,2),4,2)

# The genotype matrix
set.seed(100)
Gen1<-c(rep(0,150),sample(c(0,1),25,replace=TRUE),rep(0,25))
Gen2<-c(rep(0,175),sample(c(0,1),25,replace=TRUE))
Gen3<-c(rep(0,50),sample(c(0,1),100,replace=TRUE),rep(1,50))
Gen4<-c(sample(c(0,1),100,replace=TRUE),rep(0,50),rep(1,50))
gMat<-matrix(c(Gen1,Gen2,Gen3,Gen4),200,4)
gMat<-cbind(1,gMat)

# Preparing the inputs of the ct-cbn method
bc <- Spock$new(
     poset = DAG,
     numMutations = 4,
     genotypeMatrix =gMat
)

# Running the ct-cbn model
ResultsC<-ctcbn_single(bc)

## -----------------------------------------------------------------------------
MLposet<-ResultsC[[1]]$poset$sets

## ----fig.width=4.25, fig.height=4.25------------------------------------------
visualize_cbn_model(MLposet)

## -----------------------------------------------------------------------------
MLlmbda<-ResultsC[[1]]$lambda
loglikelihood<-ResultsC[[1]]$summary[4]

## -----------------------------------------------------------------------------
temp<-gMat[,2:5]
temp_mut<-GenotypeMatrix_Mutator(temp,0.3,0.2)
gMat_mut<-cbind(1,temp_mut)

## -----------------------------------------------------------------------------
# The poset
DAG<-matrix(c(3,3,4,4,1,2,1,2),4,2)
# Preparing the inputs of the ct-cbn method
bc <- Spock$new(
     poset = DAG,
     numMutations = 4,
     genotypeMatrix =gMat_mut
)
# Running the ct-cbn model
ResultsC_mut<-ctcbn_single(bc)

## -----------------------------------------------------------------------------
MLposet_mut<-ResultsC_mut[[1]]$poset$sets
visualize_cbn_model(MLposet_mut)

## -----------------------------------------------------------------------------
example_path <- get_examples()[1]
bc <- Spock$new(
     poset = read_poset(example_path)$sets,
     numMutations = read_poset(example_path)$mutations,
     genotypeMatrix = read_pattern(example_path)
)
ResultsC2<-ctcbn_single(bc)

## -----------------------------------------------------------------------------
Posets <- readRDS(system.file("extdata", "Posets.rds", package = "CBN2Path"))

bc <- Spock$new(
     poset = Posets,
     numMutations = 4,
     genotypeMatrix =gMat
)
ResultsC3<-ctcbn(bc)

## -----------------------------------------------------------------------------
LogLik<-numeric(219)
for (i in 1:219){
  LogLik[i]<-ResultsC3[[i]]$summary[4]
}

## -----------------------------------------------------------------------------
INDX<-which.max(LogLik)
MLposet2<-Posets[[INDX]]
identical(MLposet2,MLposet)

## -----------------------------------------------------------------------------
# The poset
DAG<-matrix(c(3,3,4,4,1,2,1,2),4,2)

# Preparing the inputs of the h-cbn method
bc <- Spock$new(
     poset = DAG,
     numMutations = 4,
     genotypeMatrix =gMat
)
# Running the h-cbn model
ResultsH<-hcbn_single(bc)

## -----------------------------------------------------------------------------
MLlmbdaH<-ResultsH[[1]]$lambda
loglikelihoodH<-ResultsH[[1]]$summary[4]

## -----------------------------------------------------------------------------
# The poset
DAG<-matrix(c(3,3,4,4,1,2,1,2),4,2)
# Preparing the inputs of the h-cbn method
bc <- Spock$new(
     poset = DAG,
     numMutations = 4,
     genotypeMatrix =gMat_mut
)
# Running the h-cbn model
ResultsH_mut<-hcbn_single(bc)

## -----------------------------------------------------------------------------
example_path <- get_examples()[1]
bc <- Spock$new(
     poset = read_poset(example_path)$sets,
     numMutations = read_poset(example_path)$mutations,
     genotypeMatrix = read_pattern(example_path)
)
ResultsH2<-hcbn_single(bc)

## -----------------------------------------------------------------------------
Posets <- readRDS(system.file("extdata", "Posets.rds", package = "CBN2Path"))

bc <- Spock$new(
     poset = Posets,
     numMutations = 4,
     genotypeMatrix =gMat
)
ResultsH3<-hcbn(bc)

## -----------------------------------------------------------------------------
LogLikH<-numeric(219)
for (i in 1:219){
  LogLikH[i]<-ResultsH3[[i]]$summary[4]
}

## -----------------------------------------------------------------------------
MLposetH<-ResultsH2[[1]]$poset$sets

INDX<-which.max(LogLikH)
MLposetH2<-Posets[[INDX]]
identical(MLposetH2,MLposetH)

## ----eval=FALSE---------------------------------------------------------------
# lambdaC<-as.numeric(ResultsC2[[1]]$lambda)
# lambdaH<-as.numeric(ResultsH2[[1]]$lambda)
# dagC<-ResultsC2[[1]]$poset$sets
# dagH<-ResultsH2[[1]]$poset$sets
# 
# ProbC<-PathProb_CBN(dagC,lambdaC,10)
# ProbH<-PathProb_CBN(dagH,lambdaH,10)

## -----------------------------------------------------------------------------
gMat2<-gMat[,2:5]
gMat2_mut<-gMat_mut[,2:5]

ProbC1<-PathProb_Quartet_CTCBN(gMat2)
ProbC2<-PathProb_Quartet_CTCBN(gMat2_mut)

ProbH1<-PathProb_Quartet_HCBN(gMat2)
ProbH2<-PathProb_Quartet_HCBN(gMat2_mut)

## -----------------------------------------------------------------------------
visualize_probabilities(ProbC1)
visualize_probabilities(ProbC2)

## -----------------------------------------------------------------------------
gene_names<-c("KRAS", "TP53", "CDKN2A", "RREB1")
visualize_probabilities(ProbC2,geneNames=gene_names)

## -----------------------------------------------------------------------------
visualize_probabilities(ProbH1)
visualize_probabilities(ProbH2)

## -----------------------------------------------------------------------------
JSD_C<-Jensen_Shannon_Divergence(ProbC1,ProbC2)
JSD_H<-Jensen_Shannon_Divergence(ProbH1,ProbH2)
JSD_C
JSD_H

## -----------------------------------------------------------------------------
Pred_C1<-Predictability(ProbC1,4)
Pred_C2<-Predictability(ProbC2,4)
Pred_C1
Pred_C2
Pred_C1-Pred_C2

## -----------------------------------------------------------------------------
Pred_H1<-Predictability(ProbH1,4)
Pred_H2<-Predictability(ProbH2,4)
Pred_H1
Pred_H2
Pred_H1-Pred_H2

## -----------------------------------------------------------------------------
PathwayC1<-Pathway_Compatibility_Quartet(gMat2)
PathwayC2<-Pathway_Compatibility_Quartet(gMat2_mut)

## -----------------------------------------------------------------------------
RhoC1<-cor(PathwayC1,ProbC1,method="spearman")
RhoC2<-cor(PathwayC2,ProbC2,method="spearman")
RhoC1
RhoC2
RhoH1<-cor(PathwayC1,ProbH1,method="spearman")
RhoH2<-cor(PathwayC2,ProbH2,method="spearman")
RhoH1
RhoH2

## -----------------------------------------------------------------------------
gMat2<-gMat[,2:5]
gMat2_mut<-gMat_mut[,2:5]

ProbR1<-PathProb_Quartet_RCBN(gMat2)
ProbR2<-PathProb_Quartet_RCBN(gMat2_mut)

## -----------------------------------------------------------------------------
visualize_probabilities(ProbR1)
visualize_probabilities(ProbR2)

## -----------------------------------------------------------------------------
JSD_R<-Jensen_Shannon_Divergence(ProbR1,ProbR2)
JSD_R

## -----------------------------------------------------------------------------
Pred_R1<-Predictability(ProbR1,4)
Pred_R2<-Predictability(ProbR2,4)
Pred_R1
Pred_R2
Pred_R1-Pred_R2

## -----------------------------------------------------------------------------
RhoR1<-cor(PathwayC1,ProbR1,method="spearman")
RhoR2<-cor(PathwayC2,ProbR2,method="spearman")
RhoR1
RhoR2

## ----warning=FALSE, results='hide'--------------------------------------------
gMat2<-gMat[,2:5]
gMat2_mut<-gMat_mut[,2:5]

ProbB1<-PathProb_Quartet_BCBN(gMat2)
ProbB2<-PathProb_Quartet_BCBN(gMat2_mut)

## -----------------------------------------------------------------------------
visualize_probabilities(ProbB1)
visualize_probabilities(ProbB2)

## -----------------------------------------------------------------------------
JSD_B<-Jensen_Shannon_Divergence(ProbB1,ProbB2)
JSD_B

## -----------------------------------------------------------------------------
Pred_B1<-Predictability(ProbB1,4)
Pred_B2<-Predictability(ProbB2,4)
Pred_B1
Pred_B2
Pred_B1-Pred_B2

## -----------------------------------------------------------------------------
RhoB1<-cor(PathwayC1,ProbB1,method="spearman")
RhoB2<-cor(PathwayC2,ProbB2,method="spearman")
RhoB1
RhoB2

## -----------------------------------------------------------------------------
F<-c(0,0.1,0.2,0.1,0.2,0.4,0.3,0.2,0.2,0.1,0,0.6,0.4,0.3,0.2,1)
G<-generate_matrix_genotypes(4)

## -----------------------------------------------------------------------------
visualize_fitness_landscape(F)

## -----------------------------------------------------------------------------
Prob_W<-PathProb_SSWM(F,4)

## -----------------------------------------------------------------------------
visualize_probabilities(Prob_W)

## ----fig.width=4.25, fig.height=4.25------------------------------------------
Pred_W<-Predictability(Prob_W,4)

