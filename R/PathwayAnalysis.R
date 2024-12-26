

#' combinations
#'
#' @param n total number of elements in the set
#' @param r subset size
#' @param v 
#' @param set 
#' @param repeats.allowed 
#'
#' @return
#' @export
#'
#' @examples
combinations<-function (n, r, v = 1:n, set = TRUE, repeats.allowed = FALSE) 
{
  if (mode(n) != "numeric" || length(n) != 1 || n < 1 || (n%%1) != 
      0) 
    stop("bad value of n")
  if (mode(r) != "numeric" || length(r) != 1 || r < 1 || (r%%1) != 
      0) 
    stop("bad value of r")
  if (!is.atomic(v) || length(v) < n) 
    stop("v is either non-atomic or too short")
  if ((r > n) & repeats.allowed == FALSE) 
    stop("r > n and repeats.allowed=FALSE")
  if (set) {
    v <- unique(sort(v))
    if (length(v) < n) 
      stop("too few different elements")
  }
  v0 <- vector(mode(v), 0)
  if (repeats.allowed) 
    sub <- function(n, r, v) {
      if (r == 0) 
        v0
      else if (r == 1) 
        matrix(v, n, 1)
      else if (n == 1) 
        matrix(v, 1, r)
      else rbind(cbind(v[1], Recall(n, r - 1, v)), Recall(n - 
                                                            1, r, v[-1]))
    }
  else sub <- function(n, r, v) {
    if (r == 0) 
      v0
    else if (r == 1) 
      matrix(v, n, 1)
    else if (r == n) 
      matrix(v, 1, n)
    else rbind(cbind(v[1], Recall(n - 1, r - 1, v[-1])), 
               Recall(n - 1, r, v[-1]))
  }
  sub(n, r, v[1:n])
}



#' permutations
#'
#' @param n total number of elements in the set
#' @param r subset size
#' @param v 
#' @param set 
#' @param repeats.allowed 
#'
#' @return
#' @export
#'
#' @examples
permutations<-function (n, r, v = 1:n, set = TRUE, repeats.allowed = FALSE) 
{
  if (mode(n) != "numeric" || length(n) != 1 || n < 1 || (n%%1) != 
      0) 
    stop("bad value of n")
  if (mode(r) != "numeric" || length(r) != 1 || r < 1 || (r%%1) != 
      0) 
    stop("bad value of r")
  if (!is.atomic(v) || length(v) < n) 
    stop("v is either non-atomic or too short")
  if ((r > n) & repeats.allowed == FALSE) 
    stop("r > n and repeats.allowed=FALSE")
  if (set) {
    v <- unique(sort(v))
    if (length(v) < n) 
      stop("too few different elements")
  }
  v0 <- vector(mode(v), 0)
  if (repeats.allowed) 
    sub <- function(n, r, v) {
      if (r == 1) 
        matrix(v, n, 1)
      else if (n == 1) 
        matrix(v, 1, r)
      else {
        inner <- Recall(n, r - 1, v)
        cbind(rep(v, rep(nrow(inner), n)), matrix(t(inner), 
                                                  ncol = ncol(inner), nrow = nrow(inner) * n, 
                                                  byrow = TRUE))
      }
    }
  else sub <- function(n, r, v) {
    if (r == 1) 
      matrix(v, n, 1)
    else if (n == 1) 
      matrix(v, 1, r)
    else {
      X <- NULL
      for (i in 1:n) X <- rbind(X, cbind(v[i], Recall(n - 
                                                        1, r - 1, v[-i])))
      X
    }
  }
  sub(n, r, v[1:n])
}



#' generate_matrix_genotypes
#'
#' @param g genotype length
#'
#' @return a genotype matrix with ncol=g and nrow=2^g
#' @export
#'
#' @examples
generate_matrix_genotypes<-function(g) 
{
  if (g > 20) 
    stop("This would generate more than one million genotypes")
  f1 <- function(n) {
    lapply(seq.int(n), function(x) combinations(n = n, r = x))
  }
  genotNums <- f1(g)
  list.of.vectors <- function(y) {
    lapply(unlist(lapply(y, function(x) {
      apply(x, 1, list)
    }), recursive = FALSE), function(m) m[[1]])
  }
  genotNums <- list.of.vectors(genotNums)
  v <- rep(0, g)
  mat <- matrix(unlist(lapply(genotNums, function(x) {
    v[x] <- 1
    return(v)
  })), ncol = g, byrow = TRUE)
  mat <- rbind(rep(0, g), mat)
  colnames(mat) <- LETTERS[1:g]
  return(mat)
}






#' Pathway_Feasibility
#'
#' @param genotypes the full set of potential binary genotypes of a given length.
#' @param DAG matrix representing the DAG of restrictions.
#' @param x the number of mutations considered.
#'
#' @return a binary vector, which indicates feasibility or infeasibility of a set of pathways
#' @export
#'
#' @examples
Pathway_Feasibility<-function(genotypes,DAG,x){
  vec<-matrix(1,nrow=(2^x),ncol=1)
  D<-dim(DAG)[1]
  if (D>0){
    for (i in 1:(2^x)){
      for (j in 1:D){
        xx<-DAG[j,1]
        yy<-DAG[j,2]
        if ((genotypes[i,yy]==1)&&(genotypes[i,xx]==0)){vec[i]<-0}## if the mutation ordering is not respected the genotype is labelled as infeasible.
      }
    }
  }
  return(vec)
}







#' PathProb_CBN: quantifies pathway probabilities using the output of CT-CBN or H-CBN 
#'
#' @param DAG matrix representing the DAG of restrictions.
#' @param LAMBDA the lambda values, which are produced by the CBN model. Note that LAMBDA is a column vector of size x+1 and the first row always equals 1.
#' @param x the number of mutations considered.
#'
#' @return vector of probabilities assigned to a set of pathways
#' @export
#'
#' @examples
PathProb_CBN<-function(DAG,LAMBDA,x){
  
  ### Step1: genotypes
  genotypes=generate_matrix_genotypes(x)## generates the genotype space[requires "OncoSimulR" package--two lines above]
  indx<-matrix(0,nrow=2^x,ncol=1)## indexing the genotypes for easier retrival
  for (k in 1:(2^x)){for (j in 1:x){indx[k,1]=indx[k,1]+2^(j-1)*genotypes[k,j]}}
  
  ### Step2: Allowed genotypes according to the DAG of restrictions (DAG)
  allowed_set<-Pathway_Feasibility(genotypes,DAG,x)
  
  ### Step3: Pathway Probabilities
  PERM<-permutations(x,x)## all x! possible permutations (mutational pathways)
  Prob<-numeric(dim(PERM)[1])## pathway probabilities
  TOT<-0 ## (the normalization factor)
  for (i1 in 1:dim(PERM)[1]){
    TEMP1<-1 #temporarily stores the probability (later needs to be normalized)
    vec<-PERM[i1,]# (i1)-th pathway
    GENO<-matrix(0,nrow=(x+1),ncol=x) # constructing the possible x+1 genotypes, each corresponding to a given step of the given mutational pathway.  
    for (j1 in 1:x){for (k1 in (j1+1):(x+1)){GENO[k1,(vec[j1])]<-1}}
    
    flag=1;
    for (j1 in 1:(x+1)){
      index=0 # index of the (j1-th) genotype in the mutational pathway.
      for (k1 in 1:x){index<-index+2^(k1-1)*GENO[j1,k1]} 
      final_index<-which(indx==index) #finding the index of the genotype in the genotype space
      if (allowed_set[final_index]==0){flag<-0}# to check whether each genotype visited is allowed. 
    }
    
    for (j1 in 1:x){
      set1<-which(GENO[j1,]==1)
      set2<-which(GENO[j1+1,]==1)
      indx_lambda<-setdiff(set2,set1)# the index of the (j1-th) genotype in the pathway to retrieve the corresponding Lambda value
      SN<-which(GENO[j1,]==0)# the set of genotypes with one additional mutation than the current genotype
      ###################################################
      SN2<-numeric(length(SN))
      for (kaka in 1:length(SN)){
        TEMP_index=0;
        for (kk in 1:x){TEMP_index<-TEMP_index+2^(kk-1)*GENO[j1,kk]}
        TEMP_index<-TEMP_index+2^(SN[kaka]-1)
        FINAL_index<-which(indx==TEMP_index);
        
        SN2[kaka]<-allowed_set[FINAL_index]
      }
      SNN<-SN[which(SN2==1)]# THE EXIT SET: the set of (allowed) genotypes with one additional mutation than the current genotype
      T<-sum(LAMBDA[(SNN+1),1]);# sum of the lambdas of the exit set [The denominator of the equation (10) in the main text]
      ###################################################
      S<-LAMBDA[(indx_lambda+1),1]# the lambda of the (j1-th mutation) [The numerator of the equation (10) in the main text]
      TEMP1<-TEMP1*(S/T)# The multiplication in the equation (10) in the main text
    }
    if (flag==0){TEMP1<-0}# If the pathway is infeasible, its probability will be zero.
    Prob[i1]<-TEMP1 #pathway probability
  }
  
  return(Prob)
}



#' PathProb_BCBN: quantifies pathway probabilities using the output of B-CBN 
#'
#' @param MAT transition probability matrix returned by B-CBN model
#'
#' @return vector of probabilities assigned to a set of pathways
#' @export
#'
#' @examples
PathProb_BCBN<-function(MAT){
  ### Step1: Enumerating all potential pathways
  x<-dim(MAT)[1]  
  PERM<-permutations(x,x)## all x! possible permutations (mutational pathways)

  ### Step2: Quantifying pathway probabilities
  Prob<-numeric(dim(PERM)[1])## pathway probabilities
  for (i in 1:dim(PERM)[1]){
    TEMP<-1
    vec<-PERM[i,]
    for (j in 1:(x-1)){
      TEMP<-TEMP*MAT[vec[j],vec[(j+1)]];
    }
    Prob[i]<-TEMP
  }
  ### Step3: Normalizing pathway probabilities
  TOT<-sum(Prob,na.rm=TRUE)
  Prob<-Prob/TOT
  ###
  return(Prob)
}



#' Title
#'
#' @param Prob1 The first (discrete) probability distribution (vector)
#' @param Prob2 The second (discrete) probability distribution (vector)
#'
#' @return Jensen Shannon Divergence between the two (discrete) probability distributions
#' @export
#'
#' @examples
Jensen_Shannon_Divergence<-function(Prob1,Prob2){
  #Prob1: the first probability distribution
  #Prob2: the second probability distribution
  D<-0
  for (i in 1:length(Prob1)){
    if (Prob1[i]>0){D<-D+Prob1[i]*log2(Prob1[i]/(0.5*Prob1[i]+0.5*Prob2[i]))}
    if (Prob2[i]>0){D<-D+Prob2[i]*log2(Prob2[i]/(0.5*Prob1[i]+0.5*Prob2[i]))}
  }
  Dv<-(D/2)
  return(Dv)
}


#' Predictability
#'
#' @param Prob Pathway probability vector  
#' @param x The length of genotype vectors
#'
#' @return Predictability 
#' @export
#'
#' @examples
Predictability<-function(Prob,x){
  TOT<-0
  for (i in 1:length(Prob)){
    if (sum(Prob[i],na.rm=TRUE)>0){TOT=TOT-Prob[i]*log(Prob[i])}
  }
  Pred=1-(TOT/log(factorial(x)))## computing the predictability
  return(Pred)
}
  



    
    
    
    






