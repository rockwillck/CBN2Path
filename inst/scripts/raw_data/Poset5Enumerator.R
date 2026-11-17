######################################################################################################################################################################################################Associated_TNBC<-read.csv("/Users/rzgar/Desktop/Research/ENCORE_LIBRARY/OLD/Associated_TNBC.csv")
### The functions 

pairGenerator<-function(vec1,vec2){
  if (vec1[2]==vec2[1]){vec<-c(vec1[1],vec2[2])}
  else if (vec1[1]==vec2[2]){vec<-c(vec2[1],vec1[2])}
  else {vec<-c()}
  return(vec)
}

### DAGcomplete function:
DAGcomplete<-function(cDAG){
  x<-dim(cDAG)[1]
  if (x>1){
    FLAG<-0
    while (FLAG==0){
      x<-dim(cDAG)[1]
      d<-choose(x,2)
      Cmat<-combn(1:x,2,simplify = TRUE)
      count<-0
      for (i in 1:d){
        vec1<-cDAG[Cmat[1,i],]
        vec2<-cDAG[Cmat[2,i],]
        vec<-pairGenerator(vec1,vec2)
        if (length(vec)>0){
          if (length(intersect(which(cDAG[,1]==vec[1]),which(cDAG[,2]==vec[2])))==0){
            cDAG<-rbind(cDAG,vec)
            count<-count+1
          }
        }
      }
      cDAG<-unique(rbind(cDAG,c(0,0)))
      cDAG<-unique(cDAG)
      if (count==0){FLAG<-1}
    }
    ww<-which(apply(cDAG,1,sum)==0)
    cDAG[ww,]<-cDAG[dim(cDAG)[1],]
    cDAG[dim(cDAG)[1],]<-c(0,0)
    return(cDAG)
  }
  else{
    return(cDAG)
  }
}


### CycleDetector:
CycleDetector<-function(cDAG){
  cDAG<-DAGcomplete(cDAG)
  x<-dim(cDAG)[1]
  cycle<-0
  if (x>1){
    d<-choose(x,2)
    Cmat<-combn(1:x,2,simplify = TRUE) 
    for (i in 1:d){
      vec1<-cDAG[Cmat[1,i],]
      vec2<-cDAG[Cmat[2,i],]
      if ((vec1[1]==vec1[[2]])&&(vec2[2]==vec2[1])){cycle<-1;break;}
      if ((vec1[1]==vec2[[2]])&&(vec1[2]==vec2[1])){cycle<-1;break;}
    }
  }
  else {if (x==1){if (cDAG[,1]==cDAG[,2]){cycle<-1}}}
  return(cycle)
}


### Pathway_Enumerator:
Pathway_Enumerator<-function(Poset,x){
  library('gtools')
  d<-dim(Poset)[1]
  PERM<-permutations(x,x)
  FEASIBLE<-rep(1,factorial(x))
  if (d>1){
    for (i in 2:d){
      for (j in 1:factorial(x)){
        xx<-which(PERM[j,]==Poset[i,1])
        yy<-which(PERM[j,]==Poset[i,2])
        if (xx>yy){FEASIBLE[j]<-0}
      }
    }
  }
  return(FEASIBLE)
}



######################################################################################################################################################################################################Associated_TNBC<-read.csv("/Users/rzgar/Desktop/Research/ENCORE_LIBRARY/OLD/Associated_TNBC.csv")
######################################################################################################################################################################################################Associated_TNBC<-read.csv("/Users/rzgar/Desktop/Research/ENCORE_LIBRARY/OLD/Associated_TNBC.csv")

## Step 1: Enumerating all potential graphs of length 4 

allMats4<-expand.grid(replicate(16, 0:1, simplify = FALSE))

DAGs4<-numeric(2^16)
for (i in 1:(2^16)){
  vec<-which(allMats4[i,]==1)
  matt<-matrix(0,nrow=length(vec),ncol=2)
  count<-0
  for (j in vec){
    n1<-ceiling(j/4)
    n2<-j%%4
    if (n2==0){n2<-4}
    count<-count+1
    matt[count,1]<-n1
    matt[count,2]<-n2
  }
  DAGs4[i]<-CycleDetector(matt)
  print(i)
}


## Step 2:  Enumerating all directed acyclic graphs (DAGs) of length 5

## Step 2.1: Enumerating all potential graphs of length 5 starting from the graphs of length 4
vecINDX<-c(rep(1:4,each=4),1:4,rep(5,5),rep(1:4,4),rep(5,4),1:5)
edgeINDX<-matrix(vecINDX,25,2)

selMats4<-allMats4[which(DAGs4==0),]
allMats3<-expand.grid(replicate(9, 0:1, simplify = FALSE))
allMats5<-matrix(0,(543*512),25)

count<-0
for (i in 1:543){
  for (j in 1:512){
    count<-count+1
    allMats5[count,]<-c(as.numeric(selMats4[i,]),as.numeric(allMats3[j,]))
    print(paste(i,j))
  }
}

## Step 2.2. Finding acyclic graphs of length 5.
DAGs5<-numeric(543*512)
for (i in 1:(543*512)){
  vec<-which(allMats5[i,]==1)
  matt<-matrix(0,nrow=length(vec),ncol=2)
  count<-0
  for (j in vec){
    count<-count+1
    matt[count,1]<-edgeINDX[j,1]
    matt[count,2]<-edgeINDX[j,2]
  }
  DAGs5[i]<-CycleDetector(matt)
  print(i)
}

### Step 3: Identify the set of feasible pathways for each DAG [This will help us to idetify the set of unique transitively closed DAGs]

#Step 3.1
count<-0
POSETS5<-list()
for (i in which(DAGs5==0)){
  count<-count+1
  vec<-which(allMats5[i,]==1)
  matt<-matrix(0,nrow=length(vec),ncol=2)
  countt<-0
  for (j in vec){
    countt<-countt+1
    matt[countt,1]<-edgeINDX[j,1]
    matt[countt,2]<-edgeINDX[j,2]
  }
  matt<-rbind(c(0,0),matt)
  POSETS5[[count]]<-matt
  print(i)
}


#Step 3.2
pathways5<-matrix(0,nrow=29281,ncol=120)
for (i in 1:29281){
  pathways5[i,]<-Pathway_Enumerator(POSETS5[[i]],5)
  print(i)
}
Upathways5<-unique(pathways5)


### Step 4: Enumerating the set of the 4231 unique transitively closed DAGs.
sumpath5<-apply(pathways5,1,sum)
tcDAGs<-list()
for (i in 1:4231){
  mat<-matrix(rep(Upathways5[i,], each = 29281), nrow = 29281)
  sel<-which(apply(abs(pathways5-mat),1,sum)==0)
  tcDAGs[[i]]<-POSETS5[[sel[which.min(sumpath5[sel])]]]
  print(i)
}

### Step 5: Formatting and storing the 4231 unique transitively closed DAGs as the 4231 unique posets to be utilized in the R-CBN model.
Posets<-list()
Posets[[1]]<-matrix(0,0,0)
for (i in 2:4231){
  d<-dim(tcDAGs[[i]])[1]
  Posets[[i]]<-matrix(tcDAGs[[i]][(2:d),],(d-1),2)
}

saveRDS(Posets, file = "~/Posets5.rds")
#The final number of unique posets of size 5 according to the above code is thus 4231, which corroborates with the number reported in The On-Line Encyclopedia of Integer Sequences (OEIS): https://oeis.org/A001035

