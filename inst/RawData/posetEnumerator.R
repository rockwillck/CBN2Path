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

Pathway_Enumerator<-function(Poset){
  library('gtools')
  d<-dim(Poset)[1]
  PERM<-permutations(4,4)
  FEASIBLE<-rep(1,24)
  if (d>1){
    for (i in 2:d){
      for (j in 1:24){
        x<-which(PERM[j,]==Poset[i,1])
        y<-which(PERM[j,]==Poset[i,2])
        if (x>y){FEASIBLE[j]<-0}
      }
    }
  }
  return(FEASIBLE)
}

######################################################################################################################################################################################################Associated_TNBC<-read.csv("/Users/rzgar/Desktop/Research/ENCORE_LIBRARY/OLD/Associated_TNBC.csv")
######################################################################################################################################################################################################Associated_TNBC<-read.csv("/Users/rzgar/Desktop/Research/ENCORE_LIBRARY/OLD/Associated_TNBC.csv")
######################################################################################################################################################################################################Associated_TNBC<-read.csv("/Users/rzgar/Desktop/Research/ENCORE_LIBRARY/OLD/Associated_TNBC.csv")

## Step 1: Enumerating all potential graphs of length 4 and checking whether they are are acyclic or not.
allMats<-expand.grid(replicate(16, 0:1, simplify = FALSE))
allG<-numeric(2^16)

for (i in 1:(2^16)){
  vec<-which(allMats[i,]==1)
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
  allG[i]<-CycleDetector(matt)
  print(i)
}

DAGs_Index<-which(allG==0) ## Obtain the index of all potential acyclic graphs (543 DAGs) 


## Step 2:  Enumerating all directed acyclic graphs (DAGs)
count<-0
DAGs<-list()
for (i in DAGs_Index){
  count<-count+1
  
  vec<-which(allMats[i,]==1)
  matt<-matrix(0,nrow=length(vec),ncol=2)
  countt<-0
  for (j in vec){
    n1<-ceiling(j/4)
    n2<-j%%4
    if (n2==0){n2<-4}
    countt<-countt+1
    matt[countt,1]<-n1
    matt[countt,2]<-n2
  }
  matt<-rbind(c(0,0),matt)
  DAGs[[count]]<-matt
  print(i)
}



### Step 3: Identify the set of feasible pathways for each DAG [This will help us to idetify the set of unique transitively closed DAGs]
pathways<-matrix(0,nrow=543,ncol=24)
for (i in 1:543){
  pathways[i,]<-Pathway_Enumerator(DAGs[[i]])
}
Upathways<-unique(pathways) ### 219 uniue pathway signatures ---> 219 unique transitively closed DAGs


### Step 4: Enumerating the set of the 219 unique transitively closed DAGs.
sumpath<-apply(pathways,1,sum)
tcDAGs<-list()
for (i in 1:219){
  mat<-c()
  for (j in 1:543){mat<-rbind(mat,Upathways[i,])}
  sel<-which(apply(abs(pathways-mat),1,sum)==0)
  tcDAGs[[i]]<-DAGs[[sel[which.min(sumpath[sel])]]]
}


### Step 5: Formatting and storing the 219 unique transitively closed DAGs as the 219 unique posets to be utilized in the R-CBN model.
Posets<-list()
Posets[[1]]<-matrix(0,0,0)
for (i in 2:219){
  d<-dim(tcDAGs[[i]])[1]
  Posets[[i]]<-matrix(tcDAGs[[i]][(2:d),],(d-1),2)
}

saveRDS(Posets, file = "~/Posets.rds")






