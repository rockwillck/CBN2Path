default_data <- function() {
  
  poset <- matrix(0,3,3)
  
  poset[1,2] <-1
  poset[2,3] <-1
  
  tr<-transitiveClosure(poset)
  theta <- c(0.8, 0.7, 0.6)
  eps <- 0.1
  N <- 10
  
  generateData(tr, theta, eps, N)
}

#' B-CBN
#'
#' @param data Generated data
#' @param n_samples Number of samples <def: 25000>
#' @param theta Theta <def: 0>
#' @param epsilon Epsilon <def: 0.05>
#' @param n_chains N-Chains <def: 4>
#' @param thin Thin <def: 10>
#' @param n_cores Number of parallelized cores <def: 1>
#' @param Max_L The maximum number of iteration <def: 1000>
#'
#' @return A matrix
#' @export
#'
#' @examples
#' bcbn()
bcbn <- function(data = default_data(), n_samples = 25000, theta = 0, epsilon = 0.05, n_chains = 4,thin = 10, Max_L =1000, n_cores = 1) {
  
  if (n_chains < n_cores) {
    message(paste("Number of chains was less than number of cores. Using number of chains (", n_chains, ") as thread count.", sep = ""))
  }
  registerDoMC(cores=min(n_chains, n_cores))
  
  n<-dim(data)[2]
  n_cases <- dim(data)[1]
  
  mlist<-list()
  edgelist<-list()
  l=0
  converged = 0
  converged2<-0
  repeat {
    l=l+1
    rets <- foreach( i = 1:n_chains ) %dopar% {
      print(paste("chain:",i))
      print(theta)
      if( all(theta == 0) ) {theta=as.double(runif(n))}
      
      if(length(mlist)!=0) {
        edges_in = c(t(edgelist[[i]][n_samples][[1]]))
        theta = as.double(mlist[[i]][n_samples,1:n])
        epsilon = mlist[[i]][n_samples,n+1]
      }
      else {
        edges_in = as.integer(rep(0,n*n))
      }
      
      ret<-.C("sample_full_cbn_", theta, as.integer(n), as.double(epsilon), edges_in, as.integer(n_samples), as.integer(thin), as.integer(c(t(data))), as.integer(n_cases), theta_out=as.double(rep(0,n*n_samples)), epsilon_out=as.double(rep(0,n_samples)), edges_out=as.integer(rep(0,n_samples*n*n)), log_posterior_out=as.double(rep(0,n_samples)))
    }
    
    mlist<-list()
    mcmclist<-list()
    edgelist<-list()
    
    for( i in 1:n_chains) {
      theta_m<-matrix(rets[[i]]$theta_out,ncol=n,byrow=TRUE)
      paramatrix<-cbind(theta_m,rets[[i]]$epsilon_out,rets[[i]]$log_posterior_out)
      mlist[[i]]<-paramatrix
      mcmclist[[i]]<-mcmc(paramatrix)
      sublist<-list()
      edgesum<-0
      for(k in 1:n_samples) {
        sublist[[k]]<-matrix(rets[[i]]$edges_out[((k-1)*n*n+1):(k*n*n)],ncol=n,byrow=TRUE)
        edgesum<-edgesum+matrix(rets[[i]]$edges_out[((k-1)*n*n+1):(k*n*n)],ncol=n,byrow=TRUE)
      }
      edgelist[[i]]<-sublist
      print(summary(paramatrix))
    }
    #image(edgelist[[1]][[getmode(edgelist[[1]])]],main=paste("mode of chain 1 run ",l))
    mclist<-mcmc.list(mcmclist)
    
    terr<-try(gdiag<-gelman.diag(mclist))
    if (!inherits(terr,'try-error')) {
      print(paste0("Criterion: ",max( gdiag$psrf[,1])))
      print(gdiag)
      print("##########################################")
      if ( (max( gdiag$psrf[,1] ) < 1.1)||(l>Max_L) ) {
        converged=1
        if ( converged==1 ) {break}
      }
      print(gdiag)
    }
    print(paste("finished run:",l))
  }
  
  #elist<-c(edgelist[[1]],edgelist[[2]],edgelist[[3]],edgelist[[4]])
  elist<-c()
  for (i in 1:n_chains){elist<-c(elist,edgelist[[i]])}
  #print(gdiag)
  #marginaledges<-marginal(elist)
  #print(marginaledges)
  return(elist)
}

#' Transitive Closure
#'
#' @param poset Poset matrix
#'
#' @return Poset matrix
#' @export
#'
#' @examples
#' poset <- matrix(0, 10, 10)
#'
#' poset[1, 2] <- 1
#' poset[2, 3] <- 1
#' poset[3, 4] <- 1
#' poset[5, 4] <- 1
#' poset[6, 7] <- 1
#' poset[8, 9] <- 1
#' poset[8, 10] <- 1
#' poset[6, 9] <- 1

#' transitiveClosure(poset)
transitiveClosure <- function(poset){
  ## Returns the transitive closure of a set of relations,
  ## i.e., the minimal poset containing the relations
  p <- nrow(poset)
  for (i in 2:p){
    poset <- poset + matrixPower(poset, i)
  }
  poset <- apply(poset, 2, function(x) as.numeric(x>0))
  return(poset)
}

#' Generate Data
#'
#' @param poset Poset matrix
#' @param thetas Vector of theta values
#' @param eps Epsilon
#' @param N N
#'
#' @return A matrix
#' @export
#'
#' @examples
#' poset <- matrix(0, 10, 10)
#'
#' poset[1, 2] <- 1
#' poset[2, 3] <- 1
#' poset[3, 4] <- 1
#' poset[5, 4] <- 1
#' poset[6, 7] <- 1
#' poset[8, 9] <- 1
#' poset[8, 10] <- 1
#' poset[6, 9] <- 1
#'
#' tr <- transitiveClosure(poset)
#' theta <- c(0.8, 0.7, 0.6, 0.7, 0.4, 0.25, 0.6, 0.75, 0.5, 0.2)
#' eps <- 0.1
#' N <- 400
#'
#' generateData(tr, theta, eps, N)
generateData <- function(poset, thetas, eps, N)
{
  n <- length(thetas)
  data <- matrix(0, N, n)
  parents_list <- list()
  for(j in 1:n){parents_list[[j]] <- getParents(poset, j)}
  
  g <- graph.adjacency( poset )
  topo_sorted <- topological.sort(g)
  #print(paste("after for:", Sys.time()))
  
  for(i in 1:N){
    for(j in topo_sorted)
    {if(length(parents_list[[j]]) == 0 || prod(data[i, parents_list[[j]] ] ) == 1){data[i, j] = rbinom(1, 1, thetas[j])}}
    eps_noise <- rbinom(n, 1, eps)
    data[i, ] = (data[i, ] + eps_noise) %%2
  }
  data
}
