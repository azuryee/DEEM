#' @import stats

########################################################
########   Median bootstrap to standard error   ########
########################################################
med_se = function(x,B){
  B_median = rep(0,B)
  n = length(x)
  for (i in 1:B) {
    id = sample(1:n,n,replace=T)
    B_median[i] = median(x[id])
  }
  return(sd(B_median)/sqrt(B))
}



tnorm <- function(x){
  s=sum(as.vector(x)*as.vector(x))
  return(s)
}


cluster_err = function(K, idx, id_est) {
  # K is number of clusters
  # id_est is the estimate label
  # idx is the true label
  # return error rate
  
  perK = combinat::permn(K)
  
  n = length(idx)
  
  K_pred = perK[[1]]
  id_pred = K_pred[id_est]
  for(t in 2:length(perK)) {
    K_now = perK[[t]]
    id_now = K_now[id_est]
    if(sum(id_now!=idx)<sum(id_pred!=idx)){
      K_pred = K_now
      id_pred = id_now
    }
  }
  
  id_err = sum(id_pred!=idx)/n*100
  
  return(list(cluster_err=id_err, K_pred=K_pred, id_pred=id_pred))
}


mkronecker = function(X){
  #X is a list of matrix to do kronecker product
  
  M = length(X)
  KronX = X[[M]]
  
  if(M!=1){
    for(m in (M-1):1){
      KronX = kronecker(KronX,X[[m]])
    }
  }
  return(KronX)
}



distortion <- function(x, y, K){
  n=length(y)
  muall=array(0,dim=dim(x[[1]]))
  for (i in 1:n){
    muall=muall+x[[i]]
  }
  muall=muall/n
  
  mu=array(list(),K)
  n.fit=rep(0,K)
  for (i in 1:K){
    mu[[i]]=array(0,dim=dim(x[[1]]))
  }
  SSb=0
  for (i in 1:n){
    mu[[y[i]]]=mu[[y[i]]]+x[[i]]
    n.fit[y[i]]=n.fit[y[i]]+1
    SSb=SSb+tnorm(x[[i]]-muall)
  }
  for (i in 1:K){
    mu[[i]]=mu[[i]]/n.fit[i]
  }
  SSw=0
  for (i in 1:n){
    SSw=SSw+tnorm(x[[i]]-mu[[y[i]]])
  }	
  SSb=SSb-SSw
  dist=SSw/SSb
  
}