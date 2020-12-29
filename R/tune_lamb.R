#' @export

tune_lamb = function(X, K, seqlamb, initial=TRUE, vec_x=NULL){
  n = length(X)
  minbic = 10^8

  for (ilambda in 1:length(seqlamb)){
    obj = DEEM(X, K, lambda=seqlamb[ilambda], initial=initial, vec_x=vec_x, niter=50)
    Xm = sapply(X,as.vector)
    loglk = logMixTenGau(Xm, obj$pi, obj$gamma, obj$mu, obj$sigma)

    bic = -2*loglk + log(n)*sum(obj$beta!=0)

    if (bic < minbic){
      minbic = bic
      opt_lamb = seqlamb[ilambda]
      opt_y = obj$y
    }
  }

  return(list(opt_lamb=opt_lamb,opt_bic=minbic,opt_y=opt_y))
}



