#' @export

tune_lamb = function(X, K, seqlamb, initial=TRUE, vec_x=NULL){
  n = length(X)
  p = ncol(vec_x)
  dimen = dim(X[[1]])
  nvars = prod(dimen)
  mincrit = 10^8

  for (ilambda in 1:length(seqlamb)){
    obj=DEEM(X, K, lambda=seqlamb[ilambda], initial=initial, vec_x=vec_x, niter=50)
    Sigma=kronecker(kronecker(obj$sigma[[3]],obj$sigma[[2]]),obj$sigma[[1]])
    invSigma=ginv(Sigma)
    crit=0

    essigma1 = obj$sigma[[1]]
    essigma2 = obj$sigma[[2]]
    essigma3 = obj$sigma[[3]]
    Sigma = kronecker(kronecker(obj$sigma[[3]],obj$sigma[[2]]),obj$sigma[[1]])
    invSigma = ginv(Sigma)

    for (i in 1:n){
      tmp = 0
      logf1 = -0.5*matrix(X[[i]]-obj$mu[[1]],nrow=1) %*% invSigma %*% matrix(X[[i]]-obj$mu[[1]],ncol=1)-
        0.5*(nvars*log(2*pi)+dimen[1]*dimen[2]*log(det(essigma3))
             +dimen[2]*dimen[3]*log(det(essigma1))
             +dimen[1]*dimen[3]*log(det(essigma2)))
      for (j in 1:K){
        fkoverf1 = exp(-0.5*matrix(X[[i]]-obj$mu[[j]],nrow=1)%*%invSigma%*%matrix(X[[i]]-obj$mu[[j]],ncol=1)
                       +0.5*matrix(X[[i]]-obj$mu[[1]],nrow=1)%*%invSigma%*%matrix(X[[i]]-obj$mu[[1]],ncol=1))
        tmp = tmp+obj$pi[j]*fkoverf1
      }
      crit = crit+log(tmp)+logf1
    }
    crit=-2*crit + log(n)*sum(obj$beta!=0)

    if (crit < mincrit){
      mincrit = crit
      opt_lamb = seqlamb[ilambda]
      opt_y = obj$y
    }
  }

  return(list(opt_lamb=opt_lamb,opt_bic=mincrit,opt_y=opt_y))
}



