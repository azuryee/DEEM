\name{DEEM}
\alias{DEEM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Doubly-enhanced EM algorithm}

\description{Doubly-enhanced EM algorithm for tensor clustering}

\usage{
DEEM(X, nclass, niter = 100, lambda = NULL, dfmax = n, pmax = nvars, pf = rep(1, nvars),
eps = 1e-04, maxit = 1e+05, sml = 1e-06, verbose = FALSE, ceps = 0.1,
initial = TRUE, vec_x = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{Input tensor (or matrix) list of length \eqn{n}{n}, where \eqn{n}{n} is the number of observations. Each element of the list is a tensor or matrix. The order of tensor can be any positive integer not less than 2.}
  \item{nclass}{Number of clusters.}
  \item{niter}{Maximum iteration times for EM algorithm. Default value is 100.}
  \item{lambda}{A user-specified \code{lambda} value. \code{lambda} is the weight of L1 penalty and a smaller \code{lambda} allows more variables to be nonzero}

  \item{dfmax}{The maximum number of selected variables in the model. Default is the number of observations \code{n}.}

  \item{pmax}{The maximum number of potential selected variables during iteration. In middle step, the algorithm can select at most \code{pmax} variables and then shrink part of them such that the nubmer of final selected variables is less than \code{dfmax}.}

  \item{pf}{Weight of lasso penalty. Default is a vector of value \code{1} and length \code{p}, representing L1 penalty of length \eqn{p}{p}. Can be mofidied to use adaptive lasso penalty.}

  \item{eps}{Convergence threshold for coordinate descent difference between iterations. Default value is \code{1e-04}.}
  \item{maxit}{Maximum iteration times for coordinate descent for all lambda. Default value is \code{1e+05}.}
  \item{sml}{Threshold for ratio of loss function change after each iteration to old loss function value. Default value is \code{1e-06}.}
  \item{verbose}{Indicates whether print out lambda during iteration or not. Default value is \code{FALSE}.}
  \item{ceps}{Convergence threshold for cluster mean difference between iterations. Default value is \code{1}.}
  \item{initial}{Whether to nitialize algorithm with K-means clustering. Default value is \code{TRUE}.}
  \item{vec_x}{Vectorized tensor data. Default value is \code{NULL}}
}


\value{
  \item{pi}{Estimated cluster wight.}
  \item{mu}{A list of estimated cluster means.}
  \item{sigma}{A list of estimated covariance matrices.}
  \item{gamma}{A \code{n} by \code{K} matrix of estimated membership weights.}
  \item{y}{A vector of estimated labels.}
  \item{iter}{Number of iterations until convergence.}
  \item{df}{Average zero elements in beta over iterations.}
  \item{beta}{A matrix of vectorized \code{B_k}.}
}

\references{Mai, Q., Zhang, X., Pan, Y. and Deng, K. (2020) \emph{A Doubly-Enhanced EM Algorithm for Model-Based Tensor Clustering, Manuscript.}}

\author{Qing Mai, Xin Zhang, Yuqing Pan, Kai Deng}



\examples{
dimen = c(5,5,5)
nvars = prod(dimen)
K = 2
n = 100
sigma = array(list(),3)

sigma[[1]] = sigma[[2]] = sigma[[3]] = diag(5)

B2=array(0,dim=dimen)
B2[1:3,1,1]=2

y = c(rep(1,50),rep(2,50))
M = array(list(),K)
M[[1]] = array(0,dim=dimen)
M[[2]] = B2

vec_x=matrix(rnorm(n*prod(dimen)),ncol=n)
X=array(list(),n)
for (i in 1:n){
  X[[i]] = array(vec_x[,i],dim=dimen)
  X[[i]] = M[[y[i]]] + X[[i]]
}

myfit = DEEM(X, nclass=2, lambda=0.05)
}
