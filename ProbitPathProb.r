#' Probability of path choices using probit model.
#'
#' This function calculates probabilities for selection of each path on the network based on set of path disutilities (costs). The function uses a probit random utility model. Calculations are exact for two route networks. In other case the probabilities are approximated by simulation.
#' @param u Vector of utlities on each path
#' @param ODpair Vector indicating the OD pair serviced by each route.
#' @param A Path-link incidence matrix
#' @param theta Vector of standard deviations for link cost errors. Defaults to 1 in each case.
#' @param Niter Number of iterations when computing probabilities by simulation. Defaults to 1000.
#' @return Vector of route choice probabilities
#' @keywords path route probability
#' @examples
#' A <- diag(2)
#' ODpair <- c(1,1)
#' u <- c(15,17)
#' LogitPathProb(u,ODpair,A,theta=1)
#' @export

ProbitPathProb <- function(u,ODpair,A,theta,Niter=1000){
	exact <- is.matrix(A) && dim(A) == 2 && all(A == diag(2))
	if (exact){
		p <- pnorm(-u[1]+u[2],sd=theta*sqrt(2))
		p <- c(p,1-p)
	}
	if (!exact){
		n.paths <- ncol(A)
		link.errors <- matrix(rnorm(Niter*nrow(A),mean=0,sd=theta),nrow=nrow(A))
		path.errors <- t(A)%*%link.errors
		path.u <- u + path.errors
		p <- numeric(n.paths)
		for (i in unique(ODpair)){
			working.u <- path.u[ODpair==i,]
			working.u <- cbind(working.u,diag(nrow(working.u)))  # Add artifical utilities to ensure every route included in tabulation of results
			p[ODpair==i] <- table(apply(working.u,2,which.min)) - rep(1,nrow(working.u)) # Includes removal of results from artificial utilities		
		}
		p <- p/Niter
	}
	p
}
