#' Probability of path choices using logit model.
#'
#' This function calculates probabilities for selection of each path on the network based on set of path disutilities (costs). The function uses a logit random utility model.
#' @param u Vector of utlities on each path
#' @param ODpair Vector indicating the OD pair serviced by each route.
#' @param theta The logit parameter. Defaults to 1.
#' @return Vector of route choice probabilities
#' @keywords path route probability
#' @examples
#' A <- matrix(c(0,1,0,0,1,0,0,0,0,1,1,0,0,0,1,0,0,0,0,1,1,0,0,0,0,0,0,1),ncol=4,byrow=T)
#' ODpair <- c(1,1,2,2)
#' Alpha <- rep(10,7)
#' Beta <- rep(2,7)
#' pow <- rep(4,7)
#' path_flow <- c(10,20,15,15)
#' x <- A%*%path_flow
#' u <- PathCost(x,A,Alpha,Beta)
#' LogitPathProb(u,ODpair,theta=0.7)
#' @export

LogitPathProb <- function(u,ODpair,theta=1){
	p1 <- exp(-theta*u)
	p.sums <- tapply(p1,ODpair,sum)
	p.sums <- rep(p.sums,unname(table(ODpair)))
	p1/p.sums
}
