#' Probability of path choices
#'
#' This function calculates probabilities for selection of each path on the network based on set of path disutilities (costs). The function uses either a logit or probit random utility model.
#' @param u Vector of utlities on each path
#' @param ODpair Vector indicating the OD pair serviced by each route (ordered by columns of the path-link incidence matrix, A).
#' @param A Path-link incidence matrix
#' @param RUM Choice of random utility model. Can be "logit" (the default) or "probit".
#' @param theta A dispersion parameter. For the logit model, a single value specifying the logit parameter. For the probit model, a vector of standard deviations for the individual link cost errors. Defaults to 1.
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
#' PathProb(u,ODpair,A,RUM="logit",theta=0.7)
#' @export

PathProb <- function(u,ODpair,A,RUM="logit",theta=1){
	if(RUM=="logit"){
		p <- LogitPathProb(u=u,ODpair=ODpair,theta=theta)		
	}
	if(RUM=="probit"){
		p <- ProbitPathProb(u=u,ODpair=ODpair,A=A,theta=theta)			
	}
	p
}