#' Generalized Stochastic User Equilibrium
#'
#' This function calculates generalized stochastic user equilibrium (GSUE) using the method of successive averages. 
#' @param ODdemand Vector of origin-destination (OD) travel demands
#' @param ODpair Vector indicating the OD pair serviced by each route (ordered by columns of the path-link incidence matrix, A).
#' @param A Path-link incidence matrix
#' @param Alpha Vector of free flow travel time parameters for each link
#' @param Beta Vector of capacity parameters for each link
#' @param pow Polynomial order of cost function for each link. Defaults to 4.
#' @param RUM Choice of random utility model. Can be "logit" (the default) or "probit".
#' @param theta A dispersion parameter. For the logit model, a single value specifying the logit parameter. For the probit model, a vector of standard deviations for the individual link cost errors. Defaults to 1.
#' @param prob.model Route choice probability model. One of "ProductMN" (product multinomial, the default), or "Poisson".
#' @param tol Tolerance for convergence assessment (measured as route mean squared difference between current flow vector and the search direction). Defaults of 1e-4.
#' @param verbose Should progress of algorithm be printed out? Defaults to FALSE.
#' @return Output is a list of SUE route flows x and corresponding path costs u.
#' @keywords GSUE
#' @examples
#' A <- matrix(c(0,1,0,0,1,0,0,0,0,1,1,0,0,0,1,0,0,0,0,1,1,0,0,0,0,0,0,1),ncol=4,byrow=T) 
#' Alpha <- rep(10,7)
#' Beta <- rep(50,7)
#' pow <- rep(4,7)
#' ODpair <- c(1,1,2,2)
#' ODdemand <- c(50,50)
#' theta <- 0.7
#' GSUE(ODdemand,ODpair,A=A,Alpha=Alpha,Beta=Beta,pow=pow,theta=theta,verbose=F)
#' @export

GSUE <- function(ODdemand,ODpair,A,Alpha,Beta,pow=4,RUM="logit",theta=1,prob.model="ProductMN",tol=1e-4,verbose=F){
  	ODdemand.rep <- rep(ODdemand,unname(table(ODpair)))
  	x <- rep(ODdemand/unname(table(ODpair)),unname(table(ODpair)))
  	gap <- 1
  	iter <- 1
  	while(gap > tol){
		u <- ExpectedPathCost(x=x,ODdemand=ODdemand,ODpair=ODpair,A=A,Alpha=Alpha,Beta=Beta,pow=pow,prob.model=prob.model)
    		y <- PathProb(u=u,ODpair=ODpair,A=A,RUM=RUM,theta=theta)*ODdemand.rep
    		gap <- sqrt(mean(((x-y)/ODdemand.rep)^2))
   		if(verbose) cat(gap,"\n")
    		if(gap > tol){
      		x <- x + (y-x)/iter
     			iter <- iter + 1
		}  
 	}
	list(x=unname(x),u=u)
}