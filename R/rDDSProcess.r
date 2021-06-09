#' Simulations from day-to-day stochastic process models for traffic 
#'
#' This functions generartes sequences of traffic flow patterns drawn from one of a number of day-to-day stochastic process models. 
#' @param ODdemand Vector of origin-destination (OD) travel demands
#' @param ODpair Vector indicating the OD pair serviced by each route (ordered by columns of the path-link incidence matrix, A).
#' @param A Path-link incidence matrix
#' @param N.days Number of days to simulate. Defaults to 100.
#' @param N.sim Number of parallel simulation runs. Defaults to 1.
#' @param x.start Initial route flow pattern. Defaults to NA, when the initial flow is set to SUE. 
#' @param Alpha Vector of free flow travel time parameters for each link
#' @param Beta Vector of capacity parameters for each link
#' @param pow Polynomial order of cost function for each link. Defaults to 4.
#' @param RUM Choice of random utility model. Can be "logit" (the default) or "probit".
#' @param theta A dispersion parameter. For the logit model, a single value specifying the logit parameter. For the probit model, a vector of standard deviations for the individual link cost errors. Defaults to 1.
#' @param psi Recency parameter, controlling weight assigned to most recent path costs when updating disutility. Defaults to 0.05.
#' @param process Defaults to "Standard", when a time-homogeneous Markoc process is employed in which route choice decisions are based on a utility which is a convex combination of the previous disutility and route costs. Alternatives are "TIDDS1" and "TIDDS2" time-inhomogeneous processes.
#' @variable.demand If FALSE (the default), travel demand is static through time. If TRUE, then travel demand is subject to Poisson variation with mean specified by \code{ODdemand}
#' @param prob.model Route choice probability model. One of "ProductMN" (product multinomial, the default), "Poisson", and "Approximate" (a quick and dirty approximation to the others).
#' @param tol Tolerance for convergence assessment for SUE. USed if \code{x.start} not specified. Defaults to 1e-4.
#' @param verbose Should progress of SUE algorithm be printed out? Defaults to FALSE.
#' @return The output is a list with components x and u, containing traffic path flows and disutilies respectively. Both x and u are 3-dimensional arrays. The first dimension indexes simulation run (if multiple parallel runs are required); the second indexes network path; and the third indexes day.
#' @keywords Day-to-day stochastic traffic assignment process
#' @examples
#' A <- matrix(c(0,1,0,0,1,0,0,0,0,1,1,0,0,0,1,0,0,0,0,1,1,0,0,0,0,0,0,1),ncol=4,byrow=T) 
#' Alpha <- rep(10,7)
#' Beta <- rep(50,7)
#' pow <- rep(4,7)
#' ODpair <- c(1,1,2,2)
#' ODdemand <- c(50,50)
#' theta <- 0.7
#' N.sim <- 1
#' N.days <- 100
#' rDDSProcess(ODdemand,ODpair,A,N.days=N.days,N.sim=N.sim,x.start=c(20,0),Alpha=Alpha,Beta=Beta,pow=pow,RUM="logit",theta=theta,process="TIDDS1",variable.demand=F)
#' rDDSProcess(ODdemand,ODpair,A,N.days=N.days,N.sim=N.sim,x.start=c(20,0),Alpha=Alpha,Beta=Beta,pow=pow,RUM="logit",theta=theta,process="TIDDS2",variable.demand=F)
#' rDDSProcess(ODdemand,ODpair,A,N.days=N.days,N.sim=N.sim,x.start=c(20,0),Alpha=Alpha,Beta=Beta,pow=pow,RUM="logit",theta=theta,process="Standard",variable.demand=F)
#' @export

rDDSProcess <- function(ODdemand,ODpair,A,N.days=100,N.sim=1,x.start=NA,Alpha,Beta,pow=4,RUM="logit",theta=1,psi=0.05,process="Standard",variable.demand=F,tol=1e-4,verbose=F){
	psi.change <- process=="TIDDS1" | process=="TIDDS2"
	N.routes <- ncol(A)
  	N.OD <- length(ODdemand)
  	N.routes.per.OD <- unname(table(ODpair))
	if (variable.demand) mean.ODdemand <- ODdemand
  	ODdemand.rep <- rep(ODdemand,N.routes.per.OD)
  	indices <- c(0,cumsum(N.routes.per.OD))
  	x.sim <- array(0,dim=c(N.sim,N.routes,N.days))
  	u.sim <- array(0,dim=c(N.sim,N.routes,N.days))
    	if(any(is.na(x.start))){
		if (variable.demand) ODdemand <- rpois(N.OD,mean.ODdemand)
		x.start <- SUE(ODdemand=ODdemand,ODpair=ODpair,A=A,Alpha=Alpha,Beta=Beta,pow=pow,RUM=RUM,theta=theta,tol=tol,verbose=verbose)$x
    	}
	u.start <- PathCost(x.start,A=A,Alpha=Alpha,Beta=Beta,pow=pow)
	for (ii in 1:N.sim){
 		x.sim[ii,,1] <- x.start
		u.sim[ii,,1] <- u.start
		if (process=="TIDDS1") xbar <- x.start
		for (tt in 2:(N.days)){
			if (variable.demand) ODdemand <- rpois(N.OD,mean.ODdemand)
			if(psi.change) psi <- 1/tt
			if (process!="TIDDS1"){
	      		cost <- PathCost(x.sim[ii,,tt-1],A=A,Alpha=Alpha,Beta=Beta,pow=pow)   
      			u <- psi*cost + (1-psi)*u.sim[ii,,tt-1]
			}
			if (process=="TIDDS1"){
				xbar <- (1-psi)*xbar + psi*x.sim[ii,,tt-1]
				u <- PathCost(xbar,A=A,Alpha=Alpha,Beta=Beta,pow=pow) 
			}
      		pr <- PathProb(u,ODpair=ODpair,A=A,RUM=RUM,theta=theta)
      		for (i in 1:N.OD){
        			indx <- (indices[i]+1):(indices[i+1]) 
        			x.sim[ii,indx,tt] <- rmultinom(1,ODdemand[i],pr[indx])
      		}
			u.sim[ii,,tt] <- u
		}
	}
	list(x=x.sim,u=u.sim)
}
