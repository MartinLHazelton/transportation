#' Stochastic User Equilibrium
#'
#' This function calculates stochastic user equilibrium (SUE) using the method of successive weighted averages. Output is a list of SUE route flows x and corresponding path costs u.
#' @param ODdemand Vector of origin-destination (OD) travel demands
#' @param ODpair Vector indicating the OD pair serviced by each route (ordered by columns of the path-link incidence matrix, A).
#' @param A Path-link incidence matrix
#' @param Alpha Vector of free flow travel time parameters for each link
#' @param Beta Vector of capacity parameters for each link
#' @param pow Polynomial order of cost function for each link. Defaults to 4.
#' @param RUM Choice of random utility model. Can be "logit" (the default) or "probit".
#' @param x.ini Initialization route flows for SUE algorithm (optional). 
#' @param theta A dispersion parameter. For the logit model, a single value specifying the logit parameter. For the probit model, a vector of standard deviations for the individual link cost errors. Defaults to 1.
#' @param d Tuning parameter; d=0 corresponds to the basic method of successive averages. Defaults to 1.
#' @param tol Tolerance for convergence assessment (measured as route mean squared difference between current flow vector and the search direction). Defaults of 1e-4.
#' @param verbose Should progress of algorithm be printed out? Defaults to FALSE.
#' @return Output is a list of SUE route flows x and corresponding path costs u.
#' @keywords SUE
#' @examples
#' A <- matrix(c(0,1,0,0,1,0,0,0,0,1,1,0,0,0,1,0,0,0,0,1,1,0,0,0,0,0,0,1),ncol=4,byrow=T) 
#' Alpha <- rep(10,7)
#' Beta <- rep(50,7)
#' pow <- rep(4,7)
#' ODpair <- c(1,1,2,2)
#' ODdemand <- c(50,50)
#' theta <- 0.7
#' SUE(ODdemand,ODpair,A=A,Alpha=Alpha,Beta=Beta,pow=pow,theta=theta,verbose=F)
#' @export

SUE <- function(ODdemand,ODpair,A,Alpha,Beta,pow=4,RUM="logit",x.ini = NULL,theta=1,tol=1e-4,verbose=F){
    if (any(ODdemand <= 0)) {
        stop("All OD demands must be strictly positive")
    }
    ODdemand.rep <- rep(ODdemand, unname(table(ODpair)))
    if (is.null(x.ini)) 
        x.ini <- rep(ODdemand/unname(table(ODpair)), unname(table(ODpair)))
    x <- x.ini
    gap <- 1
    iter <- 1
    Gam <- 0
    while (gap > tol) {
        u <- PathCost(x = x, A = A, Alpha = Alpha, Beta = Beta, 
            pow = pow)
        y <- PathProb(u = u, ODpair = ODpair, A = A, RUM = RUM, 
            theta = theta) * ODdemand.rep
        gap <- sqrt(mean(((x - y)/ODdemand.rep)^2))
        if (verbose) 
            cat(gap, "\n")
	  Gam <- Gam + iter^tune.d
	  Alp <- iter^tune.d/Gam
        if (gap > tol) {
            x <- x + Alp*(y - x)
            iter <- iter + 1
        }
    }
    list(x = unname(x), u = u)
}