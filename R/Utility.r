#' Path disutility updater
#'
#' This function updates the path disutility at time t+1 based on disutility and path costs at time t. The path costs are computed from pattern of link flows.
#' @param x Traffic flow on all links of network
#' @param A Link-path incidence matrix
#' @param psi Recency parameter, controlling weight assigned to most recent path costs when updating disutility. Defaults to 1.
#' @param Alpha Vector of free flow travel time parameters for each link
#' @param Beta Vector of capacity parameters for each link
#' @param pow Polynomial order of cost function for each link. Defaults to 4.
#' @param ut Vector of path disutilities prior to update
#' @return Vector of path disutilities 
#' @keywords disutility
#' @examples
#' A <- matrix(c(0,1,0,0,1,0,0,0,0,1,1,0,0,0,1,0,0,0,0,1,1,0,0,0,0,0,0,1),ncol=4,byrow=T)
#' Alpha <- rep(10,7)
#' Beta <- rep(2,7)
#' pow <- rep(4,7)
#' path_flow <- c(10,20,15,15)
#' x <- A%*%path_flow
#' ut <- c(10,10,10,10)
#' Utility(x,A,psi=0.5,Alpha,Beta,pow,ut)
#' @export

Utility <- function(x,A=diag(length(x)),psi=1,Alpha,Beta,pow=4,ut=0){
  	psi*PathCost(x,A,Alpha=Alpha,Beta=Beta,pow=pow) + (1-psi)*ut
}