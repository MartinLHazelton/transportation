#' Path costs based on BPR link cost function
#'
#' This function calculates the path costs as a function of link flows over network. BPR link cost functions are assumed. 
#' @param x Traffic flow on all links of network
#' @param A link-path incidence matrix
#' @param Alpha Vector of free flow travel time parameters for each link
#' @param Beta Vector of capacity parameters for each link
#' @param pow Polynomial order of cost function for each link. Defaults to 4.
#' @return Vector of path costs
#' @keywords path route cost
#' @examples
#' A <- matrix(c(0,1,0,0,1,0,0,0,0,1,1,0,0,0,1,0,0,0,0,1,1,0,0,0,0,0,0,1),ncol=4,byrow=T)
#' Alpha <- rep(10,7)
#' Beta <- rep(2,7)
#' pow <- rep(4,7)
#' path_flow <- c(10,20,15,15)
#' x <- A%*%path_flow
#' PathCost(x,A,Alpha,Beta)
#' @export

PathCost <- function(x,A=diag(length(x)),Alpha,Beta,pow=4){
 	c(t(A)%*%LinkCost(A%*%x,Alpha,Beta,pow=pow))
}