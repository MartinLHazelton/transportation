#' BPR link cost function
#'
#' This function calculates the travel costs as a function of link flow using the classic separable BPR cost function. Specifically, c(x) = Alpha*(1+0.15 (x/Beta)^pow) for pow > 0; c(x) = Alpha for pow = 0.
#' @param x Traffic volume on link
#' @param Alpha Free flow travel time parameter
#' @param Beta Link capacity parameter
#' @param pow Polynomial order of cost function. Defaults to 4.
#' @return Vector of link costs
#' @keywords BPR link cost
#' @examples
#' x <- 1:25
#' costs <- LinkCost(x,Alpha=5,Beta=20)
#' plot(x,costs)
#' @export

LinkCost <- function(x,Alpha,Beta,pow=4){
	Beta[pow==0] <- 1
	Alpha*(1+(pow>0)*0.15*(x/Beta)^pow)
}