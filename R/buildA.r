#' Build a Link-path Incidence Matrix
#'
#' This function builds a link-path incidence matrix, in which feasible routes are determined using stochastic network loading based on user-specified link costs. The routes are based on a probit model, implemented via simulation.
#' @param Links For a network with N nodes, Links is an NxN matrix. The i,j entry of Links is the link number (ID) for the node connecting node i to node j. Entries should be set to NA if there is no link between the node pair in question.
#' @param Costs A 2-column matrix, the first column of which is the link number (ID) and the second the link travel cost. Costs should be strictly positive.
#' @param theta Standard deviation of link costs when apply probit-based stochastic network loading to identify routes. Can be vector of length equal to the number of links, or a scalar (which is then replicated if necessary). Default value is theta=0, so that only shortest paths between OD pairs are generated.
#' @param nsim Number of simulation runs to find routes. In each run, link costs are modified by adding normal error, and the shortest paths computed with respect to these modified costs. nsim defaults to 100.
#' @param orig Set of nodes specified as possible origins of travel. Defaults to all nodes for which there is at least one feasible route to another node.
#' @param dest Set of nodes specified as possible destinations of travel. Defaults to all nodes that are reachable from elsewere by a feasible path.
#' @loops Include paths with loops? Defaults to FALSE.
#' @return A list containing the link-path incidence matrix A, a vector O specifying the origin for each path, and a vector D specifying the destination for each path.
#' @keywords link-path incidence matrix
#' @export
#' @examples
#' Links <- matrix(0,nrow=9, ncol=9)
#' Links[1,7] <- 1
#' Links[2,7] <- 2
#' Links[1,5] <- 3
#' Links[7,8] <- 4
#' Links[2,6] <- 5
#' Links[5,8] <- 6
#' Links[8,5] <- 7
#' Links[8,6] <- 8
#' Links[6,8] <- 9
#' Links[5,3] <- 10
#' Links[8,9] <- 11
#' Links[6,4] <- 12
#' Links[9,3] <- 13
#' Links[9,4] <- 14
#' Costs <- cbind(1:14,rep(1,14))
#' buildA(Links,Costs)
#' buildA(links,Costs,theta=0.4,orig=c(1,2),dest=c(3,4))

buildA <- function(Links,Costs,theta=0,nsim=100,orig=NA,dest=NA,loops=F){
	require(e1071)
	r <- nrow(Links)
	n <- nrow(Costs)
	A <- D <- O <- numeric(0)
	X <- Links*NA
	theta <- rep(theta,length.out=n)
	if (is.na(orig[1])) orig <- 1:r
	if (is.na(dest[1])) dest <- 1:r
	if (theta[1] <= 0) nsim <- 1
	for (i in 1:nsim){
		for ( i in Costs[,1]){
			X[Links==Costs[i,1]] <- Costs[i,2] + rnorm(1,sd=theta[i])
		}
		X[X<0] <- 1e-6		
		SP <- allShortestPaths(X)
		for (o in orig){
			for (d in dest){
				if (o!=d){
					if(!is.na(SP$length[o,d])){
						p <- extractPath(SP,o,d)
						a <- rep(0,n)
						for (j in 1:(length(p)-1) ){
							link <- Links[p[j],p[j+1]]
							a[link] <- a[link]+1
						}		
						if ((!loops & (max(a) < 1.5)) | loops){
							A <- cbind(A,a)
							O <- c(O,o)
							D <- c(D,d)
						}
					}	
				}
			}
		}
	}
	indx <- duplicated(A,MARGIN=2)
	A <- A[,!indx]
	O <- O[!indx]
	D <- D[!indx]
	indx <- order(O,D)
	A <- A[,indx]
	O <- O[indx]
	D <- D[indx]
	list(A=unname(A),O=O,D=D)
}