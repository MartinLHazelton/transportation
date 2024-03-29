% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PathCost.r
\name{PathCost}
\alias{PathCost}
\title{Path costs based on BPR link cost function}
\usage{
PathCost(x, A = diag(length(x)), Alpha, Beta, pow = 4)
}
\arguments{
\item{x}{Traffic flow on all links of network}

\item{A}{link-path incidence matrix}

\item{Alpha}{Vector of free flow travel time parameters for each link}

\item{Beta}{Vector of capacity parameters for each link}

\item{pow}{Polynomial order of cost function for each link. Defaults to 4.}
}
\value{
Vector of path costs
}
\description{
This function calculates the path costs as a function of link flows over network. BPR link cost functions are assumed.
}
\examples{
A <- matrix(c(0,1,0,0,1,0,0,0,0,1,1,0,0,0,1,0,0,0,0,1,1,0,0,0,0,0,0,1),ncol=4,byrow=T)
Alpha <- rep(10,7)
Beta <- rep(2,7)
pow <- rep(4,7)
path_flow <- c(10,20,15,15)
x <- A\%*\%path_flow
PathCost(x,A,Alpha,Beta)
}
\keyword{cost}
\keyword{path}
\keyword{route}
