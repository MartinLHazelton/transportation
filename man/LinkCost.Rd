% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LinkCost.r
\name{LinkCost}
\alias{LinkCost}
\title{BPR link cost function}
\usage{
LinkCost(x, Alpha, Beta, pow = 4)
}
\arguments{
\item{x}{Traffic volume on link}

\item{Alpha}{Free flow travel time parameter}

\item{Beta}{Link capacity parameter}

\item{pow}{Polynomial order of cost function. Defaults to 4.}
}
\value{
Vector of link costs
}
\description{
This function calculates the travel costs as a function of link flow using the classic separable BPR cost function. Specifically, c(x) = Alpha*(1+0.15 (x/Beta)^pow) for pow > 0; c(x) = Alpha for pow = 0.
}
\examples{
x <- 1:25
costs <- LinkCost(x,Alpha=5,Beta=20)
plot(x,costs)
}
\keyword{BPR}
\keyword{cost}
\keyword{link}
