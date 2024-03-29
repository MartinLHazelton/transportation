% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rDDSProcess.r
\name{rDDSProcess}
\alias{rDDSProcess}
\title{Simulations from day-to-day stochastic process models for traffic}
\usage{
rDDSProcess(
  ODdemand,
  ODpair,
  A,
  N.days = 100,
  N.sim = 1,
  x.start = NA,
  Alpha,
  Beta,
  pow = 4,
  RUM = "logit",
  theta = 1,
  psi = 0.05,
  process = "Standard",
  variable.demand = F,
  tol = 1e-04,
  verbose = F
)
}
\arguments{
\item{ODdemand}{Vector of origin-destination (OD) travel demands}

\item{ODpair}{Vector indicating the OD pair serviced by each route (ordered by columns of the path-link incidence matrix, A).}

\item{A}{Path-link incidence matrix}

\item{N.days}{Number of days to simulate. Defaults to 100.}

\item{N.sim}{Number of parallel simulation runs. Defaults to 1.}

\item{x.start}{Initial route flow pattern. Defaults to NA, when the initial flow is set to SUE.}

\item{Alpha}{Vector of free flow travel time parameters for each link}

\item{Beta}{Vector of capacity parameters for each link}

\item{pow}{Polynomial order of cost function for each link. Defaults to 4.}

\item{RUM}{Choice of random utility model. Can be "logit" (the default) or "probit".}

\item{theta}{A dispersion parameter. For the logit model, a single value specifying the logit parameter. For the probit model, a vector of standard deviations for the individual link cost errors. Defaults to 1.}

\item{psi}{Recency parameter, controlling weight assigned to most recent path costs when updating disutility. Defaults to 0.05.}

\item{process}{Defaults to "Standard", when a time-homogeneous Markoc process is employed in which route choice decisions are based on a utility which is a convex combination of the previous disutility and route costs. Alternatives are "TIDDS1" and "TIDDS2" time-inhomogeneous processes.}

\item{tol}{Tolerance for convergence assessment for SUE. USed if \code{x.start} not specified. Defaults to 1e-4.}

\item{verbose}{Should progress of SUE algorithm be printed out? Defaults to FALSE.}

\item{prob.model}{Route choice probability model. One of "ProductMN" (product multinomial, the default), "Poisson", and "Approximate" (a quick and dirty approximation to the others).}
}
\value{
The output is a list with components x and u, containing traffic path flows and disutilies respectively. Both x and u are 3-dimensional arrays. The first dimension indexes simulation run (if multiple parallel runs are required); the second indexes network path; and the third indexes day.
}
\description{
This functions generartes sequences of traffic flow patterns drawn from one of a number of day-to-day stochastic process models.
}
\examples{
A <- matrix(c(0,1,0,0,1,0,0,0,0,1,1,0,0,0,1,0,0,0,0,1,1,0,0,0,0,0,0,1),ncol=4,byrow=T) 
Alpha <- rep(10,7)
Beta <- rep(50,7)
pow <- rep(4,7)
ODpair <- c(1,1,2,2)
ODdemand <- c(50,50)
theta <- 0.7
N.sim <- 1
N.days <- 100
rDDSProcess(ODdemand,ODpair,A,N.days=N.days,N.sim=N.sim,x.start=c(20,0),Alpha=Alpha,Beta=Beta,pow=pow,RUM="logit",theta=theta,process="TIDDS1",variable.demand=F)
rDDSProcess(ODdemand,ODpair,A,N.days=N.days,N.sim=N.sim,x.start=c(20,0),Alpha=Alpha,Beta=Beta,pow=pow,RUM="logit",theta=theta,process="TIDDS2",variable.demand=F)
rDDSProcess(ODdemand,ODpair,A,N.days=N.days,N.sim=N.sim,x.start=c(20,0),Alpha=Alpha,Beta=Beta,pow=pow,RUM="logit",theta=theta,process="Standard",variable.demand=F)
}
\keyword{Day-to-day}
\keyword{assignment}
\keyword{process}
\keyword{stochastic}
\keyword{traffic}
