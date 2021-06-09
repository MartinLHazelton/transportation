#' Expected path costs
#'
#' This function calculates expected path costs as a function of mean link flows under various probability models. Its primary use is when called by the function GSUE. For the product multinomial model, a normal approximation to higher order moments is applied unless the link path incidencde matrix is diagonal.
#'
#' @param x Vector of mean link flows
#' @param ODdemand Vector of origin-destination (OD) travel demands
#' @param ODpair Vector indicating the OD pair serviced by each route (ordered by columns of the path-link incidence matrix, A).
#' @param A Path-link incidence matrix
#' @param Alpha Vector of free flow travel time parameters for each link
#' @param Beta Vector of capacity parameters for each link
#' @param pow Polynomial order of cost function for each link. Defaults to 4.
#' @param theta A dispersion parameter. For the logit model, a single value specifying the logit parameter. For the probit model, a vector of standard deviations for the individual link cost errors. Defaults to 1.
#' @param prob.model Route choice probability model. One of "ProductMN" (product multinomial, the default), or "Poisson".
#' @return Vector of expected path costs
#' @keywords expected mean path costs
#' @examples
#' ExpectedPathCost()
#' @export

ExpectedPathCost <- function(x,ODdemand,ODpair,A=diag(length(x)),Alpha,Beta,pow=4,prob.model="ProductMN"){
	Beta[pow==0] <- 1
	if(prob.model=="Poisson"){
		Omega <- A%*%diag(x)%*%t(A)
		omega <- diag(Omega)
		nu <- A%*%x
		ec <- (pow==0)*Alpha + (pow==1)*Alpha*(1 + 0.15*nu/Beta) + (pow==2)*Alpha*(1 + 0.15*((nu/Beta)^2 + omega/Beta^2)) + (pow==4)*Alpha*(1 + 0.15*((nu/Beta)^4 + 6*(nu/Beta)^2*omega/Beta^2 + 3*(omega/Beta^2)^2))
		u <- c(t(A)%*%ec)
	}
	if(prob.model=="ProductMN"){
		if(!is_diag(A)){
			N.routes <- ncol(A)
 			N.OD <- length(ODdemand)
			N.routes.per.OD <- unname(table(ODpair))
			ODdemand.rep <- rep(ODdemand,N.routes.per.OD)
			pp <- x/ODdemand.rep
  			indices <- c(0,cumsum(N.routes.per.OD))
  			Vxx <- matrix(0,ncol=N.routes,nrow=N.routes)
  			for (i in 1:N.OD){
    				indx <- (indices[i]+1):(indices[i+1])
    				Vxx[indx,indx] <- ODdemand[i]*(diag(pp[indx])-pp[indx]%*%t(pp[indx]))
  			}
			Omega <- A%*%Vxx%*%t(A)
			omega <- diag(Omega)
			nu <- A%*%x
			ec <- (pow==0)*Alpha + (pow==1)*Alpha*(1 + 0.15*nu/Beta) + (pow==2)*Alpha*(1 + 0.15*((nu/Beta)^2 + omega/Beta^2)) + (pow==4)*Alpha*(1 + 0.15*((nu/Beta)^4 + 6*(nu/Beta)^2*omega/Beta^2 + 3*(omega/Beta^2)^2))
			u <- c(t(A)%*%ec)
		}
		if(is_diag(A)){
			pp <- x/ODdemand
			nu <- x
			omega2 <- x*(1-pp)
			omega3 <- omega2*(1-2*pp)
			omega4 <- omega2*(1+3*(1-2/ODdemand)*omega2)
			nu2 <- nu^2 + omega2
			nu4 <- nu^4 + 6*omega2*nu^2 + 4*omega3*nu + omega4 
    			u <- (pow==0)*Alpha + (pow==1)*Alpha*(1 + 0.15*nu/Beta) + (pow==2)*Alpha*(1 + 0.15*nu2/Beta^2) + (pow==4)*Alpha*(1 + 0.15*nu4/Beta^4)
		}
	}
	u
}