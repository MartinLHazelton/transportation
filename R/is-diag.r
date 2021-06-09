#' Check if matrix is diagonal
#'
#' This function indicates whether the matrix supplied is diagonal.
#' @param A Matrix to be tested
#' @return Logical value indicating whether matrix is diagonal
#' @keywords diagonal matrix
#' @examples
#' A <- cbind(c(1,1),c(0,1)
#' is_matrix(A)
#' @export

is_diag <- function(A){
	if (!is.matrix(A)) return(FALSE)
  	cd <- dim(A)
  	if (cd[1] != cd[2]) return(FALSE)
	if ( all(A==diag(diag(A))) ) return(TRUE)
}
