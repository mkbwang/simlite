


#' Generate a binary adjacency matrix representing a graph
#'
#'@param dimension an integer representing the number of nodes
#'@param sparseness proportion of connected edges among all the pair of nodes
#'@returns a binary matrix
#'@export
graph_setup <- function(dimension, sparseness=0.1){

  mygraph <- matrix(0, nrow=dimension, ncol=dimension)
  upper_indices <- which(upper.tri(mygraph, diag=F))
  nz_indices <- sample(upper_indices, size=sparseness*length(upper_indices))
  mygraph[nz_indices] <- 1
  mygraph <- mygraph + t(mygraph)

  return(mygraph)
}



#' Generate a precision matrix
#'
#' @description
#' Generate a positive definite matrix with given sparseness of dependencies and condition number
#'
#'@param dimension an integer representing the number of nodes
#'@param sparseness proportion of connected edges among all the pair of nodes
#'@param cond target condition number, must be larger than one
#'@param lbound lower bound of absolute nonzero value for off diagonal terms
#'@param ubound upper bound of absolute nonzero value for off diagonal terms
#'
#'@importFrom RSpectra eigs_sym
#'
prec_setup <- function(dimension, sparseness=0.1, cond=50,
                       lbound=1, ubound=2){

  mygraph <- graph_setup(dimension, sparseness)
  myprec <- matrix(0, nrow=dimension, ncol=dimension)

  # first set up the off diagonal values
  nz_indices <- which(upper.tri(mygraph) & mygraph==1)
  pos_indices <- sample(nz_indices, size=0.5*length(nz_indices))
  neg_indices <- setdiff(nz_indices, pos_indices)

  myprec[pos_indices] <- runif(length(pos_indices), min=lbound, max=ubound)
  myprec[neg_indices] <- runif(length(neg_indices), min=-ubound, max=-lbound)
  myprec <- myprec + t(myprec)

  diag(myprec) <- max(abs(myprec))*2

  # calculate how to adjust the diagonal value
  evalue_max <- eigs_sym(myprec, k=1, which="LM")$values
  evalue_min <- eigs_sym(myprec, k=1, sigma=0)$values

  change <- (evalue_max - cond*evalue_min)/(cond-1)
  myprec <- myprec + diag(dimension)*change

  return(myprec)

}



#' Convert a precision matrix to a covariance/correlation matrix
#'
#' @param invsigma a precision matrix that is symmetric and positive definite
#' @param corr logical, whether to return a correlation matrix
prec2cov <- function(invsigma, corr=FALSE){

  # use the fact that the matrix is positive definite
  sigma <- chol2inv(chol(invsigma))
  if (corr){
    sigma <- cov2cor(sigma)
  }
  return(sigma)
}
