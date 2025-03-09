
# log likelihoods

#'sum of log likelihoods for a poisson distribution
#'@param param a named vector with one element "lambda"
#'@param x vector of observed counts
#'@return sum of log likelihoods
#'@importFrom stats dpois
logpois <- function(param, x){
  lambda <- param["lambda"]
  return(sum(dpois(x, lambda, log=T)))
}


#'sum of log likelihoods for a zero inflated poisson distribution
#'@param param a named vector two elements. One is the zero probability (pi0) and the other is the mean of nonzero part (lambda)
#'@param x vector of observed counts
#'@importFrom stats dpois
logzipois <- function(param, x){
  pi0 <- param["pi0"]
  lambda <- param["lambda"]

  x_nonzero <- x[x>0]
  loglik_nonzero <- log((1-pi0)*dpois(x_nonzero, lambda)) |> sum()
  loglik_zero <- sum(x==0) * log(pi0 + (1-pi0)*dpois(0, lambda))

  return(loglik_nonzero + loglik_zero)

}

#'sum of log likelihoods for a gamma distribution
#'@param param a named vector two elements. One is the shape parameter (alpha) and the other is the scale (beta)
#'@param x vector of observed counts
#'@importFrom stats dgamma
loggamma <- function(param, x){
  alpha <- param["shape"]
  beta <- param["scale"]
  return(sum(dgamma(x, shape=alpha, scale=beta, log=T)))
}


#'sum of log likelihoods for a zero inflated gamma distribution
#'@param param a named vector with three elements. One is the zero probability pi0.
#'The other two are shape parameter (alpha) and the other is the scale (beta)
#'@param x vector of observed counts
#'@importFrom stats dgamma
#'
logzigamma <- function(param, x){

  pi0 <- param['pi0']
  alpha <- param['shape']
  beta <- param['scale']

  x_nonzero <- x[x>0]
  loglik_nonzero <- log((1-pi0)*dgamma(x_nonzero, shape=alpha, scale=beta)) |> sum()
  loglik_zero <- sum(x==0) * log(pi0)

  return(loglik_nonzero + loglik_zero)
}


#' sum of log likelihoods for negative binomial distribution
#' @param param a named vector with two elements. One is the mean (mu) and the other is the size parameter (size)
#' @param x vector of observed counts
#' @importFrom stats dnbinom
lognb <- function(param, x){
  mu <- param["mu"]
  size <- param["size"]
  return(sum(dnbinom(x, mu=mu, size=size, log=T)))
}


#' sum of log likelihoods for zero inflated negative binomial distribution
#' @param param a named vector with three elements. The first one is the zero probability (pi0).
#' The other two are  the mean (mu) and the size parameter (size) of negative binomial distribution.
#' @param x vector of observed counts
#' @importFrom stats dnbinom
logzinb <- function(param, x){

  pi0 <- param["pi0"]
  mu <- param["mu"]
  size <- param["size"]

  x_nonzero <- x[x>0]
  loglik_nonzero <- log((1-pi0)*dnbinom(x_nonzero, mu=mu, size=size)) |> sum()
  loglik_zero <- sum(x==0) * log(pi0 + (1-pi0)*dbinom(0, mu=mu, size=size))

  return(loglik_nonzero + loglik_zero)

}


#' sum of log likelihoods for log normal distribution
#' @param param a named vector with three elements. The first one is the meanlog (mu) and the other is the sdlog (sigma)
#' @param x vector of observed counts
#' @importFrom stats dlnorm
loglnorm <- function(param, x){
  mu <- param["mu"]
  sigma <- param["sigma"]
  return(sum(dlnorm(x, meanlog=mu, sdlog=sigma, log=T)))
}


#' sum of log likelihoods for zero inflated log normal distribution
#' @param param a named vector with three elements. The first one is the meanlog (mu) and the other is the sdlog (sigma)
#' @param x vector of observed counts
#' @importFrom stats dlnorm
ziloglnorm <- function(param, x){
  pi0 <- param["pi0"]
  mu <- param["mu"]
  sigma <- param["sigma"]

  x_nonzero <- x[x>0]
  loglik_nonzero <- log((1-pi0)*dlnorm(x_nonzero, meanlog=mu, sdlog=sigma)) |> sum()
  loglik_zero <- sum(x==0) * log(pi0)

  return(loglik_nonzero + loglik_zero)

}






# transitioning between parameters and moments


#' convert parameters to mean and variances for a given distribution
#'
#' @param param a named vector of distribution parameters
#'@param dist the name of distribution, can be poisson(pois), gamma(gamma), negative binomial(nb), log normal (lnorm)
#' and their zero inflated versions(zipois, zigamma, zinb, zilnorm)
#'
#' @details
#' For poisson distribution, the name need to be "lambda". For negative binomial distribution,
#' the name need to be "mu" and "size". For gamma distribution,  the name needs to be "shape" and "scale".
#' For log normal distribution, the name needs to be "mu" and "sigma". Zero inflated distributions need to have "pi0".
#'
#'
#'@returns a named vector with mean and variances
param2mom <- function(param, dist=c("pois", "gamma", "nb", "lnorm",
                                    "zipois", "zigamma", "zinb", "zilnorm")){

  dist <- match.arg(dist)
  # helpter function for calculating mean and variance when zero inflation involved
  zidist <- function(pi0, nzmean, nzvar){
    output <- c(pi0=unname(pi0), mnz=unname(nzmean), vnz=unname(nzvar))
    meanval <- (1-pi0)*nzmean
    varval <- (1-pi0)*nzvar + pi0*(1-pi0)*(nzmean^2)
    output["m"] <- meanval
    output["v"] <- varval
    return(output)
  }

  if(dist == "pois" | dist == "zipois"){
    stopifnot("lambda" %in% names(param))
    nzmeanval <- param["lambda"] |> unname()
    nzvarval <- param["lambda"] |> unname()
    if(dist == "zipois"){
      stopifnot("pi0" %in% names(param))
      output <- zidist(pi0=param["pi0"], nzmean=nzmeanval, nzvar=nzvarval)
    } else{
      output <- c(m=nzmeanval, v=nzvarval)
    }
  } else if (dist=="gamma" | dist == "zigamma"){
    stopifnot("shape" %in% names(param) & "scale" %in% names(param))
    nzmeanval <- unname(param["shape"] * param["scale"])
    nzvarval <- unname(param["shape"] * (param["scale"]^2))
    if(dist == "zigamma"){
      stopifnot("pi0" %in% names(param))
      output <- zidist(pi0=param["pi0"], nzmean=nzmeanval, nzvar=nzvarval)
    } else{
      output <- c(m=nzmeanval, v=nzvarval)
    }
  } else if (dist == "nb" | dist == "zinb"){
    stopifnot("mu" %in% names(param) & "size" %in% names(param))
    nzmeanval <- unname(param["mu"])
    nzvarval <- unname(param["mu"] + (param["mu"]^2) / param["size"])
    if(dist == "zinb"){
      stopifnot("pi0" %in% names(param))
      output <- zidist(pi0=param["pi0"], nzmean=nzmeanval, nzvar=nzvarval)
    } else{
      output <- c(m=nzmeanval, v=nzvarval)
    }
  } else if (dist == "lnorm" | dist == "zilnorm"){
    stopifnot("mu" %in% names(param) & "sigma" %in% names(param))
    nzmeanval <- exp(param["mu"] + 0.5 * param["sigma"]^2) |> unname()
    nzvarval <- unname((exp(param["sigma"]^2) - 1) * exp(2*param["mu"] + param["sigma"]^2))
    if(dist == "zilnorm"){
      stopifnot("pi0" %in% names(param))
      output <- zidist(pi0=param["pi0"], nzmean=nzmeanval, nzvar=nzvarval)
    } else{
      output <- c(m=nzmeanval, v=nzvarval)
    }
  }

  return(output)
}


#' Convert mean and variances to parameters of a distribution
#'@param moments a named vector with mean, variance and zero inflation probability
#'@param dist distribution name, can be one of "pois", "gamma", "nb", "lnorm", "zipois", "zigamma", "zinb", "zilnorm"
#'@returns a named vector with the distribution parameters
#'@details
#' If the distribution does not have zero inflation components, the names in the vector are "m" and "v".
#' Otherwise the names need to have "pi0", "mnz" and "vnz". Poisson distribution has mean parameter "lambda".
#' Gamma distribution has two parameters "shape" and "scale". Negative binomial distribution has two parameters "mu" and "size".
#' lognormal distribution has two parameters "mu" and "sigma"
mom2param <- function(moments, dist=c("pois", "gamma", "nb", "lnorm",
                                    "zipois", "zigamma", "zinb", "zilnorm")){

  dist <- match.arg(dist)

  if (grepl("zi", dist)){
    stopifnot("pi0" %in% names(moments))
    key_mean <- "mnz"
    key_var <- "vnz"
  } else{
    key_mean <- "m"
    key_var <- "v"
  }


  if (dist == "pois" | dist =="zipois"){
    stopifnot(key_mean %in% names(moments) & key_var %in% names(moments))
    params <- c(lambda=unname(moments[key_mean]))
  } else if (dist == "gamma" | dist=="zigamma"){
    stopifnot(key_mean %in% names(moments) & key_var %in% names(moments))
    shape <- unname(moments[key_mean]^2 / moments[key_var])
    scale <- unname(moments[key_var] / moments[key_mean])
    params <- c(shape=shape, scale=scale)
  } else if (dist == "nb" | dist == "zinb"){
    stopifnot(key_mean %in% names(moments) & key_var %in% names(moments))
    mu <- unname(moments[key_mean])
    size <- unname(moments[key_mean]^2 / (moments[key_var] - moments[key_mean]))
    params <- c(mu=mu, size=size)
  } else if (dist == "lnorm" | dist == "zilnorm"){
    stopifnot(key_mean %in% names(moments) & key_var %in% names(moments))
    sigma2 <- unname(log(moments[key_var]/(moments[key_mean]^2)+1))
    mu <- unname(log(moments[key_mean]) - sigma2/2)
    params <- c(mu=mu, sigma=sqrt(sigma2))
  }

  if (grepl("zi", dist)){
    params["pi0"] <- moments["pi0"]
  }
  return(params)
}



#TODO: fit distribution function


#TODO: given quantiles, get distribution values


