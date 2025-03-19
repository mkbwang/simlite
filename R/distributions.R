

expit <- function(y){
  exp(y)/(1+exp(y))
}

logit <- function(y){
  log(y / (1-y))
}

# log likelihoods

#'sum of log likelihoods for a poisson distribution
#'@param param a named vector with one element "loglambda"
#'@param x vector of observed counts
#'@return sum of log likelihoods
#'@importFrom stats dpois
logpois <- function(param, x){
  lambda <- exp(param["loglambda"])
  return(sum(dpois(x, lambda, log=T)))
}


#'sum of log likelihoods for a zero inflated poisson distribution
#'@param param a named vector two elements. One is the zero probability (logitpi0) and the other is the mean of nonzero part (loglambda)
#'@param x vector of observed counts
#'@importFrom stats dpois
logzipois <- function(param, x){
  pi0 <- expit(param["logitpi0"])
  lambda <- exp(param["loglambda"])

  x_nonzero <- x[x>0]
  loglik_nonzero <- sum(log(1-pi0)+dpois(x_nonzero, lambda, log=T))
  loglik_zero <- sum(x==0) * log(pi0 + (1-pi0)*dpois(0, lambda))

  return(loglik_nonzero + loglik_zero)

}

#'sum of log likelihoods for a gamma distribution
#'@param param a named vector two elements. One is the shape parameter (logalpha) and the other is the scale (logbeta)
#'@param x vector of observed counts
#'@importFrom stats dgamma
loggamma <- function(param, x){
  alpha <- exp(param["logshape"])
  beta <- exp(param["logscale"])
  return(sum(dgamma(x, shape=alpha, scale=beta, log=T)))
}


#'sum of log likelihoods for a zero inflated gamma distribution
#'@param param a named vector with three elements. One is the zero probability logitpi0.
#'The other two are shape parameter (logalpha) and the other is the scale (logbeta)
#'@param x vector of observed counts
#'@importFrom stats dgamma
#'
logzigamma <- function(param, x){

  pi0 <- expit(param['logitpi0'])
  alpha <- exp(param['logshape'])
  beta <- exp(param['logscale'])

  x_nonzero <- x[x>0]
  loglik_nonzero <- sum(log(1-pi0)+dgamma(x_nonzero, shape=alpha, scale=beta, log=T))
  loglik_zero <- sum(x==0) * log(pi0)

  return(loglik_nonzero + loglik_zero)
}


#' sum of log likelihoods for negative binomial distribution
#' @param param a named vector with two elements. One is the mean (logmu) and the other is the size parameter (logsize)
#' @param x vector of observed counts
#' @importFrom stats dnbinom
lognb <- function(param, x){
  mu <- exp(param["logmu"])
  size <- exp(param["logsize"])
  return(sum(dnbinom(x, mu=mu, size=size, log=T)))
}


#' sum of log likelihoods for zero inflated negative binomial distribution
#' @param param a named vector with three elements. The first one is the zero probability (logitpi0).
#' The other two are  the mean (logmu) and the size parameter (logsize) of negative binomial distribution.
#' @param x vector of observed counts
#' @importFrom stats dnbinom
logzinb <- function(param, x){

  pi0 <- expit(param["logitpi0"])
  mu <- exp(param["logmu"])
  size <- exp(param["logsize"])

  x_nonzero <- x[x>0]
  loglik_nonzero <- sum(log(1-pi0)+dnbinom(x_nonzero, mu=mu, size=size, log=T))
  loglik_zero <- sum(x==0) * log(pi0 + (1-pi0)*dnbinom(0, mu=mu, size=size))

  return(loglik_nonzero + loglik_zero)

}


#' sum of log likelihoods for log normal distribution
#' @param param a named vector with three elements. The first one is the meanlog (logmu) and the other is the sdlog (logsigma)
#' @param x vector of observed counts
#' @importFrom stats dlnorm
loglnorm <- function(param, x){
  mu <- param["mu"]
  sigma <- exp(param["logsigma"])
  return(sum(dlnorm(x, meanlog=mu, sdlog=sigma, log=T)))
}


#' sum of log likelihoods for zero inflated log normal distribution
#' @param param a named vector with three elements. The first one is the meanlog (mu) and the other is the sdlog (sigma)
#' @param x vector of observed counts
#' @importFrom stats dlnorm
logzilnorm <- function(param, x){
  pi0 <- expit(param["logitpi0"])
  mu <- param["mu"]
  sigma <- exp(param["logsigma"])

  x_nonzero <- x[x>0]
  loglik_nonzero <- sum(log(1-pi0)+dlnorm(x_nonzero, meanlog=mu, sdlog=sigma, log=T))
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
#'@export
param2mom <- function(param, dist=c("pois", "gamma", "nb", "lnorm",
                                    "zipois", "zigamma", "zinb", "zilnorm")){

  dist <- match.arg(dist)
  # helper function for calculating mean and variance when zero inflation involved
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
#'@export
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
    # just in case variance is smaller than mean
    if (size < 0) size <- 100
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



#' Given a vector of counts, fit a distribution
#'@param x vector of observed values
#'@param dist distribution. Can be chosen from pois, gamma, nb, lnorm, zipois, zigamma, zinb and zilnorm
#'@param cutoff values larger than certain quantile are excluded from fitting the model, default 97%
#'@importFrom stats optim var
#'@export
fitdistr <- function(x, dist=c("pois", "gamma", "nb", "lnorm",
                               "zipois", "zigamma", "zinb", "zilnorm"),
                     cutoff=0.97){

  # remove extremely large values
  cutoff_value <- quantile(x, cutoff)
  x <- x[x <= cutoff_value]

  dist <- match.arg(dist)
  control = list(fnscale = -1)

  # set up initial parameters
  if(grepl("zi", dist)){
    observed_pi0 <- mean(x == 0)*0.5
    observed_mnz <- sum(x) / (length(x) - 0.5*sum(x==0))
    observed_vnz <- var(x)
    observed_mom <- c(pi0=observed_pi0, mnz=observed_mnz, vnz=observed_vnz)
  } else{
    observed_m <- mean(x)
    observed_v <- var(x)
    observed_mom <- c(m=observed_m, v=observed_v)
  }

  if (dist == "pois"){
    observed_param <- mom2param(observed_mom, dist="pois")
    t_observed_param <- log(observed_param["lambda"])
    names(t_observed_param) <- "loglambda"
    result <- optim(par=t_observed_param, fn=logpois, x=x, control=control, method="BFGS")
    fitted_param <- result$par
    fitted_param <- c(exp(fitted_param["loglambda"]))
    names(fitted_param) <- "lambda"
  } else if (dist == "zipois") {
    observed_param <- mom2param(observed_mom, dist="zipois")
    logit_pi0 <- max(logit(observed_param["pi0"]), -20)
    log_lambda <- log(observed_param["lambda"])
    t_observed_param <- c(logit_pi0, log_lambda)
    names(t_observed_param) <- c("logitpi0", "loglambda")
    result <- optim(par=t_observed_param, fn=logzipois, x=x, control=control, method="BFGS")
    fitted_param <- result$par
    fitted_param <- c(expit(fitted_param["logitpi0"]), exp(fitted_param["loglambda"]))
    names(fitted_param) <- c("pi0", "lambda")
  } else if (dist == "gamma"){
    stopifnot("Data points contain zeros, but Gamma distribution does not allow zeros." = (sum(x <= 0) == 0))
    observed_param <- mom2param(observed_mom, dist="gamma")
    log_shape <- log(observed_param["shape"])
    log_scale <- log(observed_param["scale"])
    t_observed_param <- c(log_shape, log_scale)
    names(t_observed_param) <- c("logshape", "logscale")
    result <- optim(par=t_observed_param, fn=loggamma, x=x, control=control, method="BFGS")
    fitted_param <- result$par
    fitted_param <- c(exp(fitted_param["logshape"]), exp(fitted_param["logscale"]))
    names(fitted_param) <- c("shape", "scale")
  } else if (dist == "zigamma"){
    observed_param <- mom2param(observed_mom, dist="zigamma")
    logit_pi0 <- max(logit(observed_param["pi0"]), -20)
    log_shape <- log(observed_param["shape"])
    log_scale <- log(observed_param["scale"])
    t_observed_param <- c(logit_pi0, log_shape, log_scale)
    names(t_observed_param) <- c("logitpi0", "logshape", "logscale")
    result <- optim(par=t_observed_param, fn=logzigamma, x=x, control=control, method="BFGS")
    fitted_param <- result$par
    fitted_param <- c(expit(fitted_param["logitpi0"]), exp(fitted_param["logshape"]), exp(fitted_param["logscale"]))
    names(fitted_param) <- c("pi0", "shape", "scale")
  } else if (dist == "nb"){
    observed_param <- mom2param(observed_mom, dist="nb")
    log_mu <- log(observed_param["mu"])
    log_size <- log(observed_param["size"])
    t_observed_param <- c(mu, log_size)
    names(t_observed_param) <- c("logmu", "logsize")
    result <- optim(par=t_observed_param, fn=lognb, x=x, control=control, method="BFGS")
    fitted_param <- result$par
    fitted_param <- c(exp(fitted_param["logmu"]), exp(fitted_param["logsize"]))
    names(fitted_param) <- c("mu", "size")
  } else if (dist == "zinb"){
    observed_param <- mom2param(observed_mom, dist="zinb")
    logit_pi0 <- max(logit(observed_param["pi0"]), -20)
    log_mu <- log(observed_param["mu"])
    log_size <- log(observed_param["size"])
    t_observed_param <- c(logit_pi0, log_mu, log_size)
    names(t_observed_param) <- c("logitpi0", "logmu", "logsize")
    result <- optim(par=t_observed_param, fn=logzinb, x=x, control=control, method="BFGS")
    fitted_param <- result$par
    fitted_param <- c(expit(fitted_param["logitpi0"]), exp(fitted_param["logmu"]), exp(fitted_param["logsize"]))
    names(fitted_param) <- c("pi0", "mu", "size")
  } else if (dist == "lnorm"){
    stopifnot("Data points contain zeros, but lognormal distribution does not allow zeros." = (sum(x <= 0) == 0))
    observed_param <- mom2param(observed_mom, dist="lnorm")
    # log_mu <- log(observed_param["mu"])
    log_sigma <- log(observed_param["sigma"])
    t_observed_param <- c(observed_param["mu"], log_sigma)
    names(t_observed_param) <- c("mu", "logsigma")
    result <- optim(par=t_observed_param, fn=loglnorm, x=x, control=control, method="BFGS")
    fitted_param <- result$par
    fitted_param <- c(fitted_param["mu"], exp(fitted_param["logsigma"]))
    names(fitted_param) <- c("mu", "sigma")
  } else if (dist == "zilnorm"){
    observed_param <- mom2param(observed_mom, dist="zilnorm")
    logit_pi0 <- max(logit(observed_param["pi0"]), -20)
    # log_mu <- log(observed_param["mu"])
    log_sigma <- log(observed_param["sigma"])
    t_observed_param <- c(logit_pi0, observed_param["mu"], log_sigma)
    names(t_observed_param) <- c("logitpi0", "mu", "logsigma")
    result <- optim(par=t_observed_param, fn=logzilnorm, x=x, control=control, method="BFGS")
    fitted_param <- result$par
    fitted_param <- c(expit(fitted_param["logitpi0"]), fitted_param["mu"], exp(fitted_param["logsigma"]))
    names(fitted_param) <- c("pi0", "mu", "sigma")
  }

  return(fitted_param)
}



#' fit distributions to multiple features in a count matrix at the same time
#'
#' @param count_mat nfeatures*nsamples
#' @param dist distribution name among "pois", "gamma", "nb", "lnorm", "zipois", "zigamma", "zinb", "zilnorm"
#' @param cutoff values larger than certain quantile are excluded from fitting the model, default 97%
#' @importFrom parallel detectCores makeCluster stopCluster clusterExport
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @export
fitdistr_mat <- function(count_mat, dist=c("pois", "gamma", "nb", "lnorm",
                                           "zipois", "zigamma", "zinb", "zilnorm"),
                         cutoff=0.97){

  dist <- match.arg(dist)
  nfeatures <- nrow(count_mat)
  numCores <- detectCores() - 1
  cl <- makeCluster(numCores)
  registerDoParallel(cl)
  i <- 1

  clusterExport(cl, varlist=c("fitdistr", "mom2param", "logpois",
                              "logzipois", "loggamma", "logzigamma", "lognb",
                              "logzinb", "loglnorm", "logzilnorm", "expit", "logit"))

  params <- foreach(i=1:nfeatures, .combine=rbind) %dopar%{
    fitted_param <- fitdistr(count_mat[i, ], dist=dist, cutoff=cutoff)
    fitted_param
  }

  stopCluster(cl)
  return(params)

}


