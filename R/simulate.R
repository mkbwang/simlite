

#' simulate from multivariate normal distributions
#'
#' @param n sample size
#' @param mu mean vector, by default zero
#' @param Sigma covariance matrix whose diagonals are ones. Be default identity matrix
#' @importFrom MASS mvrnorm
gausscopula <- function(n=10, Sigma=diag(nrow=20)){

  stopifnot(all(diag(Sigma) == 1))
  gaussrv <- mvrnorm(n=n, mu=rep(0, nrow(Sigma)), Sigma=Sigma)
  cdf <- pnorm(gaussrv)
  return(cdf)

}



#' generate counts based on the given distribution and the quantiles
#' @param q vector of quantiles
#' @param param distribution parameters. A named vector
#' @param dist distribution name. Choose between "pois", "gamma", "nb", "lnorm", "zipois", "zigamma", "zinb", "zilnorm"
#' @importFrom stats qpois qgamma qnbinom qlnorm
simcount <- function(q, param, dist=c("pois", "gamma", "nb", "lnorm",
                                       "zipois", "zigamma", "zinb", "zilnorm")){

  dist <- match.arg(dist)
  if (grepl("zi", dist)){ # if there is zero inflation
    stopifnot("pi0" %in% names(param))
    q_rescaled <- (q-param["pi0"]) / (1-param["pi0"])
    zero_mask <- 1*(q_rescaled > 0)
    q_rescaled[q_rescaled < 0] <- 0
  } else{
    zero_mask <- 1
    q_rescaled <- q
  }

  if (dist == "pois" | dist == "zipois"){
    stopifnot("lambda" %in% names(param))
    counts <- qpois(q_rescaled, lambda=param["lambda"]) * zero_mask
  } else if (dist == "gamma" | dist == "zigamma"){
    stopifnot("shape" %in% names(param) & "scale" %in% names(param))
    counts <- qgamma(q_rescaled, shape=param["shape"], scale=param["scale"]) * zero_mask
  } else if (dist == "nb" | dist == "zinb"){
    stopifnot("mu" %in% names(param) & "size" %in% names(param))
    counts <- qnbinom(q_rescaled, mu=param["mu"], size=param["size"]) * zero_mask
  } else if (dist == "lnorm" | dist == "zilnorm"){
    stopifnot("mu" %in% names(param) & "sigma" %in% names(param))
    counts <- qlnorm(q_rescaled, meanlog=param["mu"], sdlog = param["sigma"]) * zero_mask
  }

  return(counts)
}



#' simulate count matrix given the quantiles and the distribution type and parameters.
#' @param copula a matrix of quantile values, nsample * nfeature
#' @param params distribution parameters for each feature
#' @param dist distributions
simcountmat <- function(copula, params, dist=c("pois", "gamma", "nb", "lnorm",
                                               "zipois", "zigamma", "zinb", "zilnorm")){

  simulated_counts <- mapply(simcount, asplit(copula, 2), asplit(params, 2),
                             MoreArgs=list(dist=dist))

  return(simulated_counts)

}



