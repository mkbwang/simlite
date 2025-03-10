test_that("Poisson", {

  example_pois <- c(lambda=4)
  pois_mom <- param2mom(example_pois, dist="pois")
  pois_param <- mom2param(pois_mom, dist="pois")
  expect_equal(pois_param["lambda"], example_pois["lambda"], tolerance=1e-4)

})


test_that("Zero Inflated Poisson", {

  example_zipois <- c(pi0=0.3, lambda=3)
  zipois_mom <- param2mom(example_zipois, dist="zipois")
  zipois_param <- mom2param(zipois_mom, dist="zipois")
  expect_equal(zipois_param["pi0"], example_zipois["pi0"], tolerance=1e-4)
  expect_equal(zipois_param["lambda"], example_zipois["lambda"], tolerance=1e-4)

})


test_that("Gamma", {

  example_gamma <- c(shape=3, scale=2)
  gamma_mom <- param2mom(example_gamma, dist="gamma")
  gamma_param <- mom2param(gamma_mom, dist="gamma")
  expect_equal(gamma_param["shape"], example_gamma["shape"], tolerance=1e-4)
  expect_equal(gamma_param["scale"], example_gamma["scale"], tolerance=1e-4)

})


test_that("Zero Inflated Gamma", {

  example_gamma <- c(pi0=0.2, shape=3, scale=2)
  gamma_mom <- param2mom(example_gamma, dist="zigamma")
  gamma_param <- mom2param(gamma_mom, dist="zigamma")
  expect_equal(gamma_param["pi0"], example_gamma["pi0"], tolerance=1e-4)
  expect_equal(gamma_param["shape"], example_gamma["shape"], tolerance=1e-4)
  expect_equal(gamma_param["scale"], example_gamma["scale"], tolerance=1e-4)

})


test_that("Negative Binomial", {

  example_nb <- c(mu=10, size=1)
  nb_mom <- param2mom(example_nb, dist="nb")
  nb_param <- mom2param(nb_mom, dist="nb")
  expect_equal(nb_param["mu"], example_nb["mu"], tolerance=1e-4)
  expect_equal(nb_param["size"], example_nb["size"], tolerance=1e-4)

})



test_that("zero inflated negative binomial", {

  example_nb <- c(pi0=0.2, mu=10, size=1)
  nb_mom <- param2mom(example_nb, dist="zinb")
  nb_param <- mom2param(nb_mom, dist="zinb")
  expect_equal(nb_param["pi0"], example_nb["pi0"], tolerance=1e-4)
  expect_equal(nb_param["mu"], example_nb["mu"], tolerance=1e-4)
  expect_equal(nb_param["size"], example_nb["size"], tolerance=1e-4)

})


test_that("Log Normal", {

  example_lnorm <- c(mu=2, sigma=0.5)
  lnorm_mom <- param2mom(example_lnorm, dist="lnorm")
  lnorm_param <- mom2param(lnorm_mom, dist="lnorm")
  expect_equal(lnorm_param["mu"], example_lnorm["mu"], tolerance=1e-4)
  expect_equal(lnorm_param["sigma"], example_lnorm["sigma"], tolerance=1e-4)

})


test_that("Zero inflated Log Normal", {

  example_lnorm <- c(pi0=0.2, mu=2, sigma=0.5)
  lnorm_mom <- param2mom(example_lnorm, dist="zilnorm")
  lnorm_param <- mom2param(lnorm_mom, dist="zilnorm")
  expect_equal(lnorm_param["pi0"], example_lnorm["pi0"], tolerance=1e-4)
  expect_equal(lnorm_param["mu"], example_lnorm["mu"], tolerance=1e-4)
  expect_equal(lnorm_param["sigma"], example_lnorm["sigma"], tolerance=1e-4)

})


test_that("Fit Poisson", {
  set.seed(2025)
  nsample <- 400
  lambda <- 89
  counts <- rpois(n=nsample, lambda=lambda)
  fitted_params <- fitdistr(x=counts, dist="pois")
  expect_equal(unname(fitted_params["lambda"]), lambda, tolerance=0.1)
})


test_that("Fit Zero Inflated Poisson", {
  set.seed(2025)
  nsample <- 400
  pi0 <- 0.2
  lambda <- 2
  counts <- rbinom(n=nsample, size=1, prob=1-pi0)
  counts <- counts * rpois(n=nsample, lambda=lambda)
  fitted_params <- fitdistr(x=counts, dist="zipois")
  expect_equal(unname(fitted_params["pi0"]), pi0, tolerance=0.2)
  expect_equal(unname(fitted_params["lambda"]), lambda, tolerance=0.2)

})

test_that("Fit Gamma", {
  set.seed(2025)
  nsample <- 400
  shape <- 3
  scale <- 0.6
  counts <- rgamma(n=nsample, shape=shape, scale=scale)
  fitted_params <- fitdistr(x=counts, dist="gamma")
  expect_equal(unname(fitted_params["shape"]), shape, tolerance=0.1)
  expect_equal(unname(fitted_params["scale"]), scale, tolerance=0.2)
})

test_that("Fit Zero Inflated Gamma", {
  set.seed(2025)
  nsample <- 400
  pi0 <- 0.2
  shape <- 3
  scale <- 0.6
  counts <- rbinom(n=nsample, size=1, prob=1-pi0)
  counts <- counts * rgamma(n=nsample, shape=shape, scale=scale)
  fitted_params <- fitdistr(x=counts, dist="zigamma")
  expect_equal(unname(fitted_params["pi0"]), pi0, tolerance=0.2)
  expect_equal(unname(fitted_params["shape"]), shape, tolerance=0.2)
  expect_equal(unname(fitted_params["scale"]), scale, tolerance=0.2)
})


test_that("Fit Negative Binomial", {

  set.seed(2025)
  nsample <- 400
  mu <- 20
  size <- 1.5
  counts <- rnbinom(n=nsample, size=size, mu=mu)
  fitted_params <- fitdistr(x=counts, dist="nb")
  expect_equal(unname(fitted_params["mu"]), mu, tolerance=0.1)
  expect_equal(unname(fitted_params["size"]), size, tolerance=0.2)

})


test_that("Fit Zero Inflated Negative Binomial", {

  set.seed(2024)
  nsample <- 400
  pi0 <- 0.2
  mu <- 4
  size <- 1
  counts <- rbinom(n=nsample, size=1, prob=1-pi0)
  counts <- counts * rnbinom(n=nsample, size=size, mu=mu)
  fitted_params <- fitdistr(x=counts, dist="zinb")
  expect_equal(unname(fitted_params["pi0"]), pi0, tolerance=0.2)
  expect_equal(unname(fitted_params["mu"]), mu, tolerance=0.2)
  expect_equal(unname(fitted_params["size"]), size, tolerance=0.3)

})


test_that("Fit log normal", {

  set.seed(2024)
  nsample <- 400
  mu <- 4
  sigma <- 0.2
  counts <- rlnorm(n=nsample, meanlog=mu, sdlog=sigma)
  fitted_params <- fitdistr(x=counts, dist="lnorm")
  expect_equal(unname(fitted_params["mu"]), mu, tolerance=0.2)
  expect_equal(unname(fitted_params["sigma"]), sigma, tolerance=0.1)

})



test_that("Fit zero inflated log normal", {

  set.seed(2025)
  nsample <- 400
  pi0 <- 0.2
  mu <- 1.5
  sigma <- 0.2
  counts <- rbinom(n=nsample, size=1, prob=1-pi0)
  counts <- counts * rlnorm(n=nsample, meanlog=mu, sdlog=sigma)
  fitted_params <- fitdistr(x=counts, dist="zilnorm")
  expect_equal(unname(fitted_params["pi0"]), pi0, tolerance=0.2)
  expect_equal(unname(fitted_params["mu"]), mu, tolerance=0.1)
  expect_equal(unname(fitted_params["sigma"]), sigma, tolerance=0.2)

})
