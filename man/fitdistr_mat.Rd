% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distributions.R
\name{fitdistr_mat}
\alias{fitdistr_mat}
\title{fit distributions to multiple features in a count matrix at the same time}
\usage{
fitdistr_mat(
  count_mat,
  dist = c("pois", "gamma", "nb", "lnorm", "zipois", "zigamma", "zinb", "zilnorm"),
  cutoff = 0.97,
  ncores = 1
)
}
\arguments{
\item{count_mat}{nfeatures*nsamples}

\item{dist}{distribution name among "pois", "gamma", "nb", "lnorm", "zipois", "zigamma", "zinb", "zilnorm"}

\item{cutoff}{values larger than certain quantile are excluded from fitting the model, default 97\%}

\item{ncores}{number of cores for parallel processing, by default 1}
}
\description{
fit distributions to multiple features in a count matrix at the same time
}
\examples{
data("gut_experiment_subset")
count_data <- SummarizedExperiment::assay(gut_experiment_subset)
subset_genes <- sample(nrow(count_data), 30)
subset_count_data <- round(count_data[subset_genes, ])
pois_params <- fitdistr_mat(subset_count_data, dist="pois", ncores=4)
}
