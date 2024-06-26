% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weights_mcmc_lib.R
\name{GetOLSMCMCWeights}
\alias{GetOLSMCMCWeights}
\title{Get MrPaw weights for the logistic MCMC estimator.}
\usage{
GetOLSMCMCWeights(lin_post, survey_df, pop_df, pop_w = NULL)
}
\arguments{
\item{lin_post}{The output of \code{brm(..., survey_df, family=gaussian())}}

\item{survey_df}{The survey dataframe}

\item{pop_df}{The population dataframe}

\item{pop_w}{Optional.  The weight given to each row of pop_df.  Defaults to ones.}
}
\value{
Draws from the MrP estimate, and the weight vector
whose n-th entry is d E\link{MrP | X, Y} / d y_n.
}
\description{
Get MrPaw weights for the logistic MCMC estimator.
}
