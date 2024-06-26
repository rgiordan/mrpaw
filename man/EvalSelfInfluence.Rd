% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diagnostics_lib.R
\name{EvalSelfInfluence}
\alias{EvalSelfInfluence}
\title{Evaluate the self-influence matrix dyhat / dy}
\usage{
EvalSelfInfluence(post, survey_df)
}
\arguments{
\item{post}{The output of \code{brm(..., survey_df, family=binomial(link="logit"))}}

\item{survey_df}{The survey dataframe}
}
\description{
Evaluate the self-influence matrix dyhat / dy
}
