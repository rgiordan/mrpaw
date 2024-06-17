library(testthat)
library(devtools)
library(mrpaw)

#setwd("/home/rgiordan/Documents/git_repos/SurveyWeighting/mrpaw/mrpaw")
#source("tests/testthat/helper.R")
setwd("~/Documents/git_repos/mrpaw/")

devtools::load_all()


post_filename <- test_path("mcmc_cache/logit_post_test.rds")
file.exists(post_filename)

getwd()
#testthat::test_file("tests/testthat/test_mcmc.R")
#testthat::test_local("tests/testthat/test_mcmc.R")
testthat::test_file("tests/testthat/test_mcmc.R")



method <- "ols_response_name"

if (method %in% c("ols", "ols_response_name")) {
  MrPawFunction <- GetOLSMCMCWeights
} else if (method %in% c("logit", "logit_response_name")) {
  MrPawFunction <- GetLogitMCMCWeights
} else {
  expect_true(FALSE, sprintf("Unknown method %s", method))
}

rds_load <- SafeLoadPosterior(method)
post <- rds_load$post
sim_data <- rds_load$sim_data
agg_list <- rds_load$agg_list

y_col <- f_lhs(as.formula(formula(post)))
agg_list <- AggregateSimulationData(sim_data, y_col)

# Test that this runs and produces weights of the correct length.
mcmc_mrp <- MrPawFunction(
  post, 
  sim_data$survey_df, 
  agg_list$pop_agg_df, 
  pop_w=agg_list$pop_agg_df$w)

expect_true(length(mcmc_mrp$w) == nrow(sim_data$survey_df))