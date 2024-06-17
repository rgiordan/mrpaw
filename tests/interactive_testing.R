library(testthat)
library(devtools)
library(mrpaw)

setwd("/home/rgiordan/Documents/git_repos/SurveyWeighting/mrpaw/mrpaw")
#source("tests/testthat/helper.R")

devtools::load_all()


post_filename <- test_path("mcmc_cache/logit_post_test.rds")
file.exists(post_filename)

getwd()
#testthat::test_file("tests/testthat/test_mcmc.R")
#testthat::test_local("tests/testthat/test_mcmc.R")
testthat::test_file("tests/testthat/test_mcmc.R")



method <- "ols"


print(sprintf("Testing MCMC for method %s", method))
if (method == "ols") {
  MrPawFunction <- GetOLSMCMCWeights
} else if (method == "logit") {
  MrPawFunction <- GetLogitMCMCWeights
} else {
  expect_true(FALSE, "This should never happen")
}

rds_load <- SafeLoadPosterior(method)
post <- rds_load$post
sim_data <- rds_load$sim_data
agg_list <- AggregateSimulationData(sim_data)

# Sanity check that I'm computing the log likelihood correctly
# (I'm not taking into account the prior so there will be some small mismatch)
ols_ll_draws <- GetOLSLikelihoodComponentDraws(post, sim_data$survey_df)

sigma_draws <- ols_ll_draws$sigma_draws
resid_draws <- ols_ll_draws$resid_draws
lp_draws_check <- post %>% spread_draws(lp__) %>% pull(lp__)
lp_mat <- -0.5 * (resid_draws^2) / (sigma_draws^2) - log(sigma_draws)
lp_draws <- apply(lp_mat, FUN=sum, MARGIN=1)
expect_true(cor(lp_draws, lp_draws_check) > 0.99)
