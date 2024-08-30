library(tidyverse)
library(brms)
library(mrpaw)


# 1. Read in cleaned data ------------------------------------------------------

survey_df <- readRDS("survey_df.RDS")
poststrat_df <- readRDS("poststrat_df.RDS") 


# 2. Specify model -------------------------------------------------------------

cat_regressors <- c("gender", "age_group", "stateicp", "education_level")

regressor_string <- " ~ gender + (1 | age_group) + (1 | stateicp) + (1 | education_level)"

model_string <- paste0("vote_biden", regressor_string)

logit_prior <- c(prior(normal(0, 5), class = b))



# 3. Run MCMC ------------------------------------------------------------------

num_draws <- 3000 
force_rerun <- F
brm_logit_filename <- "brms_logit_fit.Rdata"

if (!file.exists(brm_logit_filename) || force_rerun) {
  stan_time <- Sys.time()
  logit_post <- brm(formula(model_string), 
                    family = bernoulli(link = "logit"),
                    data = survey_df,
                    prior = logit_prior,
                    control = list(adapt_delta = 0.99),
                    chains = 4, cores = 4, seed = 1543, 
                    warmup = 500, iter = num_draws)
  stan_time <- Sys.time() - stan_time
  print(stan_time)
  save(stan_time, logit_post, file = brm_logit_filename)
} else {
  load(brm_logit_filename)
}

print(logit_post)


# 4. MrP and implicit weights --------------------------------------------------

force_rerun <- F
mrp_aw_filename <- "mrp_aw.Rdata"

if (!file.exists(mrp_aw_filename) || force_rerun) {
  mrp_aw <- GetLogitMCMCWeights(logit_post, 
                                survey_df = survey_df, 
                                pop_df = poststrat_df, 
                                pop_w = poststrat_df$prop) 
  save(mrp_aw, file = mrp_aw_filename)
} else {
  load(mrp_aw_filename)
}

