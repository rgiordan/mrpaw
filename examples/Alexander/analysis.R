library(tidyverse)
library(brms)
library(mrpaw)


# 1. Read in cleaned data ------------------------------------------------------

survey_df <- readRDS("survey_df.RDS")
poststrat_df <- readRDS("poststrat_df.RDS") 


# 2. Specify model -------------------------------------------------------------

cat_regressors <- c("age_group", "educ_group", "state_name", "decade_married")

regressor_string <- paste0(" ~ (1 | age_group) + (1 | educ_group) + ",
                           "(1 | state_name) + (1 | decade_married)")

model_string <- paste0("kept_name", regressor_string)

# logit_prior <-


# 3. Run MCMC ------------------------------------------------------------------

num_draws <- 2000 
force_rerun <- F
brm_logit_filename <- "brms_logit_fit.Rdata"

if (!file.exists(brm_logit_filename) || force_rerun) {
  stan_time <- Sys.time()
  logit_post <- brm(formula(model_string), 
                    family = bernoulli(link = "logit"),
                    data = survey_df,
                    #prior = logit_prior,
                    control = list(adapt_delta = 0.99, 
                                   max_treedepth = 11),
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

force_rerun <- T
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




# EXTRA (negative weights and non-monotonicity) --------------------------------

# Since there are many negative weights, we check that the MrP estimate is not 
# monotonically decreasing in the number of responses = 0.

# Flip to 0 the response of individuals that chose response 1 and that have 
# negative weights that are large in magnitude
survey_df_flip = survey_df
survey_df_flip[which(mrp_aw$w < -0.0003 & survey_df$kept_name == 1), 
               'kept_name'] = 0


# Run MCMC with perturbed dataset
num_draws <- 2000 
force_rerun <- F
brm_logit_filename <- "brms_logit_fit_flip.Rdata"

if (!file.exists(brm_logit_filename) || force_rerun) {
  stan_time <- Sys.time()
  logit_post_flip <- brm(formula(model_string), 
                    family = bernoulli(link = "logit"),
                    data = survey_df_flip,
                    #prior = logit_prior,
                    control = list(adapt_delta = 0.99, 
                                   max_treedepth = 11),
                    chains = 4, cores = 4, seed = 1543, 
                    warmup = 500, iter = num_draws)
  stan_time <- Sys.time() - stan_time
  print(stan_time)
  save(stan_time, logit_post_flip, file = brm_logit_filename)
} else {
  load(brm_logit_filename)
}

# print(logit_post_flip)


# MrP and implicit weights for perturbed dataset
force_rerun <- F
mrp_aw_filename <- "mrp_aw_flip.Rdata"

if (!file.exists(mrp_aw_filename) || force_rerun) {
  mrp_aw_flip <- GetLogitMCMCWeights(logit_post_flip, 
                                survey_df = survey_df, 
                                pop_df = poststrat_df, 
                                pop_w = poststrat_df$prop)
  save(mrp_aw_flip, file = mrp_aw_filename)
} else {
  load(mrp_aw_filename)
}

# Compare old and new MrP
cat(paste0('------------------------------------------------------------------------',
'\n', 'MrP estimate on original data: ',  round(mean(mrp_aw$mrp_draws), 3), 
'\n', 'MrP estimate on perturbed data (increased number of 0 responses): ',
round(mean(mrp_aw_flip$mrp_draws), 3), '\n',
'------------------------------------------------------------------------'))


