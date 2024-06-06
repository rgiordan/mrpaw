library(testthat)


GetTestData <- function() {
  set.seed(42)

  n_groups <- 2
  degree <- 2
  n_obs <- 1000
  n_obs_pop <- 10 * n_obs

  sim_data <- SimulateSurveyData(n_groups, n_obs, n_obs_pop, degree=degree)
  return(sim_data)
}


AggregateSimulationData <- function(sim_data) {

  ##########################################################
  # Aggregate across groups to check 

  group_effects <- sim_data$group_effects
  g_cols <- group_effects$g_cols

  survey_agg_df <-
    sim_data$survey_df %>%
    group_by(pick(all_of(c(g_cols, "s")))) %>%
    summarize(count=n(), ybar=mean(y), ey=mean(ey), .groups="drop") %>%
    mutate(w=count / sum(count))

  pop_agg_df <-
    sim_data$pop_df %>%
    group_by(pick(all_of(c(g_cols, "s")))) %>%
    summarize(count=n(), ybar=mean(y), ey=mean(ey), .groups="drop") %>%
    mutate(w=count / sum(count))

  # Compare the averages and weights directly from the population
  # to the sample.
  joint_df <- left_join(
    pop_agg_df,
    survey_agg_df,
    by=c(g_cols, "s"),
    suffix=c("_pop", "_sur")) %>%
    mutate(w_opt=w_pop / count_sur)

  return(list(joint_df=joint_df, survey_agg_df=survey_agg_df, pop_agg_df=pop_agg_df))
}


AssertNearlyEqual <- function(x, y, tol=1e-9, desc=NULL) {
  diff_norm <- max(abs(x - y))
  if (is.null(desc)) {
    info_str <- sprintf("%e > %e", diff_norm, tol)
  } else {
    info_str <- sprintf("%s: %e > %e", desc, diff_norm, tol)
  }
  expect_true(diff_norm < tol, info=info_str)
}


AssertNearlyZero <- function(x, tol=1e-15, desc=NULL) {
  x_norm <- max(abs(x))
  if (is.null(desc)) {
    info_str <- sprintf("%e > %e", x_norm, tol)
  } else {
    info_str <- sprintf("%s: %e > %e", desc, x_norm, tol)
  }
  expect_true(x_norm < tol, info=info_str)
}






################################################################
# Load saved posterior samples (or run if needed)

RunAndCachePosterior <- function(
  sim_data, reg_form, post_filename, family, num_draws=100) {
    stan_time <- Sys.time()
    post <- brm(formula(reg_form), sim_data$survey_df, family=family,
                chains=4, cores=4, seed=1543, 
                iter=num_draws)
    stan_time <- Sys.time() - stan_time
    print(stan_time)
    saveRDS(list(
      post=post, 
      sim_data=sim_data),
      file=post_filename)
}


SafeLoadPosterior <- function(method) {
  if (method == "ols") {
    family <- gaussian()
    post_filename <- test_path("mcmc_cache/ols_post_test.rds")
  } else if (method == "logit") {
    family <- bernoulli(link="logit")
    post_filename <- test_path("mcmc_cache/logit_post_test.rds")
  } else {
    expect_true(FALSE, sprintf("Unknown method %s", method))
  }

  print(post_filename)

  if (!file.exists(post_filename)) {
    print(sprintf("Could not locate file %s", post_filename))
    print("Re-running MCMC and saving result to rds file")
    sim_data <- GetTestData()
    group_effects <- sim_data$group_effects
    g_cols <- sim_data$group_effects$g_cols
    g_sum <- paste(group_effects$g_cols, collapse=" + ")
    reg_form <-sprintf( "y ~ 1 + (%s)^%d", g_sum, degree=2)
    RunAndCachePosterior(sim_data, reg_form, post_filename, family=family, num_draws=100)
  } else {
        print(sprintf("Found cached posterior samples in %s", post_filename))
  }

  print(sprintf("Reading posterior from file %s", post_filename))
  rds_load <- readRDS(post_filename)
  return(rds_load)
}

