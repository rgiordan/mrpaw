###################################################
# Various ways to use the MrP weights

library(einsum)
library(tidyverse)
library(brms)
library(tidybayes)



SafeGetYhatDraws <- function(mrpaw_list, post, survey_df) {
    if ("yhat_draws" %in% names(mrpaw_list)) {
        yhat_draws <- mrpaw_list$yhat_draws
    } else {
        yhat_draws <- posterior_epred(post, newdata=survey_df)
    }
    return(yhat_draws)
}


SafeGetEtaDraws <- function(mrpaw_list, post, survey_df) {
    if ("eta_draws" %in% names(mrpaw_list)) {
        eta_draws <- mrpaw_list$eta_draws
    } else {
        eta_draws <- posterior_linpred(post, newdata=survey_df)
    }
    return(eta_draws)
}


#' Evaluate importance sampling MCMC weights for new y values.
#' @param logit_post The output of `brm(..., survey_df, family=binomial(link="logit"))`
#' @param survey_df The survey dataframe
#' @param y_new A vector of new responses, the same length as the original response
#'
#' @export
GetLogitImportanceWeights <- function(logit_post, survey_df, y_new) {
  stopifnot(class(logit_post) == "brmsfit")
  CheckLogitFamily(logit_post)
  linpred <- posterior_linpred(logit_post, newdata=survey_df)
  y <- GetResponse(survey_df)
  stopifnot(ncol(linpred) == length(y))

  # Transpose to take advantage of column broadcasting
  ll_diff <- (t(linpred) * (y_new - y)) %>% t() %>%
    apply(X=., MARGIN=1, FUN=sum)
  
  w_is <- exp(ll_diff)
  ess_is <- sum(w_is)^2 / sum(w_is^2)

  return(list(ll_diff=ll_diff, w_is=w_is, ess=ess_is))
}



#' Evaluate importance sampling MCMC weights for new y values.
#' @param lin_post The output of `brm(..., survey_df, family=gaussian())`
#' @param survey_df The survey dataframe
#' @param y_new A vector of new responses, the same length as the original response
#'
#' @export
GetOLSImportanceWeights <- function(lin_post, survey_df, y_new) {
  stopifnot(class(lin_post) == "brmsfit")
  CheckOLSFamily(lin_post)
  
  ols_draws <- GetOLSLikelihoodComponentDraws(lin_post=lin_post, survey_df=survey_df)
  mu_draws <- ols_draws$yhat_draws
  sigma_draws <- ols_draws$sigma_draws

  # Transpose to take advantage of column broadcasting
  y1 <- y_new
  y0 <- GetResponse(lin_post)
  y1_ressq_draws <- (y1 - t(mu_draws))^2 %>% t() %>% apply(X=., MARGIN=1, FUN=sum)
  y0_ressq_draws <- (y0 - t(mu_draws))^2 %>% t() %>% apply(X=., MARGIN=1, FUN=sum)

  ll_diff <- -0.5 * (y1_ressq_draws - y0_ressq_draws) / (sigma_draws^2)
  
  w_is <- exp(ll_diff)
  ess_is <- sum(w_is)^2 / sum(w_is^2)
  
  return(list(ll_diff=ll_diff, w_is=w_is, ess=ess_is))
}



#' Evaluate the constant in the local expansion.
#' @param mrpaw_list The output of one of the Get*Weights functions
#' @param y The original vector of responses
#'
#' @export
EvalOffset <- function(mrpaw_list, y) {
    mrp_base <- mean(mrpaw_list$mrp_draws)
    offset <- sum(mrpaw_list$w * y) - mrp_base
    return(offset)
}



#' Evaluate the self-influence matrix dyhat / dy
#' @param post The output of `brm(..., survey_df, family=binomial(link="logit"))`
#' @param survey_df The survey dataframe
#'
#' @export
EvalSelfInfluence <- function(post, survey_df) {
    yhat_draws <- posterior_epred(post, newdata=survey_df)
    eta_draws <- posterior_linpred(post, newdata=survey_df)
    dyhat_dy <- cov(yhat_draws, eta_draws)
    return(dyhat_dy)
}



#' Evaluate the linear approximation to the parametric bootstrap bias and variance.
#' @param mrpaw_list The output of one of the Get*MCMCWeights functions
#' @param post The output of `brm(..., survey_df, family=binomial(link="logit"))`
#' @param survey_df The survey dataframe
#'
#' @export
EvalParametricBootstrap <- function(mrpaw_list, post, survey_df) {
    # TODO: fork based on the class so the same function works for both
    stopifnot(class(post) == "brmsfit")
    CheckLogitFamily(post)

    yhat_draws <- SafeGetYhatDraws(mrpaw_list, post, survey_df)

    yhat <- colMeans(yhat_draws)
    y <- GetResponse(post)
    boot_bias <- sum(mrpaw_list$w * (yhat - y))

    vhat <- yhat * (1 - yhat)
    boot_var <- sum(vhat * (mrpaw_list$w^2))
    return(list(bias=boot_bias, var=boot_var))
}



#' Evaluate the nonparametric influence function and IJ variance.
#' @param mrpaw_list The output of one of the Get*MCMCWeights functions
#' @param post The output of `brm(..., survey_df, family=binomial(link="logit"))`
#' @param survey_df The survey dataframe
#'
#' @export
EvalInfluenceFunction <- function(mrpaw_list, post, survey_df) {
    # TODO: fork based on the class so the same function works for both
    stopifnot(class(post) == "brmsfit")
    CheckLogitFamily(post)

    eta_draws <- SafeGetEtaDraws(mrpaw_list, post, survey_df)

    y <- GetResponse(post)
    lp_mat <- (y * t(eta_draws) - log(1 + exp(t(eta_draws)))) %>% t()
    lp_draws <- apply(lp_mat, FUN=sum, MARGIN=1)
    infl_vec <- cov(mrpaw_list$mrp_draws, lp_mat) %>% as.numeric()

    # With N datapoints, Var(N * infl) \approx Var(\sqrt{N} * MrP); the
    # factor of n_obs gives an estimate of Var(MrP).
    n_obs <- length(y)
    ij_var <- n_obs * var(infl_vec)

    return(list(infl_vec=infl_vec, ij_var=ij_var))
}



#' Evaluate the Fisher consistency of the coefficients.
#' Assume the data are drawn from the correct parametric model, and the
#' posterior means are consistent.  Then if we use importance sampling to
#' weight the log likelihood with parameter theta_1 versus theta_0, then
#' the change in the estimate of theta should be (theta_1 - theta_0), and
#' the derivative with respect to this change should be 1.
#'
#' @param mrpaw_list The output of one of the Get*MCMCWeights functions
#' @param post The output of `brm(..., survey_df, family=binomial(link="logit"))`
#' @param survey_df The survey dataframe
#'
#' @export
EvalFisherConsistency <- function(mrpaw_list, post, survey_df) {
    # TODO: fork based on the class so the same function works for both
    stopifnot(class(post) == "brmsfit")
    CheckLogitFamily(post)

    yhat_draws <- SafeGetYhatDraws(mrpaw_list, post, survey_df)
    eta_draws <- SafeGetEtaDraws(mrpaw_list, post, survey_df)

    # Draws of the coefficients
    draw_df <- spread_draws(post, `b_.*`, regex = TRUE)
    x <- model.matrix(as.formula(formula(post)), data=survey_df)

    # Reweight the log likelihood contrbutions with importance sampling
    # weights for theta_1 versus theta_0.
    # The estimate of d E[theta] / d theta_1 | theta_0 is given by
    # cov(theta, \sum_n \partial \log p(y_n | \theta) / \partial \theta)
    y <- GetResponse(post)
    resid_draws <- (y - t(yhat_draws)) %>% t()
    ll_grad_draws <- resid_draws %*% x

    draws_mat <- as_draws_matrix(post)
    draw_colnames <- colnames(x) %>%
        str_replace("\\(", "") %>%
        str_replace("\\)", "") %>%
        paste0("b_", .) 

    stopifnot(all(draw_colnames %in%  colnames(draws_mat)))
    coeff_draws <- draws_mat[, draw_colnames]
    fisher_mat <- cov(coeff_draws, ll_grad_draws)
    return(list(
        ll_grad_draws=ll_grad_draws, 
        draws_mat=draws_mat, 
        fisher_mat=fisher_mat,
        x=x))
}





#' Estimate the coefficient were extra regressors to be included.
#' @param x0 A matrix whose coefficients are assumed to be zero in the original fit
#' @param mrpaw_list The output of one of the Get*MCMCWeights functions
#' @param survey_df The survey dataframe
#'
#' @return The estimated change (implicitly from zero) in the estimated
#' coefficients using a single Newton step at each draw.
#'
#' @export
EvalAlphaHat <- function(x0, mrpaw_list, post, survey_df) {
    # TODO: fork based on the class so the same function works for both
    stopifnot(class(post) == "brmsfit")
    CheckLogitFamily(post)

    n_obs <- nrow(survey_df)
    n_draws <- nrow(yhat_draws)
    stopifnot(nrow(x0) == n_obs)

    yhat_draws <- SafeGetYhatDraws(mrpaw_list, post, survey_df)
    y <- GetResponse(post)
    resid_draws <- (y - t(yhat_draws)) %>% t()
    var_draws <- yhat_draws * (1 - yhat_draws)
    
    EvalLLGradDraws <- einsum_generator("mn,np -> mp")
    EvalLLHessDraws <- einsum_generator("mn,ni,nj -> mij")

    ll_grad_draws <- EvalLLGradDraws(resid_draws, x0)
    ll_hess_draws <- EvalLLHessDraws(var_draws, x0, x0)

    alpha_hat_draws <- sapply(
        1:n_draws, function(m) {
            solve(ll_hess_draws[m,,], ll_grad_draws[m,])
        }) %>% t()
        dim(alpha_hat_draws)
    alpha_hat <- colMeans(alpha_hat_draws)
    names(alpha_hat) <- colnames(x0)
    return(list(alpha_hat=alpha_hat, ll_grad_draws=ll_grad_draws, ll_hess_draws=ll_hess_draws))
}