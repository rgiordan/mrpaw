



GetPopulationWeights <- function(pop_df, pop_w=NULL) {
    if (is.null(pop_w)) {
        pop_w <- rep(1, nrow(pop_df)) / nrow(pop_df)
    }

    weight_sum <- sum(pop_w)
    if (abs(weight_sum - 1) > 1e-6) {
        warning(sprintf("The population weights do not sum to one: %f", weight_sum))
    }
    return(pop_w)
}


#' Get MrPaw weights for the OLS estimator.
#' @param lm_fit The output of `lm`
#' @param survey_df The survey dataframe
#' @param pop_df The population dataframe
#'
#' @return The MrP OLS estimate for sanity checking, and the weight vector
#' whose n-th entry is d MrP / d y_n.
#'
#'@export
GetOLSWeights <- function(lm_fit, survey_df, pop_df, pop_w=NULL) {
    stopifnot(class(lm_fit) == "lm")
    pop_w <- GetPopulationWeights(pop_df, pop_w)
    reg_form <- formula(lm_fit)
    x_ols <- model.matrix(reg_form, survey_df)

    # The formula may expect a response called y, so add it in.
    x_pop <- model.matrix(reg_form, pop_df %>% mutate(y=0))

    # The OLS model weights have a closed form.  (The logistic ones do
    # too, but I'll just put this in for now)
    # d \hat\beta / dy_n = (X^T X)^{-1} x_n
    # MrP = \sum_s w_s x_s^T \hat\beta

    betahat <- coefficients(lm_fit)
    yhat_pop <- x_pop %*% betahat
    xtx <- t(x_ols) %*% x_ols
    mrp_ols <- t(pop_w) %*% yhat_pop %>% as.numeric()
    w_ols <- t(pop_w) %*% x_pop %*% solve(xtx, t(x_ols)) %>% as.numeric()
    return(list(mrp=mrp_ols, w=w_ols))
}



CheckLogitFamily <- function(logit_fit) {
    logit_family <- family(logit_fit)
    stopifnot(logit_family$family %in% c("binomial", "bernoulli"))
    stopifnot(logit_family$link == "logit")
}


# TODO: combine this and the above function with a class checker?

#' Get MrPaw weights for the logistic estimator.
#' @param logit_fit The output of `glm(..., survey_df, family=binomial(link="logit"))`
#' @param survey_df The survey dataframe
#' @param pop_df The population dataframe
#'
#' @return The MrP Logistic estimate for sanity checking, and the weight vector
#' whose n-th entry is d MrP / d y_n.
#'
#'@export
GetLogitWeights <- function(logit_fit, survey_df, pop_df, pop_w=NULL) {
    stopifnot(class(logit_fit) == c("glm", "lm"))
    CheckLogitFamily(logit_fit)
    pop_w <- GetPopulationWeights(pop_df, pop_w)

    reg_form <- formula(logit_fit)

    phat <- predict(logit_fit, survey_df, type="response")
    phat_pop <- predict(logit_fit, pop_df, type="response")
    vhat <- phat * (1 - phat)
    vhat_pop <- phat_pop * (1 - phat_pop)

    reg_form <- formula(logit_fit)
    x_ols <- model.matrix(reg_form, survey_df)

    # The formula may expect a response called y, so add it in.
    x_pop <- model.matrix(reg_form, pop_df %>% mutate(y=0))

    # Mrp = w^T expit(x_pop beta)
    # d Mrp / d betahat = w^^T (v . x_pop) 
    # hess = v_survey . x_ols
    # dbetahat_dw <- solve(hess, t(x_ols)) =>
    # d Mrp / d w = w^^T (v . x_pop) hess^{-1} x_survey
    hess <- t(x_ols * vhat) %*% x_ols
    mrp_chain_rule_term <- colSums(pop_w * vhat_pop * x_pop)
    w_logit <- t(mrp_chain_rule_term) %*% solve(hess, t(x_ols)) %>% as.numeric()

    mrp_logit <- sum(pop_w * phat_pop)

    return(list(mrp=mrp_logit, w=w_logit))
}












##############################################
# MCMC

#' Get MrPaw weights for the logistic MCMC estimator.
#' @param logit_post The output of `brm(..., survey_df, family=binomial(link="logit"))`
#' @param survey_df The survey dataframe
#' @param pop_df The population dataframe
#' @param pop_w Optional.  The weight given to each row of pop_df.  Defaults to ones.
#' @param save_preds Optional.  If true, save the posterior predictions for re-use.
#'
#' @return Draws from the MrP estimate, and the weight vector
#' whose n-th entry is d E[MrP | X, Y] / d y_n.
#'
#'@export
GetLogitMCMCWeights <- function(logit_post, survey_df, pop_df, pop_w=NULL, save_preds=FALSE) {
    stopifnot(class(logit_post) == "brmsfit")

    CheckLogitFamily(logit_post)

    pop_w <- GetPopulationWeights(pop_df, pop_w)

    # posterior_epred should be yhat.
    # posterior_linpred should be theta^T x_n.  
    # Draws are in rows and observations in columns.
    yhat_pop <- posterior_epred(logit_post, newdata=pop_df)
    mrp_draws_logit <- yhat_pop %*% pop_w

    # Get the influence scores for the logit model
    # The log likelihood derivative for the n^th datapoint is just the theta^T x_n
    # TODO: optionally return the linpred for further diagnostics
    ll_grad_draws_logit <- posterior_linpred(logit_post, newdata=survey_df)
    w_logit <- cov(mrp_draws_logit, ll_grad_draws_logit)[1,]

    result_list <- list(
        mrp_draws=mrp_draws_logit,
        w=w_logit
    )

    if (save_preds) {
        yhat_draws <- posterior_epred(logit_post, newdata=survey_df)
        eta_draws <- posterior_linpred(logit_post, newdata=survey_df)

        result_list$yhat_pop <- yhat_pop
        result_list$yhat_draws <- yhat_draws
        result_list$eta_draws <- eta_draws
    }
    return(result_list)
}



GetOLSLikelihoodComponentDraws <- function(lin_post, survey_df) {
    # posterior_epred should be yhat.
    # posterior_linpred should be theta^T x_n.  
    # Draws are in rows and observations in columns.

    # get_variables(lin_post)
    sigma_draws <- lin_post %>% spread_draws(sigma) %>% pull(sigma)
    yhat_draws <- posterior_linpred(lin_post, newdata=survey_df)
    stopifnot(ncol(yhat_draws) == nrow(survey_df))
    resid_draws <- (survey_df$y - t(yhat_draws)) %>% t()
    return(list(resid_draws=resid_draws, sigma_draws=sigma_draws, yhat_draws=yhat_draws))
}




CheckOLSFamily <- function(lin_post) {
    post_family <- family(lin_post)
    stopifnot(post_family$family == "gaussian")
    stopifnot(post_family$link == "identity")
}



#' Get MrPaw weights for the logistic MCMC estimator.
#' @param lin_post The output of `brm(..., survey_df, family=gaussian())`
#' @param survey_df The survey dataframe
#' @param pop_df The population dataframe
#' @param pop_w Optional.  The weight given to each row of pop_df.  Defaults to ones.
#'
#' @return Draws from the MrP estimate, and the weight vector
#' whose n-th entry is d E[MrP | X, Y] / d y_n.
#'
#'@export
GetOLSMCMCWeights <- function(lin_post, survey_df, pop_df, pop_w=NULL) {
    stopifnot(class(lin_post) == "brmsfit")
    CheckOLSFamily(lin_post)

    # TODO: here and elsewhere check that "y" is actually the response in the
    # formula, or else extract it.
    stopifnot("y" %in% names(survey_df))
    print("Warning: the response variable must be y for this function to work!")

    pop_w <- GetPopulationWeights(pop_df, pop_w)


    yhat_pop <- posterior_epred(lin_post, newdata=pop_df)

    mrp_draws_lin <- yhat_pop %*% pop_w

    ols_ll_draws <- GetOLSLikelihoodComponentDraws(lin_post, survey_df)

    # The log likelihood derivative for the n^th datapoint is
    # sigma^{-2} (y_n - \hat{y}_n)
    ll_grad_draws_lin <- -1 * ols_ll_draws$resid_draws / (ols_ll_draws$sigma_draws^2)
    w_lin <- cov(mrp_draws_lin, ll_grad_draws_lin)[1,]

    return(list(
        mrp_draws=mrp_draws_lin,
        w=w_lin
    ))
}