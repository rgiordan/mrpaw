

# Use pop_w for weights if specfieid, otherwise use 
# a vector of ones as long as pop_df.
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




CheckLogitFamily <- function(logit_fit) {
    logit_family <- family(logit_fit)
    stopifnot(logit_family$family %in% c("binomial", "bernoulli"))
    stopifnot(logit_family$link == "logit")
}





CheckOLSFamily <- function(lin_post) {
    post_family <- family(lin_post)
    stopifnot(post_family$family == "gaussian")
    stopifnot(post_family$link == "identity")
}



#' Get the response variable (y) from the posterior.
#' I don't see this use clearly documented, so I want to factor it out
#' for testing.
GetResponse <- function(post) {
    return(standata(post)$Y)
}

