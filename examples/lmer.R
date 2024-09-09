library(tidyverse)
library(lme4)
library(mrpaw)
library(gridExtra)

AEqBLine <- function() {
  geom_abline(aes(slope=1, intercept=0))
}

##########################################
# Simulate some categories
# 

set.seed(25338)

# induce group probabilities with a truncated normal
n_groups <- 2
degree <- 2
n_obs <- 1000
n_obs_pop <- 100000

# Simulate some data.
sim_data <- SimulateSurveyData(n_groups, n_obs, n_obs_pop, degree=degree)

# Survey data:
survey_df <- sim_data$survey_df %>% mutate(x=runif(n()))

# Population data:
pop_df <- sim_data$pop_df %>% mutate(x=runif(n()))

mrp_true <- with(pop_df, mean(ey))
print(mrp_true)

# If there is imbalance, this should differ from the true mrp.
with(survey_df, mean(ey))

# Accumulate data within group.  This helps speed up MCMC prediction.
agg_list <- AggregateSimulationData(sim_data)
survey_agg_df <- agg_list$survey_agg_df
pop_agg_df <- agg_list$pop_agg_df

# The optimal weights are the ratio of the population and sample
# weights within the group.  Get the optimal weight for each row
# of the survey for comparison with implicit weights.
joint_df <- agg_list$joint_df
w_opt <-
  survey_df %>%
  inner_join(select(joint_df, s, w_opt), by="s") %>%
  pull(w_opt)


#####################################
# Run logistic and OLS regression


re_terms <- sprintf("(1 + x | %s)", sim_data$group_effects$g_cols) %>% paste(collapse=" + ")
reg_form <-sprintf( "y ~ 1 +%s", re_terms)

lin_fit <- lmer(formula(reg_form), survey_df)
logit_fit <- glmer(formula(reg_form), survey_df, family=binomial(link="logit"))

reg_terms <- lFormula(formula(reg_form), survey_df)

colnames(reg_terms$X)
names(reg_terms$reTrms$Ztlist)

summary(lin_fit)
class(lin_fit)

# This should be everything you need
fixef(lin_fit) %>% names()
ranef(lin_fit) %>% names()
sigma(lin_fit)
VarCorr(lin_fit) %>% names()

Sigma_g1 <- VarCorr(lin_fit)[["g1"]]
sd_g1 <- attr(Sigma_g1, "stddev")
sqrt(diag(Sigma_g1)) / sd_g1 # Sanity check

getME(lin_fit, "Z") %>% dim()
getME(lin_fit, "Ztlist") %>% names()
getME(lin_fit, "Lambda") %>% dim()


# ??merMod

# From ?mkLmerDevfun
lmod <- lFormula(formula(reg_form), survey_df)
devfun <- do.call(mkLmerDevfun, lmod)
ls(environment(devfun))

# ?optimizeLmer
opt <- optimizeLmer(devfun)

# Inside optimzeLmer
optimizeLmer
rho <- environment(devfun)
rho$pp$theta # Apparently the actual optimization parameter

opt[1:3]
mkMerMod(environment(devfun), opt, lmod$reTrms, fr=lmod$fr)

# Heeeey
# https://www.alexejgossmann.com/Dissect_lmer_part1/
# https://www.alexejgossmann.com/Dissect_lmer_part2/
# https://www.alexejgossmann.com/Dissect_lmer_part3/
# Though this is pretty tied to the REML algorithm, and the mathematical
# details are pretty spare

# This might help
# https://github.com/lme4/lme4pureR/blob/c25bd6ad183ce3959ff6fd66dd7f571456928687/R/pls.R#L107

###################################
# Try to get the point prediction

# y = X beta + \sum_k Z_k gamma_K + eps
reg_terms <- lFormula(formula(reg_form), survey_df)
beta <- fixef(lin_fit)
eps_var <- sigma(lin_fit)^2
gamma_covs <- VarCorr(lin_fit)
gamma_means <- ranef(lin_fit)

# What we want is called "b", which is also lambda u. Here we go:
bhat <- getME(lin_fit, "b")
gamma_means
# unlist(gamma_means) - bhat # Nope

# Now to use our own code we want to map the elements of "b" back
# onto the levels of the random effects.

# maybe this 
# https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf

# Hmmm
# see ?mkReTrms.
# How to map the rows of Zt onto the elements of gamma_means?
reg_terms$reTrms$Zt %>% dim()
reg_terms$reTrms$Ztlist %>% names() # Not helpful



{ 
  # Copied from 
  # getAnywhere(ranef.merMod)
  object <- lin_fit
  lin_fit@pp$b(1)
  class(lin_fit@pp)
  class(lin_fit@pp$b)
  ans <- object@pp$b(1)
  fl <- object@flist
  levs <- lapply(fl, levels)
  asgn <- attr(fl, "assign")
  cnms <- object@cnms
  nc <- lengths(cnms)
  nb <- diff(object@Gp)
  nbseq <- rep.int(seq_along(nb), nb)
  ml <- split(ans, nbseq)
  for (i in seq_along(ml)) {
    ml[[i]] <- matrix(
      ml[[i]], ncol = nc[i],
      byrow = TRUE, dimnames = list(NULL, cnms[[i]]))
  }
  ans <- lapply(seq_along(fl), function(i) {
    m <- ml[asgn == i]
    b2 <- vapply(m, nrow, numeric(1))
    ub2 <- unique(b2)
    if (length(ub2) > 1) 
      stop("differing numbers of b per group")
    rnms <- if (ub2 == length(levs[[i]])) 
      levs[[i]]
    else seq(ub2)
    data.frame(do.call(cbind, m), row.names = rnms, check.names = FALSE)
  })
  names(ans) <- names(fl)
  
  whichel <- names(ans)
  stopifnot(is(whichel, "character"))
  whchL <- names(ans) %in% whichel
  ans <- ans[whchL]
}





reg_terms$reTrms$cnms
reg_terms$reTrms$nl


class(reg_terms)

nrow(survey_df)

stopifnot(names(reg_terms$X) == names(beta))
#stopifnot(names(reg_terms$reTrms$Zt) == names(gamma_means))

beta <- getME(lin_fit, "beta")
bhat <- getME(lin_fit, "b")
Z <- getME(lin_fit, "Z")
X <- getME(lin_fit, "X")
fixef(lin_fit) - beta
y <- survey_df$y
yhat <- X %*% beta + Z %*% bhat
print(mean((yhat - y)^2))

