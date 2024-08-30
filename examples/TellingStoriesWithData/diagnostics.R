source('analysis.R')

# observed outcome 
y_obs <- survey_df$vote_biden  #data dependent
y_obs_factor <- as.factor(y_obs)

n <- nrow(survey_df)

# offset and MrP estimate
cat(sprintf("Offset = %0.3f \t MrP = %0.3f (%0.3f)\n",  
            round(mean(mrp_aw$mrp_draws), 3) - sum(mrp_aw$w * y_obs),
            round(mean(mrp_aw$mrp_draws), 3), 
            round(sd(mrp_aw$mrp_draws), 3)))

# are there negative weights?
cat(sprintf("There are %i negative weights, which are %0.2f%% of the total number of weights", 
            sum(mrp_aw$w < 0),
            round(mean(mrp_aw$w < 0) * 100, 2)))

# weights sum
cat(sprintf("The weights sum to %0.3f", 
            round(sum(mrp_aw$w), 3)))

# dataframe with survey data and weights
weights_df <- cbind(survey_df, w = mrp_aw$w)

# force weights to sum to 1?
normalize = T
if (normalize == T) {weights_df$w <- weights_df$w / sum(weights_df$w)}

# histogram of weights
weights_df %>% 
  ggplot() +
  geom_histogram(aes(w), bins = 100) +
  geom_vline(aes(xintercept = 1 / n)) +
  ggtitle('Weights distribution by outcome') +
  theme_bw()

# density of weights by outcome
weights_df %>% 
  ggplot() +
  geom_density(aes(w, fill = y_obs_factor), alpha = .7) +
  geom_vline(aes(xintercept = 1 / n)) +
  theme_bw()

# comparing the covariate distribution of outliers (defined as respondents with 
# weight > 95% quantile of weight distribution) to the covariate distribution of 
# sample and population

# NOTE: almost all of the outliers belong to the oldest age group and to the 
# lowest educational group which are highly underrepresented in the survey 
# sample in comparison to the population

for (x in cat_regressors) { 
  outliers_df <- rbind(  
    poststrat_df %>%
      group_by(!!sym(x)) %>%
      summarize(n = sum(n)) %>% 
      mutate(tot = sum(n), 
             type = 'population'),
    
    weights_df %>%
      group_by(!!sym(x)) %>%
      count() %>%
      ungroup() %>%
      mutate(tot = sum(n),
             type = 'survey'), 
    
    weights_df %>%
      filter(w > quantile(weights_df$w, probs = .95)) %>%
      group_by(!!sym(x)) %>%
      count() %>%
      ungroup() %>%
      mutate(tot = sum(n),
             type = 'outliers')
  )
  
  print(outliers_df %>%
          mutate(fraction = n/tot) %>% 
          ggplot() +
          geom_point(aes(!!sym(x), fraction, col = type, 
                         shape = type, size = type)) +
          scale_size_manual(values = c(3, 4, 3)) +
          scale_shape_manual(values = c(3, 16, 16)) +
          scale_color_manual(values = c('black', 'red', 'orange')) +
          ggtitle('Fraction of outliers, population, and survey respondents by covariate level') +
          theme_bw())
}


# weights' distribution by covariate levels

for (x in cat_regressors) {
  print(weights_df %>%
          ggplot() +
          geom_boxplot(aes(!!sym(x), w)) +
          ggtitle("Weights' distribution by covariate levels") +
          theme_bw())
}

# covariate balance checks
for (x in cat_regressors) { 
  cov_balance_df <- rbind(
    poststrat_df %>%
      group_by(!!sym(x)) %>%
      summarize(n = sum(n)) %>% 
      mutate(type = 'population',
             Frequency = n / sum(n)) %>%
      dplyr::select(! n),
    
    weights_df %>%
      group_by(!!sym(x)) %>%
      count() %>%
      ungroup() %>%
      mutate(type = 'survey',
             Frequency = n / sum(n)) %>%
      dplyr::select(!n),
    
    weights_df %>%
      group_by(!!sym(x)) %>%
      summarize(type = 'weighted',
                Frequency = sum(w))
  )
  
  print(ggplot(cov_balance_df) +
          geom_point(aes(!!sym(x), Frequency, col = type, 
                         shape = type, size = type)) +
          scale_size_manual(values = c(4, 3, 3)) +
          scale_shape_manual(values = c(16, 16, 3)) +
          scale_color_manual(values = c('red', 'orange', 'black')) +
          ggtitle('Covariate balance checks') +
          theme_bw())
}


# measure of ESS
weights_df %>% 
  summarize(effective_sample_size = sum(abs(w))^2 / sum(w^2),
            sample_size = n())

