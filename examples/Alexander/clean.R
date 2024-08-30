########## MRP example using martial name change survey #################

# Modified code from https://github.com/MJAlexander/marriage-name-change
# Reference: https://www.monicaalexander.com/posts/2019-08-07-mrp/ 

library(tidyverse)
library(haven)


# 1. Read in data ---------------------------------------------

# path to MNCS data
# this file can be downloaded here: https://osf.io/uzqdn/
mncs_file <- "MNCS-PV2.dta"

# path to ACS data
# downloaded from IPUMS USA (https://usa.ipums.org/usa/index.shtml) as a .dta 
# file. download AGE, SEX, STATEFIP, MARST, YRMARR, EDUC for the ACS 2017 5 year 
# sample. 
acs_file <- "usa_00003.dta.gz"

# read in MNCS
d <- read_dta(mncs_file)

# read in the ACS data
dc <- read_dta(gzfile(acs_file))  


# 2. Tidy up data -------------------------------------------------------------

# MNCS

# keep variables of interest and only those responses with a year and age

d <- d %>% 
  select(yrmar,
         agemar,
         agemarc,
         genmar,
         spgenmar,
         namechg,
         ednow,
         state) %>% 
  filter(!is.na(agemar), !is.na(yrmar))

# keep only hetero women, make age group, education and decade married variables

d <- d %>% 
  filter(genmar==2&spgenmar==1) %>% 
  mutate(kept_name = as.numeric(namechg==1),
         state_name = str_to_lower(
           as.character(factor(state, 
                               levels = attributes(d$state)$labels, 
                               labels = names(attributes(d$state)$labels)))),
         age = agemar + (2019-yrmar),
         age_group = (as.character(cut(age, 
                                       breaks = c(seq(15, 80, by = 5), Inf),
                                       labels = seq(15, 80, by = 5), 
                                       right = FALSE
         ))),
         decade_married = (as.character(cut(yrmar, 
                                            breaks = c(seq(1969, 2019, by = 10), 
                                                       Inf),
                                            labels = seq(1969, 2019, by = 10), 
                                            right = FALSE
         ))),
         educ_group = case_when(
           ednow<5 ~ "<BA",
           ednow==5 ~ "BA",
           ednow>5 ~ ">BA",
           TRUE ~ "NA"
         )) %>%
  filter(educ_group != "NA" & !is.na(state_name))


# keep only what we want, tidy up and save
# only looked at 25+ year olds to give people a chance to have finished a BA. 
# 80+ and 1969 decade filtered out as these only have one observation

survey_df <- d %>% 
  select(kept_name, state_name, age_group, decade_married, educ_group) %>% 
  filter(age_group>20, age_group<80, decade_married>1969)

saveRDS(survey_df, "survey_df.RDS")

# ACS

# create age groups, decade married and education levels

d_acs <- dc %>% 
  filter(sex==2, age>14, marst!=6, yrmarr>1968) %>% 
  mutate(age_group = (as.character(cut(age, 
                                       breaks = c(seq(15, 80, by = 5), Inf),
                                       labels = seq(15, 80, by = 5), 
                                       right = FALSE
  ))),
  decade_married = (as.character(cut(yrmarr, 
                                     breaks = c(seq(1969, 2019, by = 10), Inf),
                                     labels = seq(1969, 2019, by = 10), 
                                     right = FALSE
  ))),
  educ_group = case_when(
    educ<10 ~ "<BA",
    educ==10~"BA",
    educ==11 ~">BA",
    TRUE~ "NA"
  )) 


# 3. Calculate counts and proportions from ACS to post-stratify on -------------


# calculate cell counts and save

poststrat_df <- d_acs %>% 
  group_by(age_group, statefip, educ_group, decade_married) %>% 
  summarise(n = sum(perwt))%>% 
  mutate(state_name = str_to_lower(
    as.character(factor(statefip, 
                        levels = attributes(d_acs$statefip)$labels, 
                        labels = names(attributes(d_acs$statefip)$labels))))) %>% 
  ungroup() %>% 
  select(-statefip) %>% 
  select(state_name, age_group, decade_married, educ_group, n) %>% 
  filter(age_group>20, age_group<80, decade_married>1969) %>%
  mutate(prop = n / sum(n))

saveRDS(poststrat_df, "poststrat_df.RDS")
