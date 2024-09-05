########## MRP example: FOrecasting 2020 elections  #################

# Modified code from https://tellingstorieswithdata.com/ (Chapters: 6, 8, 16)

library(haven)
library(labelled)
library(tidyverse)


# 1. Read in data ---------------------------------------------

# path to natioscape data
# this file can be downloaded from https://www.voterstudygroup.org/data/nationscape 
# Further details in Chapter 8
nationscape_file <- "ns20200625.dta"

# path to ACS data
# downloaded from IPUMS USA (https://usa.ipums.org/usa/index.shtml) as a .dta 
# file. download STATEICP, AGE, SEX, EDUC for the ACS 2019 sample. 
# Further details in Chapter 6.
acs_file <- "usa_00015.dta"

# read in Nationscape data
raw_nationscape_data <- read_dta(nationscape_file)

# read in ACS data
ipums_extract <- read_dta(acs_file)


# 2. Tidy up data -------------------------------------------------------------

# The Stata format separates labels so reunite those
raw_nationscape_data <-
  to_factor(raw_nationscape_data)

# Just keep relevant variables
nationscape_data <-
  raw_nationscape_data |>
  dplyr::select(vote_2020, gender, education, state, age)

nationscape_data <-
  nationscape_data |>
  filter(vote_2020 %in% c("Joe Biden", "Donald Trump")) |>
  mutate(vote_biden = if_else(vote_2020 == "Joe Biden", 1, 0)) |>
  dplyr::select(-vote_2020)

nationscape_data <-
  nationscape_data |>
  mutate(
    age_group = case_when(
      age <= 29 ~ "18-29",
      age <= 44 ~ "30-44",
      age <= 59 ~ "45-59",
      age >= 60 ~ "60+",
      TRUE ~ "Trouble"
    ),
    gender = case_when(
      gender == "Female" ~ "female",
      gender == "Male" ~ "male",
      TRUE ~ "Trouble"
    ),
    education_level = case_when(
      education %in% c(
        "3rd Grade or less",
        "Middle School - Grades 4 - 8",
        "Completed some high school",
        "High school graduate"
      ) ~ "High school or less",
      education %in% c(
        "Other post high school vocational training",
        "Completed some college, but no degree"
      ) ~ "Some post sec",
      education %in% c(
        "Associate Degree",
        "College Degree (such as B.A., B.S.)",
        "Completed some graduate, but no degree"
      ) ~ "Post sec +",
      education %in% c("Masters degree",
                       "Doctorate degree") ~ "Grad degree",
      TRUE ~ "Trouble"
    )
  ) |>
  dplyr::select(-education,-age)


# Format state names to match IPUMS
states_names_and_abbrevs <-
  tibble(stateicp = state.name, state = state.abb)

nationscape_data <-
  nationscape_data |>
  left_join(states_names_and_abbrevs, by = "state")

rm(states_names_and_abbrevs)

# Make lowercase to match IPUMS data
nationscape_data <-
  nationscape_data |>
  mutate(stateicp = tolower(stateicp))

# Replace NAs with DC
nationscape_data$stateicp <-
  replace_na(nationscape_data$stateicp, "district of columbia")

# Tidy the class
survey_df <-
  nationscape_data |>
  mutate(across(c(gender, stateicp, education_level, age_group), 
                as_factor)) |>
  mutate(age_group =
           factor(age_group, levels = c("18-29", "30-44", "45-59", "60+")),
         education_level =
           factor(education_level, levels = c("High school or less", 
                                              "Some post sec", 
                                              "Post sec +", 
                                              "Grad degree")))

saveRDS(survey_df, "survey_df.RDS")


# 3. Calculate counts and proportions from ACS to post-stratify on -------------

ipums_extract <- 
  ipums_extract |>
  dplyr::select(stateicp, sex, age, educd) |>
  to_factor()

cleaned_ipums <-
  ipums_extract |>
  mutate(age = as.numeric(age)) |>
  filter(age >= 18) |>
  rename(gender = sex) |>
  mutate(
    age_group = case_when(
      age <= 29 ~ "18-29",
      age <= 44 ~ "30-44",
      age <= 59 ~ "45-59",
      age >= 60 ~ "60+",
      TRUE ~ "Trouble"
    ),
    education_level = case_when(
      educd %in% c(
        "nursery school, preschool", "kindergarten", "grade 1",
        "grade 2", "grade 3", "grade 4", "grade 5", "grade 6",
        "grade 7", "grade 8", "grade 9", "grade 10", "grade 11",
        "12th grade, no diploma", "regular high school diploma",
        "ged or alternative credential", "no schooling completed"
      ) ~ "High school or less",
      educd %in% c(
        "some college, but less than 1 year",
        "1 or more years of college credit, no degree"
      ) ~ "Some post sec",
      educd  %in% c("associate's degree, type not specified",
                    "bachelor's degree") ~ "Post sec +",
      educd %in% c(
        "master's degree",
        "professional degree beyond a bachelor's degree",
        "doctoral degree"
      ) ~ "Grad degree",
      TRUE ~ "Trouble"
    )
  ) |>
  dplyr::select(gender, age_group, education_level, stateicp) |>
  mutate(across(c(
    gender, stateicp, education_level, age_group),
    as_factor)) |>
  mutate(age_group =
           factor(age_group, levels = c("18-29", "30-44", "45-59", "60+")),
         education_level =
           factor(education_level, levels = c("High school or less", 
                                              "Some post sec", 
                                              "Post sec +", 
                                              "Grad degree")))

poststrat_df <-
  cleaned_ipums |>
  count(stateicp, gender, age_group, education_level) |>
  mutate(prop = n / sum(n))

poststrat_df$stateicp <- factor(poststrat_df$stateicp)
poststrat_df$gender <- factor(poststrat_df$gender)

saveRDS(poststrat_df, "poststrat_df.RDS")
