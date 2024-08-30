########## MRP example on Gay Policy Responsiveness #################

# Modified code from 
# https://scholar.princeton.edu/sites/default/files/jkastellec/files/mrp_primer_replication_files.zip
# References: 
# https://jkastellec.scholar.princeton.edu/sites/g/files/toruqf3871/files/jkastellec/files/mrp_primer.pdf
# https://www.columbia.edu/~jrl2124/Lax_Phillips_Gay_Policy_Responsiveness_2009.pdf


library("foreign")

# 1. Read in data ---------------------------------------------

#read in megapoll
survey_df <- read.dta("gay_marriage_megapoll.dta", 
                          convert.underscore = TRUE) 

#read in state-level dataset
Statelevel <- read.dta("state_level_update.dta", convert.underscore = TRUE)
Statelevel <- Statelevel[order(Statelevel$sstate.initnum),]

#read in Census data
Census <- read.dta("poststratification 2000.dta", convert.underscore = TRUE)
Census <- Census[order(Census$cstate),]
Census$state.initnum <- match(Census$cstate, Statelevel$sstate)

poststrat_df <- Census %>% 
  rename_with(~ str_replace(., "c", ""), starts_with("c"))

# 2. Tidy up data -------------------------------------------------------------

#Create variables

#At level of megapoll

# from 1 for white males to 6 for hispanic females
survey_df$race.female <- (survey_df$female * 3) + survey_df$race.wbh
# from 1 for 18-29 with low edu to 16 for 65+ with high edu
survey_df$age.edu.cat <- 4 * (survey_df$age.cat -1) + survey_df$edu.cat
# proportion of evangelicals in respondent's state
survey_df$p.evang.full <- Statelevel$p.evang[survey_df$state.initnum]
# proportion of mormon's in respondent's state
survey_df$p.mormon.full <-Statelevel$p.mormon[survey_df$state.initnum]
# combined evangelical + mormom proportions
survey_df$p.relig.full <- survey_df$p.evang.full + survey_df$p.mormon.full
# kerry's % of 2-party vote in respondent's state in 2004
survey_df$p.kerry.full <- Statelevel$kerry.04[survey_df$state.initnum]
survey_df$state[survey_df$state == ""] <- NA
survey_df$region[survey_df$region == ""] <- NA
survey_df <- survey_df %>% 
  drop_na(p.kerry.full) %>%
  mutate(across(contains('cat'), as.factor),
         across(contains('female'), as.factor))


#At census level 

poststrat_df$race.female <- (poststrat_df$female *3) + poststrat_df$race.WBH 
poststrat_df$age.edu.cat <- 4 * (poststrat_df$age.cat-1) + poststrat_df$edu.cat 
poststrat_df$p.evang.full<-  Statelevel$p.evang[poststrat_df$state.initnum]
poststrat_df$p.mormon.full <- Statelevel$p.mormon[poststrat_df$state.initnum]
poststrat_df$p.relig.full <- poststrat_df$p.evang.full + poststrat_df$p.mormon.full
poststrat_df$p.kerry.full <-  Statelevel$kerry.04[poststrat_df$state.initnum]
poststrat_df = poststrat_df %>% 
  rename('n' = .freq, 
         'n.state' = freq.state) %>%
  mutate(prop = n/sum(n),
         across(contains('cat'), as.factor),
         across(contains('female'), as.factor))

saveRDS(survey_df, "survey_df.RDS")
saveRDS(poststrat_df, "poststrat_df.RDS")
