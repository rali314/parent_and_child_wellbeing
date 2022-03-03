# Load libraries
pac <- c('gridExtra', 'kableExtra', 'xtable', 'MatchIt', 'optmatch', 'tidyverse', 
         'survey', 'rbounds', 'sensitivityfull', 'mice')
lapply(pac, require, character.only = TRUE)

# Set directory
setwd('/Users/Shared/Box/ECA')

# Define vectors
###########################
cov <- c('male_parent', 'age_parent_imp', 'uae_national', 
         'region_ad', 'ba_higher_education', 'jobloss_covid_ind1')
outcomes <- c('parent_wellbeing', 'parent_stress', 'health_risk', 
              'social_isolation_p', 'Q15_2', 'Q18b_8', 'Q18b_9')
w1_cov <- c('parent_stress', 'health_risk', 'Q15_2')
outcomes_names <- c('Health and well-being (1-5)', ' ', 'Perceived parent stress (1-4)', ' ', 'Perceived health risk (%)', ' ', 
                    'Feelings of social isolation (1-5)', ' ', 'Likelihood of economic hardship (%)', ' ', 'Sleep quality (1-5)', ' ', 'Eating habits (1-5)', ' ')
outcomes_names_ns <- c('Health and well-being (1-5)', 'Perceived parent stress (1-4)', 'Perceived health risk (%)', 
                    'Feelings of social isolation (1-5)', 'Likelihood of economic hardship (%)', 'Sleep quality (1-5)', 'Eating habits (1-5)')
outcomes2 <- unlist(lapply(outcomes, function(v) {paste(v, '2', sep = '')}))
cov_match <- c('male_parent', 'age_parent_imp', 'ba_higher_education', 'uae_national')


# Import data
################
df_wide_temp <- read_csv("processed_data/parent_survey_wide_deid_with_indices.csv") %>% 
  mutate(covid_vaccine = case_when(Q502_w2 == 'Yes' ~ 1,
                                   Q502_w2 == 'No' ~ 0),
         ID_obs = as.character(ID_obs)) %>%
  dplyr::filter(!is.na(covid_vaccine))

# Fix missing values first for psych outcomes
###########################
var_vec <- unlist(lapply(outcomes, function(v) {c(paste(v,'1', sep = ''), paste(v,'2', sep = ''))}))
var_vec_w1 <- unlist(lapply(outcomes, function(v) {paste(v,'1', sep = '')}))
var_vec_w2 <- unlist(lapply(outcomes, function(v) {paste(v,'2', sep = '')}))

# first, examine missingness pattern
md.pattern(df_wide_temp)
pMiss <- function(x){sum(is.na(x))/length(x)*100} # see % of missing values in each variable

miss_p_tbl <- as.data.frame(apply(df_wide_temp,2,pMiss)) %>%
  kable("html", digits = 2, caption = "% missing values for psych variables") %>%
  kable_styling(bootstrap_options=c("striped", "hover", "condensed", "responsive"),
                full_width=FALSE)
miss_p_tbl

# use multiple imputation using predictive mean matching method
tempData <- mice(df_wide_temp,m=1,maxit=50,
                 meth='pmm',seed=500)
df_wide_temp <- complete(tempData, 1)

# check one last time for missingness
md.pattern(df_wide_temp[, var_vec])


###############
###############

# Determinants of vaccination
#####################
outcomes1 <- unlist(lapply(outcomes, function(v) {paste(v, '1', sep = '')}))
glm_fmla <- paste('covid_vaccine ~', paste(c(cov, outcomes1), collapse = '+'))

covac_mod <- glm(glm_fmla, family = binomial(), data = df_wide_temp)
m_out <- summary(covac_mod)
or <- as.data.frame(exp(coef(covac_mod)))

tibble('Predictor Variable' = c('Intercept', c('Gender: male', 'Age', 'UAE national', 'Region: Abu Dhabi', 'Education: BA and higher', 'Job loss due to covid'), 
                                outcomes_names_ns),
       # coef = as.data.frame(round(coef(covac_mod), 2)),
       OR = as.data.frame(round(exp(coef(covac_mod)), 2)),
       as.data.frame(cbind(exp(confint.default(covac_mod)))),
       # 'Std. Error' = as.data.frame(round(m_out$coef[,2], 2)),
       'P-value' = as.data.frame(round(m_out$coef[,4], 2)))%>%
  kable("html", digits = 2, caption = " Logistic regression results for factors potentially associated with COVID-19 vaccine") %>%
  kable_styling(bootstrap_options=c("striped", "hover", "condensed", "responsive"),
                full_width=FALSE)

# Checking balance before vaccinated and non 
# vaccinated participants before matching
###################################
cov_chr <- c('male_parent', 'uae_national', 
             'region_ad', 'ba_higher_education', 'jobloss_covid_ind1')

df_wide <- df_wide_temp %>%
  mutate_at(
    vars(one_of(cov_chr)),
    funs(case_when(
      . == 'Yes' ~ 1,
      . == 'No' ~ 0)
    )) %>%
  mutate_at(cov_chr, as.numeric)
bal_vec <- c('age_parent_imp', cov_chr)

covac_balance <- lapply(bal_vec, FUN = function(v) {summary(lm(as.matrix(df_wide[, v]) ~ df_wide$covid_vaccine))})

covac_balance_tbl <- tibble(Variable = c('Age', 'Gender: male', 'UAE national', 'Region: Abu Dhabi', 'Education: BA and higher', 'Job loss due to COVID-19'),
                            Unvaccinated = unlist(lapply(1:length(bal_vec), FUN = function(i) {covac_balance[[i]]$coef[1]})),
                            Vaccinated = unlist(lapply(1:length(bal_vec), FUN = function(i) {covac_balance[[i]]$coef[2]} + covac_balance[[i]]$coef[1])),
                            'P-value' = round(unlist(lapply(1:length(bal_vec), FUN = function(i) {covac_balance[[i]]$coef[8]})), 3),
) %>%
  kable("html", digits = 2, caption = "Demographic differences between vaccinated and unvaccinated participants prior to vaccination") %>%
  kable_styling(bootstrap_options=c("striped", "hover", "condensed", "responsive"),
                full_width=FALSE)

covac_balance_tbl

# Perform matching
#####################
fmps <- paste('covid_vaccine ~', paste(cov_match , collapse ="+"))
mod_match_nn <- matchit(formula(fmps),
                        caliper = 0.25,
                        method = "nearest", discard = 'control', data = df_wide) 
mod_match_nn_out <- summary(mod_match_nn)

# Check matching performance
matching_performance <- cbind(Variable = c('Distance', 'Gender: Male', 'Age', 'Education: BA and higher', 'Nationality: UAE'),
                              as.data.frame(round(mod_match_nn_out[3]$sum.all[,1:3], 2))[3],
                              as.data.frame(mod_match_nn_out[4]$sum.matched[,1:3])[3],
                              row.names = NULL)
colnames(matching_performance) <- c(' ', 'Std mean diff - original sample', 'Std mean diff - matched sample')

matching_performance

matching_performance_tbl <- tibble(matching_performance) %>%
  kable("html", digits = 2, caption = "Standard mean difference in demographic characteristics between vaccinated and unvaccinated participants (samples) at wave 1 (pre-vaccination) before and after the matching") %>% 
  kable_styling(bootstrap_options=c("striped", "hover", "condensed", "responsive"),
                full_width=FALSE)
matching_performance_tbl

# Estimate treatment effects using matched data
##########################################
dta_nn <- match.data(mod_match_nn)


trt_eff_nn1 <- summary(lm(parent_wellbeing2 ~ covid_vaccine+male_parent+age_parent_imp+uae_national+ba_higher_education+parent_wellbeing1+parent_stress1+Q15_21, data = dta_nn))
trt_eff_nn2 <- summary(lm(parent_stress2 ~ covid_vaccine+male_parent+age_parent_imp+uae_national+ba_higher_education+parent_stress1+Q15_21, data = dta_nn))
trt_eff_nn3 <- summary(lm(health_risk2 ~ covid_vaccine+male_parent+age_parent_imp+uae_national+ba_higher_education+health_risk1+parent_stress1+Q15_21, data = dta_nn))
trt_eff_nn4 <- summary(lm(social_isolation_p2 ~ covid_vaccine+male_parent+age_parent_imp+uae_national+ba_higher_education+social_isolation_p1+parent_stress1+Q15_21, data = dta_nn))
trt_eff_nn5 <- summary(lm(Q15_22 ~ covid_vaccine+male_parent+age_parent_imp+uae_national+ba_higher_education+Q15_21+parent_stress1, data = dta_nn))
trt_eff_nn6 <- summary(lm(Q18b_82 ~ covid_vaccine+male_parent+age_parent_imp+uae_national+ba_higher_education+Q18b_81+parent_stress1+Q15_21, data = dta_nn))
trt_eff_nn7 <- summary(lm(Q18b_92 ~ covid_vaccine+male_parent+age_parent_imp+uae_national+ba_higher_education+Q18b_91+parent_stress1+Q15_21, data = dta_nn))


vec_nn <- list(trt_eff_nn1, trt_eff_nn2, trt_eff_nn3, trt_eff_nn4, trt_eff_nn5, trt_eff_nn6, trt_eff_nn7)

trt_effect_nn_tbl <- tibble(Variable = outcomes_names,
                            Unvaccinated = unlist(lapply(1:length(vec_nn), function(v){c(round(vec_nn[[v]]$coefficients[1], 2), '')})),
                            Diff = c(unlist(lapply(c(1), function(v){c(round(vec_nn[[v]]$coefficients[2], 2), paste('(', round(vec_nn[[v]]$coefficients[11], 2), ')', sep = ''))})),
                                     unlist(lapply(c(2), function(v){c(round(vec_nn[[v]]$coefficients[2], 2), paste('(', round(vec_nn[[v]]$coefficients[10], 2), ')', sep = ''))})),
                                     unlist(lapply(c(3,4), function(v){c(round(vec_nn[[v]]$coefficients[2], 2), paste('(', round(vec_nn[[v]]$coefficients[11], 2), ')', sep = ''))})),
                                     unlist(lapply(c(5), function(v){c(round(vec_nn[[v]]$coefficients[2], 2), paste('(', round(vec_nn[[v]]$coefficients[10], 2), ')', sep = ''))})),
                                     unlist(lapply(c(6,7), function(v){c(round(vec_nn[[v]]$coefficients[2], 2), paste('(', round(vec_nn[[v]]$coefficients[11], 2), ')', sep = ''))}))),
                            'P-value' = c(unlist(lapply(c(1), function(v){c(round(vec_nn[[v]]$coefficients[29], 2), '')})), 
                                          unlist(lapply(c(2), function(v){c(round(vec_nn[[v]]$coefficients[26], 2), '')})), 
                                          unlist(lapply(c(3,4), function(v){c(round(vec_nn[[v]]$coefficients[29], 2), '')})), 
                                          unlist(lapply(c(5), function(v){c(round(vec_nn[[v]]$coefficients[26], 2), '')})), 
                                          unlist(lapply(c(6,7), function(v){c(round(vec_nn[[v]]$coefficients[29], 2), '')}))) 
)%>%
  kable("html", digits = 2, caption = "Estimates of COVID-19 vaccination effect on economic, health, and psychosocial factors") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                full_width = FALSE)

trt_effect_nn_tbl

################
################

# Supplementary material
######################
# balance performance comparison
mod_match_opt <- matchit(formula(fmps),
                         method = "optimal", discard = 'none', data = df_wide) 
mod_match_opt_out <- summary(mod_match_opt)

mod_match_nn_nocal <- matchit(formula(fmps),
                              method = "nearest", discard = 'none', data = df_wide) 
mod_match_nn_nocal_out <- summary(mod_match_nn_nocal)

matching_performance2 <- cbind(Variable = c('Distance', 'Gender: male', 'Age', 'Education: BA and higher', 'UAE national'),
                               as.data.frame(round(mod_match_nn_out[3]$sum.all[,1:3], 2)[3]),
                               as.data.frame(mod_match_nn_out[4]$sum.matched[,1:3])[3],
                               as.data.frame(mod_match_nn_nocal_out[4]$sum.matched[,1:3])[3],
                               as.data.frame(mod_match_opt_out[4]$sum.matched[,1:3])[3],
                               row.names = NULL)
colnames(matching_performance) <- c(' ', 'Std mean diff - original sample', 'Std mean diff - matched sample')
matching_performance2
matching_performance_tbl2 <- tibble(matching_performance2) %>%
  kable("html", digits = 2, caption = "Balance comparison before and after matching") %>% 
  kable_styling(bootstrap_options=c("striped", "hover", "condensed", "responsive"),
                full_width=FALSE)
matching_performance_tbl2

# treatment effects using pair match without caliper
dta_nn_nocal <- match.data(mod_match_nn_nocal)
tab_nn_nocal <- lapply(outcomes2, FUN = function(v) {summary(lm(as.matrix(dta_nn_nocal[, v]) ~ covid_vaccine+male_parent+age_parent_imp+uae_national+ba_higher_education+
                                                                  parent_wellbeing1+parent_stress1+health_risk1+social_isolation_p1+Q15_21+Q18b_81+Q18b_91, data = dta_nn_nocal))})

tab_nn_nocal_tbl <- tibble(Variable = outcomes_names,
                           Unvaccinated = unlist(lapply(1:length(outcomes2), FUN = function(i){c(round(tab_nn_nocal[[i]]$coef[1], 2), '')})),
                           Diff = unlist(lapply(1:length(outcomes2), FUN = function(i) {c(round(tab_nn_nocal[[i]]$coef[2],2), paste('(', round(tab_nn_nocal[[i]]$coef[15], 2), ')', sep = ''))})),
                           P_value = unlist(lapply(1:length(outcomes2), FUN = function(i) {c(round(tab_nn_nocal[[i]]$coef[41], 2), '')}))) %>%
  kable("html", digits = 2, caption = "Estimates of COVID-19 vaccination effectusing pair matching with no caliper") %>%
  kable_styling(bootstrap_options=c("striped", "hover", "condensed", "responsive"),
                full_width=FALSE)
tab_nn_nocal_tbl

# treatment effects using optimal pair match
dta_opt <- match.data(mod_match_opt)
tab_opt <- lapply(outcomes2, FUN = function(v) {summary(lm(as.matrix(dta_opt[, v]) ~ covid_vaccine+male_parent+age_parent_imp+uae_national+ba_higher_education+
                                                             parent_wellbeing1+parent_stress1+health_risk1+social_isolation_p1+Q15_21+Q18b_81+Q18b_91, data = dta_opt))})

tab_opt_tbl <- tibble(Variable = outcomes_names,
                      Unvaccinated = unlist(lapply(1:length(outcomes2), FUN = function(i){c(round(tab_opt[[i]]$coef[1], 2), '')})),
                      Diff = unlist(lapply(1:length(outcomes2), FUN = function(i) {c(round(tab_opt[[i]]$coef[2],2), paste('(', round(tab_opt[[i]]$coef[15], 2), ')', sep = ''))})),
                      P_value = unlist(lapply(1:length(outcomes2), FUN = function(i) {c(round(tab_opt[[i]]$coef[41], 2), '')}))) %>%
  kable("html", digits = 2, caption = "Estimates of COVID-19 vaccination effect: using optimal pair matching") %>%
  kable_styling(bootstrap_options=c("striped", "hover", "condensed", "responsive"),
                full_width=FALSE)
tab_opt_tbl


# matching on psych variables at baseline
cov_match2 <- c(cov_match, 'parent_stress1', 'Q15_21')
fmla2 <- paste('covid_vaccine ~', paste(cov_match2 , collapse ="+"))
m_ps2 <- glm(fmla2,
             family = binomial(), data = df_wide)
mod_match_nn2 <- matchit(formula(fmla2),
                         caliper = 0.25,
                         method = "nearest", discard = 'control', data = df_wide) 
mod_match_nn2_out <- summary(mod_match_nn2)

matching_performance3 <- cbind(Variable = c('Distance', 'Gender: male', 'Age', 'Education: BA and higher', 'UAE national', 'Parent stress (1-4)', 'Likelihood of economic hardship (%)'),
                               as.data.frame(round(mod_match_nn2_out[3]$sum.all[,1:3], 2)),
                               as.data.frame(mod_match_nn2_out[4]$sum.matched[,1:3])[3],
                               row.names = NULL)
colnames(matching_performance3) <- c(' ', 'Means treated', 'Means control', 'Std mean diff - original sample', 'Std mean diff - matched sample')
matching_performance3
matching_performance3 <- tibble(matching_performance3) %>%
  kable("html", digits = 2, caption = "Balance comparison before and after matching using psych vars") %>% 
  kable_styling(bootstrap_options=c("striped", "hover", "condensed", "responsive"),
                full_width=FALSE)
matching_performance3

dta_nn2 <- match.data(mod_match_nn2)

trt_eff_nn21 <- summary(lm(parent_wellbeing2 ~ covid_vaccine+male_parent+age_parent_imp+uae_national+ba_higher_education+parent_wellbeing1+parent_stress1+Q15_21, data = dta_nn2))
trt_eff_nn22 <- summary(lm(parent_stress2 ~ covid_vaccine+male_parent+age_parent_imp+uae_national+ba_higher_education+parent_stress1+Q15_21, data = dta_nn2))
trt_eff_nn23 <- summary(lm(health_risk2 ~ covid_vaccine+male_parent+age_parent_imp+uae_national+ba_higher_education+health_risk1+parent_stress1+Q15_21, data = dta_nn2))
trt_eff_nn24 <- summary(lm(social_isolation_p2 ~ covid_vaccine+male_parent+age_parent_imp+uae_national+ba_higher_education+social_isolation_p1+parent_stress1+Q15_21, data = dta_nn2))
trt_eff_nn25 <- summary(lm(Q15_22 ~ covid_vaccine+male_parent+age_parent_imp+uae_national+ba_higher_education+Q15_21+parent_stress1, data = dta_nn2))
trt_eff_nn26 <- summary(lm(Q18b_82 ~ covid_vaccine+male_parent+age_parent_imp+uae_national+ba_higher_education+Q18b_81+parent_stress1+Q15_21, data = dta_nn2))
trt_eff_nn27 <- summary(lm(Q18b_92 ~ covid_vaccine+male_parent+age_parent_imp+uae_national+ba_higher_education+Q18b_91+parent_stress1+Q15_21, data = dta_nn2))

vec_nn2 <- list(trt_eff_nn21, trt_eff_nn22, trt_eff_nn23, trt_eff_nn24, trt_eff_nn25, trt_eff_nn26, trt_eff_nn27)

trt_effect_nn2_tbl <- tibble(Variable = outcomes_names,
                             Unvaccinated = unlist(lapply(1:length(vec_nn2), function(v){c(round(vec_nn2[[v]]$coefficients[1], 2), '')})),
                             Diff = c(unlist(lapply(c(1), function(v){c(round(vec_nn2[[v]]$coefficients[2], 2), paste('(', round(vec_nn2[[v]]$coefficients[11], 2), ')', sep = ''))})),
                                      unlist(lapply(c(2), function(v){c(round(vec_nn2[[v]]$coefficients[2], 2), paste('(', round(vec_nn2[[v]]$coefficients[10], 2), ')', sep = ''))})),
                                      unlist(lapply(c(3,4), function(v){c(round(vec_nn2[[v]]$coefficients[2], 2), paste('(', round(vec_nn2[[v]]$coefficients[11], 2), ')', sep = ''))})),
                                      unlist(lapply(c(5), function(v){c(round(vec_nn2[[v]]$coefficients[2], 2), paste('(', round(vec_nn2[[v]]$coefficients[10], 2), ')', sep = ''))})),
                                      unlist(lapply(c(6,7), function(v){c(round(vec_nn2[[v]]$coefficients[2], 2), paste('(', round(vec_nn2[[v]]$coefficients[11], 2), ')', sep = ''))}))),
                             'P-value' = c(unlist(lapply(c(1), function(v){c(round(vec_nn2[[v]]$coefficients[29], 2), '')})), 
                                           unlist(lapply(c(2), function(v){c(round(vec_nn2[[v]]$coefficients[26], 2), '')})), 
                                           unlist(lapply(c(3,4), function(v){c(round(vec_nn2[[v]]$coefficients[29], 2), '')})), 
                                           unlist(lapply(c(5), function(v){c(round(vec_nn2[[v]]$coefficients[26], 2), '')})), 
                                           unlist(lapply(c(6,7), function(v){c(round(vec_nn2[[v]]$coefficients[29], 2), '')}))) 
)%>%
  kable("html", digits = 2, caption = "Estimates of COVID-19 vaccination effect: Controlling for status at wave 1. Standard errors reported in parentheses") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                full_width = FALSE)

trt_effect_nn2_tbl
