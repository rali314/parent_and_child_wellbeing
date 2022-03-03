/************************
Abu Dhabi Early Childhood Authority (ECA) Project
Purpose: Clean parent survey variables and create indices
Author: Rahma Ali ra3624@nyu.edu
Date created: 2/Dec/2021
************************/

* Load deidentified data
use processed_data/parent_survey_raw_deid.dta, clear


************************
* Prep control variables
************************
* Parent (respondent) sex indicator
gen male_parent = Sex_parent == 1

* Parent education; Bachelor degree and higher indicator variable
gen ba_higher_education = 0
replace ba_higher_education = 1 if inlist(Q6_corrected_w1, 6,7,8)

* Parent age; impute missing values with mean value
summarize Age_Parent
gen age_parent_imp = Age_Parent
replace age_parent_imp = r(mean) if missing(Age_Parent)

* Region; Abu Dhabi region indicator variable
gen region_ad = Q8_corr_woAlD == 1

* Nationality; UAE nationality indicator variable
gen uae_national = Q508a_w1 == 184

* Job loss because of covid
gen jobloss_covid_ind = 1 if Q27a_w1 == 1
replace jobloss_covid_ind = 0 if jobloss_covid_ind == .


* Generate a variable to capture time gap between wave 1 and wave 2 data collection
gen duration_days = StartDate_w2 - StartDate_w1
gen duration_weeks = duration_days/7
label var duration_days "Number of days between the two data collection interviews"
label var duration_weeks "Number of weeks between the two data collection interviews"

* Label indicator variables 
label define yesno 1 "Yes" 0 "No"
foreach v of varlist male_parent ba_higher_education region_ad ///
uae_national jobloss_covid_ind {
	label val `v' yesno 
}

*********************
* Produce the long data
*********************

* Prepare for reshaping
global reshape_vars ///
Q19_1_w Q19_2_w Q19_3_w Q19_4_w ///
Q18a_1_w Q18a_2_w Q18a_3_w Q18a_4_w ///
Q23_1_w Q23_2_w Q23_3_w ///
Q18a_6_w ///
Q508a_w Q27a_w Q27b_w  Q28a_1_w Q28a_2_w ///
Q37b_1_w Q45a_1_w Q45a_2_w StartDate_w  ///
Q24_1_w Q24_2_w Q24_3_w ///
Q43b_1_w Q15_1_w Q15_3_w ///
Q15_2_w Q18b_8_w Q18b_9_w


* Reshape data from wide to long format
reshape long $reshape_vars, i(ID_obs) j(wave) 

* Remove the _w suffix
rename (*_w) (*)

* Generate wave2 indicator variable
gen wave2 = wave == 2

* Label fixed effects
label var wave2 "Wave 2 of data collection"
label var uae_national "Parent is a UAE national"
label var age_parent_imp "Parent's age - after imputing missing values with the mean"
label var male_parent "Respondent parent gender: male"
label var ba_higher "Respondent parent education: BA degree and higher"
label var region_ad "Region: Abu Dhabi"

*************
* Clean outcome variables
*************
* Health outcomes
*************
egen health_risk = rowmean(Q15_1 Q15_3)
label var health_risk "Health risk index: parents' perceived risk of getting COVID-19 increased over time"

egen parent_stress = rowmean(Q19_1 Q19_2 Q19_3 Q19_4)
label var parent_stress "Perceived parenting stress"

egen parent_wellbeing = rowmean(Q18a_1 Q18a_2 Q18a_3 Q18a_4)
label var parent_wellbeing "Parent's health and well-being index"

* Social isolation
egen social_isolation_p = rowmean(Q23_1 Q23_2 Q23_3)
label var social_isolation_p "Parent's perceived social isolation index"


*label define rating 1 "Poor" 2 "Fair" 3 "Good" 4 "Very good" 5 "Excellent"
foreach v of varlist Q18a_1 Q18a_2 Q18a_3 Q18a_4 ///
{
	label values `v' rating
}


* Keep only relevant variables to the analysis
keep jobloss_covid_ind parent_wellbeing parent_stress health_risk social_isolation_p Q15_2 Q18b_8 Q18b_9 male_parent age_parent_imp uae_national region_ad ba_higher_education ID_obs wave Q502_w2 jobloss


* Reshape back to wide format
reshape wide parent_wellbeing parent_stress health_risk social_isolation_p Q15_2 Q18b_8 Q18b_9 jobloss, i(ID_obs) j(wave)

* Save clean wide data format
save "processed_data/parent_survey_wide_deid_with_indices.dta", replace

export delimited using "processed_data/parent_survey_wide_deid_with_indices.csv", replace








