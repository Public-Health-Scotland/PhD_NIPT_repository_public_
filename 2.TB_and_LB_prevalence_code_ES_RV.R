#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# RStudio Workbench is strictly for use by Public Health Scotland staff and     
# authorised users only, and is governed by the Acceptable Usage Policy https://github.com/Public-Health-Scotland/R-Resources/blob/master/posit_workbench_acceptable_use_policy.md.
#
# This is a shared resource and is hosted on a pay-as-you-go cloud computing
# platform.  Your usage will incur direct financial cost to Public Health
# Scotland.  As such, please ensure
#
#   1. that this session is appropriately sized with the minimum number of CPUs
#      and memory required for the size and scale of your analysis;
#   2. the code you write in this script is optimal and only writes out the
#      data required, nothing more.
#   3. you close this session when not in use; idle sessions still cost PHS
#      money!
#
# For further guidance, please see https://github.com/Public-Health-Scotland/R-Resources/blob/master/posit_workbench_best_practice_with_r.md.
#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


# Wed May  28 13:12 2025 ------------------------------

# Authors: Elinor Sebire and Rute Vieira 

## script adapted from original total birth and live birth prevalence scripts for publication in GitHub. To accompany the publication of:
##Total and livebirth prevalence of singleton pregnancies with Down’s Syndrome in Scotland between 2000 and 2021: a population-based study ##


###################### SET UP SCRIPT #################################

######## 1. Load packages ###########################

#library(readxl)      ##for reading in excel files
#library(writexl)     ##for writing to excel files 
library(dplyr)       ##for working with 'tidy' data
#library(stringr)     ##for working with strings
#library(janitor)     ##for rounding 0.5 upwards and clean_names function
#library(magrittr)    ##for %<>% operator
library(ggplot2)
#library(epiDisplay)
#library(reshape2)
library(tidyr)
#library(tidyverse)
library(data.table)
library(epitools)    ## for estimating 95% CIs for prevalence
library(glmmTMB)     ## for modelling the data
library(VGAM)
#library(Hmisc)
#library(MatrixModels)
#library(mvtnorm)
#library(quantreg)
library(rms)         ## for restricted cubic splines 
#library(gridExtra) 


######### 2. Read in the data ########################

SLiCCD_cohort_2 <- readRDS("/PHI_conf/NIPT_evaluation_anon/Data/cohort_2/sliccd_cohort2_2000-2021_anon_with_eurocat_groups.rds")
#updated SLICCD file with eurocat groups and month of birth 
## filepaths for storing data and outputs 

NRS_total_births <-read.csv("/PHI_conf/NIPT_evaluation_anon/Data/cohort_2/NRS_denom_data_final_totalbirths.csv")
## dataset with denominators for total birth prevalence


######## 3. Define filepath where data is held ############

#define filepath where data held
filepath1 <- "/PHI_conf/NIPT_evaluation_anon/Data/"

#e.g. define a working directly
wd_phd <- "/PHI_conf/NIPT_evaluation_anon/Data/wd_phd/"

objective1 <- "/PHI_conf/NIPT_evaluation_anon/Data/wd_phd/Objective1/"

objective2 <- "/PHI_conf/NIPT_evaluation_anon/Data/wd_phd/Objective2/"



################## 4. Format the dataset - extract cohort #####################

# sort dataset by year of preganncy end:###########

sliccd <- SLiCCD_cohort_2 %>%
  arrange(time_period)

# extract only those with singleton pregnancy and DS #######

sliccd_ds_extract <- sliccd %>%
  filter(ds_singleton_cohort == 1)



############## 5. Merge SLiCCD and NRS datasets ##################

#### Prepare SLiCCD dataset

## define age categories in SLICCD dataset
# 6 categories + unknown
sliccd_ds_extract <- sliccd_ds_extract %>%
  mutate(age_category = cut(maternal_age_at_conception,
                            breaks = c(-Inf, 20, 25, 30, 35, 40, Inf),
                            labels = c('<20', '20 - 24', '25-29', '30-34', '35-39', '40+'),
                            include.lowest = TRUE,
                            right = FALSE),
         age_category = ifelse(is.na(age_category), 'Unknown', as.character(age_category)))

# 5 categories  + unknown (merging <20 and 20-24 to equal <25)
sliccd_ds_extract <- sliccd_ds_extract%>%
  mutate(age_category = recode(age_category,
                               "<20" = "<25",
                               "20 - 24" = "<25"))

## recode infant sex in sliccd (1 = male, 2 = female, 3 = unknown)
sliccd_ds_extract <- sliccd_ds_extract %>%
  mutate(sex_2 = case_when(sex == 1 ~ 'Male',
                           sex == 2 ~ 'Female',
                           sex == 3 ~ 'Unknown',
                           TRUE ~ NA_character_))

## sort out structure and values for specific variables in the dataset - characters/factors and 'unknown' groups 
sliccd_ds_extract$age_category <- as.factor(sliccd_ds_extract$age_category)
sliccd_ds_extract$sex_2 <- as.factor(sliccd_ds_extract$sex_2)
sliccd_ds_extract$pregnancy_end_type <- as.factor(sliccd_ds_extract$pregnancy_end_type)

## aggregate (count) values by year, age category, simd quintile, and health board to calculate numerator for prevalences
sliccd_numerator_data <- sliccd_ds_extract%>%
  group_by(time_period, age_category, maternal_simd_quintile, healthboard_of_residence)%>%
  summarise(Count = n())

## replace NAs with Unknown in SIMD and Healthboard columns

# SIMD
sliccd_numerator_data$maternal_simd_quintile <- as.character(sliccd_numerator_data$maternal_simd_quintile) 
sliccd_numerator_data <- sliccd_numerator_data%>%
  mutate( maternal_simd_quintile = replace_na(maternal_simd_quintile, "Unknown"))
sliccd_numerator_data$maternal_simd_quintile <- as.factor(sliccd_numerator_data$maternal_simd_quintile)

# Healthboard 
sliccd_numerator_data <- sliccd_numerator_data%>%
  mutate(healthboard_of_residence = replace_na(healthboard_of_residence, "Unknown"))
sliccd_numerator_data$healthboard_of_residence<- as.factor(sliccd_numerator_data$healthboard_of_residence)


#### Prepare NRS dataset - denominator data for total births

# recode SIMD
NRS_total_births <- NRS_total_births%>%
  mutate(maternal_simd_quintile = recode(maternal_simd_quintile,
                                         "1 (most deprived)" = "1",
                                         "5 (least deprived)" = "5"))

### sort out structure and values for specific variables in the dataset - characters/factors 
NRS_total_births$age_category <- as.factor(NRS_total_births$age_category)
NRS_total_births$healthboard_of_residence <- as.factor(NRS_total_births$healthboard_of_residence)
NRS_total_births$maternal_simd_quintile <- as.factor(NRS_total_births$maternal_simd_quintile)


#### Merge SLiCCD and NRS  datasets
merged_data_total_births <- NRS_total_births%>%
  full_join(sliccd_numerator_data, by = c("time_period", "age_category", "maternal_simd_quintile", "healthboard_of_residence"))

# rename the count columns 
merged_data_total_births <- merged_data_total_births%>%
  rename(count_denom = Count.x,
         count_sliccd = Count.y
  )



############## 6. Descriptive statistics - total birth prevalences ##################

## total births in both sliccd and denom data
total_counts <- merged_data_total_births%>%
  summarise(total_denom = sum(count_denom, na.rm = TRUE),
            total_sliccd = sum(count_sliccd, na.rm = TRUE))

## crude total birth prevalence (= livebirths + stillbirths + terminations / total births (stillbirths + livebirths) in NRS Scotland data x 10000)
total_counts$total_sliccd/total_counts$total_denom*10000

# create function to calculate the crude prevalence and 95 CI
calculate_CI<- function(total_sliccd , total_denom){
  prevalence <- (total_sliccd / total_denom)* 10000
  ci <- pois.exact(total_sliccd, pt =  total_denom, conf.level = 0.95) * 10000
  ci_lower <- (ci$lower / total_denom)* 10000
  ci_upper <- (ci$upper / total_denom)* 10000
  list(prevalence = prevalence, ci_lower = ci$lower, ci_upper = ci$upper)
}

# CI for crude total prevalence
ci_totalPrev <- calculate_CI(total_counts$total_sliccd,total_counts$total_denom)

## crude total birth prev (95% CI) per 10000 per year in dataset
# total counts per year 
total_counts_per_year <- merged_data_total_births%>%
  group_by(time_period)%>%
  summarise(total_denom = sum(count_denom, na.rm = TRUE),
            total_sliccd = sum(count_sliccd, na.rm = TRUE))

# prevalence and 95% CI per 10000 per year
prevalence_data_per_year <- total_counts_per_year%>%
  rowwise()%>%
  mutate(ci = list(calculate_CI(total_sliccd, total_denom)),
         prevalence = ci$prevalence,
         ci_lower = ci$ci_lower,
         ci_upper = ci$ci_upper,
         births_prev_ci = paste0(round(prevalence,2), " (", round(ci_lower, 2), "," ,round(ci_upper, 2), ")")) %>%
  ungroup()
prevalence_data_per_year

# plotting prev data per year

prevalence_data_per_year%>%
  ggplot(aes(x = factor(time_period), y = prevalence, group = 1))+  
  geom_rect(aes(xmin = '2000', xmax = '2012', ymin = 0, ymax = 30),#to use for final models? ADDS IN INDICTAION OF WHEN NIPT was implemented - too busy for all graphs
            fill = "cadetblue1", alpha = 0.02)+
  annotate("text", x = '2005', y = 28, label = 'PRR = 1.00 (0.99, 1,02), p = 0.99', vjust = 1, hjust = 0.25, size = 5)+ 
  geom_rect(aes(xmin = '2013', xmax = '2015', ymin = 0, ymax = 30),#to use for final models? ADDS IN INDICTAION OF WHEN NIPT was implemented - too busy for all graphs
            fill = "cadetblue2", alpha = 0.02)+
  annotate("text", x = '2013', y = 6, label = 'PRR = 1.16 (1.07, 1.26), p <0.001', vjust = 1, hjust = 0.25, size = 5)+
  geom_rect(aes(xmin = '2016', xmax = '2021', ymin = 0, ymax = 30),#to use for final models? ADDS IN INDICTAION OF WHEN NIPT was implemented - too busy for all graphs
            fill = "cadetblue3", alpha = 0.02)+
  annotate("text", x = '2018', y = 28, label = 'PRR = 0.94 (0.92, 0.97), p <0.001', vjust = 1, hjust = 0.25, size = 5)+
  geom_point(colour = "blue4")+
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, colour = "blue4")+
  geom_line(colour = "blue4")+  
  geom_line(aes(y = mean(prevalence), group = 1), colour = "red")+
  labs(x = "Year (of pregnancy end)", y = "Total birth prevalence, per 10,000 total births (95% CI)")+
  ##add shaded area to highlight the implemenation of NIPT
  
  theme_minimal()+
  # scale_x_continuous(limits = c(2000,2021), breaks = c(2000:2021))
  ylim(0,30)+
  theme(text = element_text(size=20))


############## 7. Modelling total birth prevalences time trends ##################

## Transforming/centering year (time_period) as glmmTMB needs it to start in 0
prevalence_data_per_year <- transform(prevalence_data_per_year, time_period = time_period - min(time_period))


##### 1. Total birth prevalence (= livebirths + stillbirths + terminations / total births (stillbirths + livebirths) in NRS Scotland data x 10000)


#### a. assuming linear time trend ######
birth_prev_linear1 <- glmmTMB(total_sliccd ~ time_period + offset(log(total_denom)),
                              family = compois, data = prevalence_data_per_year)     
summary(birth_prev_linear1)

# obtain prevalence RR and 95% CI
round(exp(confint(birth_prev_linear1)), 3)

## plot linear trend for crude prevalence + 95% CI
newx <- seq(min(prevalence_data_per_year$time_period), max(prevalence_data_per_year$time_period), by = 1)

# Predict values
births_COMpoissonm1 <- data.frame(time_period = newx, total_denom = prevalence_data_per_year$total_denom)
predicted_values <- predict(birth_prev_linear1, newdata = births_COMpoissonm1, type = "link", se.fit = TRUE)

# Extract predicted values and standard errors
predicted_means <- predicted_values$fit
standard_errors <- predicted_values$se.fit

# Calculate confidence intervals
z_value <- qnorm(0.975)  # 1.96 for a 95% confidence interval
lower_bound <- exp(predicted_means - z_value * standard_errors)
upper_bound <- exp(predicted_means + z_value * standard_errors)

# display results
result_df_COM1 <- data.frame(
  time_period = births_COMpoissonm1$time_period,
  Predicted_Mean = (exp(predicted_means) / prevalence_data_per_year$total_denom) * 10000,
  Lower_Bound = (lower_bound / prevalence_data_per_year$total_denom) * 10000,
  Upper_Bound = (upper_bound / prevalence_data_per_year$total_denom) * 10000,
  birth_prev = (prevalence_data_per_year$total_sliccd / prevalence_data_per_year$total_denom) * 10000  # Include prevs in the result_df
)
print(result_df_COM1)



# Figure 1
prevalence_data_per_year$birth_prev <- prevalence_data_per_year$prevalence

scatterplot <- ggplot(prevalence_data_per_year, aes(x = time_period + 2000, y = birth_prev)) +
  geom_point(aes(y = birth_prev), shape = 1, size = 2, color = "black") 

# add fitted values and confidence intervals as blue lines
lineplot <- scatterplot +
  geom_line(data = result_df_COM1, aes(x = time_period+2000, y = Predicted_Mean), color = "blue") +
  geom_ribbon(data = result_df_COM1, aes(x = time_period+2000, ymin = Lower_Bound, ymax = Upper_Bound), fill = "blue", alpha = 0.2) +
  
  # customize plot labels and appearance
  labs(title = "Linear time trend",
       x = "Year",
       y = "Prevalence (per 10,000 total births)") +
  scale_x_continuous(
    breaks = seq(2000,2021, by = 1),
    labels = seq (2000, 2021, by = 1))+
  theme_minimal()+
  theme(text = element_text(size=20))


#### b. assuming non-linear time trend ################

##### 3 knots (RCS) ################
birth_prev_nonlin1 <- glmmTMB(total_sliccd ~ rcs(time_period,3) + offset(log(total_denom)),
                              family = compois, data = prevalence_data_per_year)
summary(birth_prev_nonlin1)

##### 4 knots (RCS) ###############
birth_prev_nonlin1.2 <- glmmTMB(total_sliccd ~ rcs(time_period,4) + offset(log(total_denom)),
                                family = compois, data = prevalence_data_per_year)
summary(birth_prev_nonlin1.2)

##### 5 knots (RCS) ###############
birth_prev_nonlin1.3 <- glmmTMB(total_sliccd ~ rcs(time_period,5) + offset(log(total_denom)),
                                family = compois, data = prevalence_data_per_year)
summary(birth_prev_nonlin1.3)

##### 4 knots but specifiying the knot placement ########

knots <- c(0,10,16,21)
birth_prev_nonlin1.4 <- glmmTMB(total_sliccd ~ rcs(time_period,knots) + offset(log(total_denom)),
                                family = compois, data = prevalence_data_per_year)
summary(birth_prev_nonlin1.4)

knots2 <- c(0,13,16,21)

birth_prev_nonlin1.5 <- glmmTMB(total_sliccd ~ rcs(time_period,knots2) + offset(log(total_denom)),
                                family = compois, data = prevalence_data_per_year)
summary(birth_prev_nonlin1.5)


#// Best model was non-linear with 4 knots (0,13,16,21). A piecewise regression will be performed to obtain the pRRs and 95% CIs //#


### Piecewise Regression Model

### SUB SECTION OF DATA FROM 2000-2012 (time period 1) and running a Poisson regression with linear trend to look at the trend in this period only.

period_1 <- prevalence_data_per_year%>%
  filter(between(time_period, 0, 12))

period_2 <- prevalence_data_per_year%>%
  filter(between(time_period, 13, 16))

period_3 <- prevalence_data_per_year%>%
  filter(between(time_period, 17, 21))

###### run separate linear regressions for these periods:

#period 1
PW_period1 <- glmmTMB(total_sliccd ~ time_period + offset(log(total_denom)),
                      family = compois, data = period_1) 
summary(PW_period1)
round(exp(confint(PW_period1)), 3)

#period 2
PW_period2 <- glmmTMB(total_sliccd ~ time_period + offset(log(total_denom)),
                      family = compois, data = period_2)
summary(PW_period2)
round(exp(confint(PW_period2)), 3)

#period 3
PW_period3 <- glmmTMB(total_sliccd ~ time_period + offset(log(total_denom)),
                      family = compois, data = period_3)
summary(PW_period3)
round(exp(confint(PW_period3)), 3)


## Figure S1 -- specified 4 knots (0,13,16,21)

# Predict values
births_COMpoissonm1 <- data.frame(time_period = newx, total_denom = prevalence_data_per_year$total_denom)
predicted_values_nl <- predict(birth_prev_nonlin1.5, newdata = births_COMpoissonm1, type = "link", se.fit = TRUE)

# Extract predicted values and standard errors
predicted_means_nl <- predicted_values_nl$fit
standard_errors_nl <- predicted_values_nl$se.fit

# Calculate confidence intervals
z_value <- qnorm(0.975)  # 1.96 for a 95% confidence interval
lower_bound_nl <- exp(predicted_means_nl - z_value * standard_errors_nl)
upper_bound_nl <- exp(predicted_means_nl + z_value * standard_errors_nl)

# Display results
result_df_COM_nl <- data.frame(
  time_period = births_COMpoissonm1$time_period,
  Predicted_Mean = (exp(predicted_means_nl) / prevalence_data_per_year$total_denom) * 10000,
  Lower_Bound = (lower_bound_nl / prevalence_data_per_year$total_denom) * 10000,
  Upper_Bound = (upper_bound_nl / prevalence_data_per_year$total_denom) * 10000,
  birth_prev = (prevalence_data_per_year$total_sliccd / prevalence_data_per_year$total_denom) * 10000  # Include prevs in the result_df
)
print(result_df_COM_nl)

# Plot
prevalence_data_per_year$birth_prev <- (prevalence_data_per_year$total_sliccd / prevalence_data_per_year$total_denom) * 10000

scatterplot <- ggplot(prevalence_data_per_year, aes(x = time_period + 2000, y = birth_prev)) +
  geom_point(aes(y = birth_prev), shape = 1, size = 2, color = "black") 

# Add fitted values and confidence intervals as blue lines
lineplot <- scatterplot +
  geom_line(data = result_df_COM_nl, aes(x = time_period+2000, y = Predicted_Mean), color = "blue") +
  geom_ribbon(data = result_df_COM_nl, aes(x = time_period+2000, ymin = Lower_Bound, ymax = Upper_Bound), fill = "blue", alpha = 0.2) +
  
  # Customize plot labels and appearance
  labs(title = "Non-linear time trend - 4 knots (specified 2000, 2013, 2016, 2021)",
       x = "Year",
       y = "Prevalence (per 10,000 total births)") +
  scale_x_continuous(
    breaks = seq(2000,2021, by = 1),
    labels = seq (2000, 2021, by = 1))+
  theme_minimal()+
  theme(text = element_text(size=20))


############## 8. Modelling total birth prevalences by maternal age ###########

## breakdown counts for sliccd and NRS by age are found in merged_data_total_births
# creating dataset with just maternal age group and prevalence per age group 
MA_group_data <- merged_data_total_births%>%
  group_by(age_category)%>%
  summarise(count_sliccd = sum(count_sliccd, na.rm = TRUE),
            count_denom = sum(count_denom, na.rm = TRUE))%>%
  ungroup()

# removing unknown category 
MA_group_data <- MA_group_data%>%
  filter(!age_category == 'Unknown')

# add total birth prevalence column to dataset
MA_group_data <- MA_group_data%>%
  mutate(totalbirth_prevalence = count_sliccd/count_denom * 10000)

str(MA_group_data)

# plot crude prevalence by age category 
ggplot(MA_group_data, aes(x=age_category, y = totalbirth_prevalence))+
  geom_line(group =1, colour = "purple")+
  geom_point(stat="identity", fill = "purple")+
  labs(title = "The total birth prevalence of babies with DS in Scotland by maternal age groups between 2000-2021", x = "Maternal age group", y = "Total birth prevalence (per 10,000 births)")+
  theme_minimal()



##### model prevalence by age model - VGLM package (no time) #########################
vglm1<-vglm(MA_group_data$count_sliccd ~ MA_group_data$age_category, poissonff(bred=TRUE), offset=log(count_denom), 
            data=MA_group_data) #when working with just categorical predictor have used vglm

## summary of model 
vglm1summary <- summary(vglm1)
show.summary.vglm(vglm1summary)

# CIs
CI <- data.frame(confint(vglm1))

## create summary table of model
# dataframe 
vglm1summary_tbl <- data.frame (coef(vglm1))

## adding CIs 
vglm1summary_tbl <- vglm1summary_tbl%>%
  mutate(Lower = CI$X2.5.., Upper = CI$X97.5.. )

## rename table columns 
colnames(vglm1summary_tbl)[1] <- "Coefficients_Est"

## add in coefficients names column 
vglm1summary_tbl <- vglm1summary_tbl%>%
  mutate(Coefficients = c("Intercept" ,"25-29", "30-34", "35-39", "40+"))

# reorder columns in data
vglm1summary_tbl <- vglm1summary_tbl[c("Coefficients", "Coefficients_Est", "Lower", "Upper")]

# add in PRR for each coefficient and CI
vglm1summary_tbl <- vglm1summary_tbl%>%
  mutate(pRR = exp(coef(vglm1)), pRR_CI_lower = exp(Lower), PRR_CI_upper = exp(Upper))

# add in p values from model
vglm1summary_tbl<- vglm1summary_tbl%>%
  mutate("P value" = coef(summary(vglm1))[,'Pr(>|z|)'])

##rounding all to 2 s.f.
vglm1summary_tbl$pRR_CI_lower <- round(vglm1summary_tbl$pRR_CI_lower, 2)
vglm1summary_tbl$pRR_CI_upper <- round(vglm1summary_tbl$PRR_CI_upper, 2)
vglm1summary_tbl$pRR <-round(vglm1summary_tbl$pRR,2)
vglm1summary_tbl$`P value` <-round(vglm1summary_tbl$`P value`,2)
vglm1summary_tbl$Coefficients_Est <-round(vglm1summary_tbl$Coefficients_Est,2)
vglm1summary_tbl$Lower <-round(vglm1summary_tbl$Lower,2)
vglm1summary_tbl$Upper <-round(vglm1summary_tbl$Upper,2)
vglm1summary_tbl <- vglm1summary_tbl%>%
  mutate(estimate_CI = paste0(Coefficients_Est," (", Lower, ", ", Upper,")"))%>%
  mutate(pRR_CI = paste0(pRR," (", vglm1summary_tbl$pRR_CI_lower, ", ", pRR_CI_upper,")"))

## AIC of model (for info - not reported)
AIC(vglm1)


## plot MA models predicted prevalence-------------------------------------

## For birth prevalence 

# Predict values
age_vglm1 <- data.frame(age_category = MA_group_data$age_category, count_sliccd = MA_group_data$count_sliccd, count_denom = MA_group_data$count_denom)
predicted_values <- predict(vglm1, newdata = age_vglm1, type = "link", se.fit = TRUE)#predicted values for each age group in model

# Extract predicted values and standard errors
predicted_means <- predicted_values$fit
standard_errors <- predicted_values$se.fit

# Calculate confidence intervals
z_value <- qnorm(0.975)  # 1.96 for a 95% confidence interval
lower_bound <- exp(predicted_means - z_value * standard_errors)
upper_bound <- exp(predicted_means + z_value * standard_errors)

# Display results
result_df_vglm1 <- data.frame(
  age_category = age_vglm1$age_category,
  Predicted_Mean = (exp(predicted_means) / MA_group_data$count_denom) * 10000,
  Lower_Bound = (lower_bound / MA_group_data$count_denom) * 10000,
  Upper_Bound = (upper_bound / MA_group_data$count_denom) * 10000,
  birth_prev = (MA_group_data$count_sliccd / MA_group_data$count_denom) * 10000  # Include prevs in the result_df
)
print(result_df_vglm1)


######## plot this data - predicted TB prev per MA 
ggplot(result_df_vglm1, aes(x = age_category, y = birth_prev)) +
  geom_line(aes(y = Predicted_Mean, group = 1), colour = "purple")+
  geom_errorbar(aes(ymin = Lower_Bound, ymax = Upper_Bound), colour = "purple", width = 0.4)+
  geom_point(aes(y = birth_prev), shape = 19, size = 2, alpha = 0.7, color = "black")+ 
  labs(title = 'Observed and predicted total birth prevalence of pregnancies with DS, by maternal age group',#this is plotting predicted from model, and observed actual TB prevalence.
       x = "Maternal age group (yrs)", y = "Total birth prevalence (per 10,000 births")+
  scale_y_continuous(breaks = seq(0, 125, by = 20))+
  theme(text = element_text(size=25))+
  theme_minimal()

### final MA model prevalence plot for manuscript (Figure 3a) ####

MA_plot1 <- ggplot(result_df_vglm1, aes(x = age_category, y = birth_prev)) +
  geom_line(aes(y = birth_prev, group = 1), colour = "steel blue")+
  geom_errorbar(aes(ymin = Lower_Bound, ymax = Upper_Bound), colour = "steel blue", width = 0.6)+
  geom_point(aes(y = birth_prev), shape = 19, size = 3 ,color = "steel blue")+  
  labs(x = "Maternal age group (yrs)", y = "Total birth prevalence per 10,000 total births", title = "a")+
  theme_minimal()+
  ylim(0,120)+
  theme(
    axis.title.x = element_text(margin = margin(t = 15)),  # Add space above x-axis title
    axis.title.y = element_text(margin = margin(r = 15))   # Add space to the right of y-axis title
  )+
  theme(text = element_text(size=20))

MA_plot1



##### model prevalence by age and time #########################

### creating time and MA dataset to use #########
MA_group_data_bytime <- merged_data_total_births%>%
  group_by(time_period, age_category)%>%
  summarise(count_sliccd = sum(count_sliccd, na.rm = TRUE),
            count_denom = sum(count_denom, na.rm = TRUE))%>%
  ungroup()

# removing unknown category 
MA_group_data_bytime <- MA_group_data_bytime%>%
  filter(!age_category == 'Unknown')

# add total birth prevalence column to dataset
MA_group_data_bytime <- MA_group_data_bytime%>%
  mutate(totalbirth_prevalence = count_sliccd/count_denom * 10000)

# change time to 0-21
MA_group_data_bytime <- transform(MA_group_data_bytime, time_period = time_period - min(time_period))

# maternal age adjusted by (4 knots) time trend 
knots2 <- c(0,13,16,21)
model_with_time <- glmmTMB(count_sliccd ~ age_category + rcs(time_period,knots2) + offset(log(count_denom)),
                           family = compois, data = MA_group_data_bytime) 
summary (model_with_time)
model_time_summ<- summary(model_with_time)
MA_model_with_time <- data.frame(model_time_summ$coefficients$cond)

# table of model output MA + TIME -----------------------------
setDT(MA_model_with_time, keep.rownames = TRUE)
MA_model_with_time$Estimate<- round(MA_model_with_time$Estimate, 3)
MA_model_with_time$Pr...z..<- round(MA_model_with_time$Pr...z.., 3)

# remove unneeded columns 
MA_model_with_time$Std..Error <- NULL
MA_model_with_time$z.value <- NULL

# confidence intervals in dataset, rounded to 2 dp
CI_model_with_time_summary <- data.frame(round(confint(model_with_time), 2))
CI_model_with_time_summary <- CI_model_with_time_summary%>%
  mutate(pRR = exp(Estimate), CI_lower = exp(X2.5..), CI_upper = exp(X97.5..))
CI_model_with_time_summary$pRR <- round(CI_model_with_time_summary$pRR, 2)
CI_model_with_time_summary$CI_lower <- round(CI_model_with_time_summary$CI_lower, 2)
CI_model_with_time_summary$CI_upper <- round(CI_model_with_time_summary$CI_upper, 2)
CI_model_with_time_summary <- CI_model_with_time_summary%>%
  mutate(estimate_CI = paste0(Estimate," (", X2.5.., ", ", X97.5..,")"))%>%
  mutate(pRR_CI = paste0(pRR," (", CI_lower, ", ", CI_upper,")"))

# add in CI column to overall summary table 
MA_model_with_time <- MA_model_with_time%>%
  mutate(Estimate_CI = CI_model_with_time_summary$estimate_CI)%>%
  mutate(pRR_CI = CI_model_with_time_summary$pRR_CI)

### plot observed and predicted prev per year - MA + time ----------------------------------------

# Predict values
age_time <- data.frame(time_period = MA_group_data_bytime$time_period, age_category = MA_group_data_bytime$age_category, count_sliccd = MA_group_data_bytime$count_sliccd, count_denom = MA_group_data_bytime$count_denom)
predicted_values <- predict(model_with_time, newdata = age_time, type = "link", se.fit = TRUE)

# Extract predicted values and standard errors
predicted_means <- predicted_values$fit
standard_errors <- predicted_values$se.fit

# Calculate confidence intervals
z_value <- qnorm(0.975)  # 1.96 for a 95% confidence interval
lower_bound <- exp(predicted_means - z_value * standard_errors)
upper_bound <- exp(predicted_means + z_value * standard_errors)

# Display results
result_df_vglm1 <- data.frame(
  time_period = age_time$time_period,
  age_category = age_time$age_category,
  Predicted_Mean = (exp(predicted_means) / MA_group_data_bytime$count_denom) * 10000,
  Lower_Bound = (lower_bound / MA_group_data_bytime$count_denom) * 10000,
  Upper_Bound = (upper_bound / MA_group_data_bytime$count_denom) * 10000,
  birth_prev = (MA_group_data_bytime$count_sliccd / MA_group_data_bytime$count_denom) * 10000  # Include prevs in the result_df
)
print(result_df_vglm1)

MA_group_data_bytime <- MA_group_data_bytime%>%
  mutate(predicted_mean = result_df_vglm1$Predicted_Mean, lower_bound = result_df_vglm1$Lower_Bound, upper_bound = result_df_vglm1$Upper_Bound)

## final MA + time plot (Figure 4a)

MA_plot2 <- ggplot(MA_group_data_bytime, aes(x = time_period+2000, y = totalbirth_prevalence, group = age_category))+
  geom_point(aes(colour = age_category),size=3)+
  geom_line(alpha = 0.2)+
  geom_line(aes(y = predicted_mean, colour = age_category))+
  geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound,fill = age_category), alpha = 0.3 )+
  
  labs(x = "Year",
       y = "Total birth prevalence per 10,000 total births", title = "a",
       colour = "Maternal age group", fill = "Maternal age group")+
  ylim(0,150)+
  theme_minimal()+
  theme(
    axis.title.x = element_text(margin = margin(t = 15)),  # Add space above x-axis title
    axis.title.y = element_text(margin = margin(r = 15))   # Add space to the right of y-axis title
  )+
  theme(text = element_text(size=20))

MA_plot2



####### 9. SIMD quintile model --------------------------------------------------

# create dataset - SIMD data in SLiCCD
SIMD_group_data <- merged_data_total_births%>%
  group_by(maternal_simd_quintile)%>%
  summarise(count_sliccd = sum(count_sliccd, na.rm = TRUE),
            count_denom = sum(count_denom, na.rm = TRUE))%>%
  ungroup()

# removing unknown category 
SIMD_group_data <- SIMD_group_data%>%
  filter(!maternal_simd_quintile == 'Unknown')

# Add total birth prevalence column to dataset 
SIMD_group_data <- SIMD_group_data%>%
  mutate(totalbirth_prevalence = count_sliccd/count_denom * 10000)
str(SIMD_group_data)

## plot crude prevalence by SIMD group --------------------------------------

ggplot(SIMD_group_data, aes(x=maternal_simd_quintile, y = totalbirth_prevalence))+
  geom_line(group =1, colour = "purple")+
  geom_point(stat="identity", fill = "purple")+
  labs(title = "The total birth prevalence of babies with DS in Scotland by maternal SIMD quintile", x = "Maternal SIMD group", y = "Total birth prevalence (per 10,000 total births)")+
  theme_minimal()+
  ylim(0,30)

dev.off()

## DS prev by SIMD COMpoisson model - vglm ---------------------------------
vglm1<-vglm(SIMD_group_data$count_sliccd ~ SIMD_group_data$maternal_simd_quintile, poissonff(bred=TRUE), offset=log(count_denom), 
            data=SIMD_group_data)

# summary of model
vglm1summary <- summary(vglm1)
show.summary.vglm(vglm1summary)

# CIs
CI <- data.frame(confint(vglm1))

# create summary table
vglm1summary_tbl <- data.frame (coef(vglm1))
vglm1summary_tbl <- vglm1summary_tbl%>%
  mutate(Lower = CI$X2.5.., Upper = CI$X97.5.. )

# rename table columns 
colnames(vglm1summary_tbl)[1] <- "Coefficients_Est"

# add in coefficients names column 
vglm1summary_tbl <- vglm1summary_tbl%>%
  mutate(Coefficients = c("Intercept" ,"2", "3", "4", "5"))

# reorder columns in data
vglm1summary_tbl <- vglm1summary_tbl[c("Coefficients", "Coefficients_Est", "Lower", "Upper")]

# add in pRR for each coefficient and CI
vglm1summary_tbl <- vglm1summary_tbl%>%
  mutate(pRR = exp(coef(vglm1)), pRR_CI_lower = exp(Lower), PRR_CI_upper = exp(Upper))

# add in p values
vglm1summary_tbl<- vglm1summary_tbl%>%
  mutate("P value" = coef(summary(vglm1))[,'Pr(>|z|)'])

# rounding to 2sf
vglm1summary_tbl$pRR_CI_lower <- round(vglm1summary_tbl$pRR_CI_lower, 2)
vglm1summary_tbl$pRR_CI_upper <- round(vglm1summary_tbl$PRR_CI_upper, 2)
vglm1summary_tbl$pRR <-round(vglm1summary_tbl$pRR,2)
vglm1summary_tbl$`P value` <-round(vglm1summary_tbl$`P value`,2)
vglm1summary_tbl$Coefficients_Est <-round(vglm1summary_tbl$Coefficients_Est,2)
vglm1summary_tbl$Lower <-round(vglm1summary_tbl$Lower,2)
vglm1summary_tbl$Upper <-round(vglm1summary_tbl$Upper,2)

# final table
vglm1summary_tbl <- vglm1summary_tbl%>%
  mutate(estimate_CI = paste0(Coefficients_Est," (", Lower, ", ", Upper,")"))%>%
  mutate(pRR_CI = paste0(pRR," (", vglm1summary_tbl$pRR_CI_lower, ", ", pRR_CI_upper,")"))


### plot predicted DS prevalence - by SIMD --------------------
## For total birth prevalence 
# Predict values
SIMD_vglm1 <- data.frame(maternal_simd_quintile = SIMD_group_data$maternal_simd_quintile, count_sliccd = SIMD_group_data$count_sliccd, count_denom = SIMD_group_data$count_denom)
predicted_values <- predict(vglm1, newdata = SIMD_vglm1, type = "link", se.fit = TRUE) #predicted values per SIMD group

# Extract predicted values and standard errors
predicted_means <- predicted_values$fit
standard_errors <- predicted_values$se.fit

# Calculate confidence intervals
z_value <- qnorm(0.975)  # 1.96 for a 95% confidence interval
lower_bound <- exp(predicted_means - z_value * standard_errors)
upper_bound <- exp(predicted_means + z_value * standard_errors)

# Display results
result_df_vglm1 <- data.frame(
  maternal_simd_quintile = SIMD_vglm1$maternal_simd_quintile,
  Predicted_Mean = (exp(predicted_means) / SIMD_group_data$count_denom) * 10000,
  Lower_Bound = (lower_bound / SIMD_group_data$count_denom) * 10000,
  Upper_Bound = (upper_bound / SIMD_group_data$count_denom) * 10000,
  birth_prev = (SIMD_group_data$count_sliccd / SIMD_group_data$count_denom) * 10000  # Include prevs in the result_df
)
print(result_df_vglm1)

ggplot(result_df_vglm1, aes(x = maternal_simd_quintile, y = birth_prev)) +
  geom_line(aes(y = birth_prev, group = 1), colour = "purple")+
  geom_errorbar(aes(ymin = Lower_Bound, ymax = Upper_Bound), colour = "purple", width = 0.4)+
  geom_point(aes(y = birth_prev), shape = 19, size = 2, color = "black", alpha = 0.7)+
  labs(x = "Maternal SIMD group", y = "Total birth prevalence per 10,000 total births")+
  ylim(0,30)+
  theme_minimal()+
  theme(text = element_text(size=20))

## final unadjusted SIMD plot - Figure 3b ----------------------------------------------------

SIMD_PLOT1<- ggplot(result_df_vglm1, aes(x = maternal_simd_quintile, y = birth_prev)) +
  geom_line(aes(y = birth_prev, group = 1), colour = "steel blue")+
  geom_errorbar(aes(ymin = Lower_Bound, ymax = Upper_Bound), colour = "steel blue", width = 0.4)+
  geom_point(aes(y = birth_prev), shape = 19, size = 3, color = "steel blue")+
  labs(x = "Maternal SIMD group (1 - most deprived, 5 - least deprived)", y = "Total birth prevalence per 10,000 total births", title = "b")+
  ylim(0,30)+
  theme_minimal()+
  theme(
    axis.title.x = element_text(margin = margin(t = 15)),  # Add space above x-axis title
    axis.title.y = element_text(margin = margin(r = 15))   # Add space to the right of y-axis title
  )+
  theme(text = element_text(size=20))

SIMD_PLOT1

## SIMD + time  -----------------------------------------------------------
# creating time and SIMD dataset to use
SIMD_group_data_bytime <- merged_data_total_births%>%
  group_by(time_period, maternal_simd_quintile)%>%
  summarise(count_sliccd = sum(count_sliccd, na.rm = TRUE),
            count_denom = sum(count_denom, na.rm = TRUE))%>%
  ungroup()

# removing unknown category 
SIMD_group_data_bytime <- SIMD_group_data_bytime%>%
  filter(!maternal_simd_quintile == 'Unknown')

# Add total birth prevalence column to dataset
SIMD_group_data_bytime <- SIMD_group_data_bytime%>%
  mutate(totalbirth_prevalence = count_sliccd/count_denom * 10000)

# change time to 0-21:
SIMD_group_data_bytime <- transform(SIMD_group_data_bytime, time_period = time_period - min(time_period))

# SIMD adjusted by (4 knots) time trend 
model_with_time <- glmmTMB(count_sliccd ~ maternal_simd_quintile + rcs(time_period, knots2) + offset(log(count_denom)),
                           family = compois, data = SIMD_group_data_bytime) 

# summary of model and table 
summary(model_with_time)
model_time_summ<- summary(model_with_time)
SIMD_model_with_time <- data.frame(model_time_summ$coefficients$cond)
setDT(SIMD_model_with_time, keep.rownames = TRUE)
SIMD_model_with_time$Estimate<- round(SIMD_model_with_time$Estimate, 2)
SIMD_model_with_time$Pr...z..<- round(SIMD_model_with_time$Pr...z.., 2)

# remove unneeded columns 
SIMD_model_with_time$Std..Error <- NULL
SIMD_model_with_time$z.value <- NULL

# CIs into dataframe
CI_model_with_time_summary <- data.frame(round(confint(model_with_time), 2))
CI_model_with_time_summary <- CI_model_with_time_summary%>%
  mutate(pRR = exp(Estimate), CI_lower = exp(X2.5..), CI_upper = exp(X97.5..))
CI_model_with_time_summary$pRR <- round(CI_model_with_time_summary$pRR, 2)
CI_model_with_time_summary$CI_lower <- round(CI_model_with_time_summary$CI_lower, 2)
CI_model_with_time_summary$CI_upper <- round(CI_model_with_time_summary$CI_upper, 2)

## final table
CI_model_with_time_summary <- CI_model_with_time_summary%>%
  mutate(estimate_CI = paste0(Estimate," (", X2.5.., ", ", X97.5..,")"))%>%
  mutate(pRR_CI = paste0(pRR," (", CI_lower, ", ", CI_upper,")"))

###add in CI column to summary table 
SIMD_model_with_time <- SIMD_model_with_time%>%
  mutate(Estimate_CI = CI_model_with_time_summary$estimate_CI)%>%
  mutate(pRR_CI = CI_model_with_time_summary$pRR_CI)

######### plot change in prevalence over time with SIMD - crude ----------------

ggplot(SIMD_group_data_bytime, aes(x = time_period+2000, y = totalbirth_prevalence, color = maternal_simd_quintile, group = maternal_simd_quintile))+
  geom_line(size=0.5)+
  geom_point(size=1)+
  labs(title = "Total birth prevalence of singleton pregnancies with DS, per maternal SIMD group over time",
       x = "Year (of pregnancy end)",
       y = "Prevalence of DS (per 10,000 total births)",
       colour = "Maternal SIMD group")+
  theme_minimal()


#### plot of predicted prev per year - Figure 4b

# Predict values
SIMD_time <- data.frame(time_period = SIMD_group_data_bytime$time_period, maternal_simd_quintile = SIMD_group_data_bytime$maternal_simd_quintile, count_sliccd = SIMD_group_data_bytime$count_sliccd, count_denom = SIMD_group_data_bytime$count_denom)
predicted_values <- predict(model_with_time, newdata = SIMD_time, type = "link", se.fit = TRUE)

# Extract predicted values and standard errors
predicted_means <- predicted_values$fit
standard_errors <- predicted_values$se.fit

# Calculate confidence intervals
z_value <- qnorm(0.975)  # 1.96 for a 95% confidence interval
lower_bound <- exp(predicted_means - z_value * standard_errors)
upper_bound <- exp(predicted_means + z_value * standard_errors)

# Display results
result_df_vglm1 <- data.frame(
  time_period = SIMD_time$time_period,
  maternal_simd_quintile = SIMD_time$maternal_simd_quintile,
  Predicted_Mean = (exp(predicted_means) / SIMD_group_data_bytime$count_denom) * 10000,
  Lower_Bound = (lower_bound / SIMD_group_data_bytime$count_denom) * 10000,
  Upper_Bound = (upper_bound / SIMD_group_data_bytime$count_denom) * 10000,
  birth_prev = (SIMD_group_data_bytime$count_sliccd / SIMD_group_data_bytime$count_denom) * 10000  # Include prevs in the result_df
)
print(result_df_vglm1)

SIMD_group_data_bytime <- SIMD_group_data_bytime%>%
  mutate(predicted_mean = result_df_vglm1$Predicted_Mean, lower_bound = result_df_vglm1$Lower_Bound, upper_bound = result_df_vglm1$Upper_Bound)

ggplot(SIMD_group_data_bytime, aes(x = time_period, y = totalbirth_prevalence, group = maternal_simd_quintile))+
  geom_line(aes(colour = maternal_simd_quintile, y = predicted_mean))+
  geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound, fill = maternal_simd_quintile), alpha = 0.3)+
  geom_line(aes(colour = maternal_simd_quintile), alpha = 0.3, 
            size=0.5)+
  geom_point(aes(colour = maternal_simd_quintile), size=3)+
  labs(title = "Observed and predicted total birth prevalence of singleton pregnancies with DS per SIMD group, adjusted for overall time trend",
       x = "Year (of pregnancy end)",
       y = "Total birth prevalence of DS (per 10,000 total births)",
       colour = "Maternal SIMD group",
       fill = "Maternal SIMD group")+
  theme_minimal()+
  theme(text = element_text(size=15))

### final SIMD plot + time for manuscript ---------------------------------------

SIMD_plot2 <- ggplot(SIMD_group_data_bytime, aes(x = time_period+2000, y = totalbirth_prevalence, group = maternal_simd_quintile))+
  geom_line(aes(colour = maternal_simd_quintile, y = predicted_mean))+
  geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound, fill = maternal_simd_quintile), alpha = 0.3)+
  geom_line(aes(colour = maternal_simd_quintile), alpha = 0.3, 
            size=0.5)+
  geom_point(aes(colour = maternal_simd_quintile), size=3)+
  labs(x = "Year",
       y = "Total birth prevalence per 10,000 total births",
       title = "b",
       colour = "Maternal SIMD group",
       fill = "Maternal SIMD group")+
  ylim(0,40)+
  theme_minimal()+
  theme(
    axis.title.x = element_text(margin = margin(t = 15)),  # Add space above x-axis title
    axis.title.y = element_text(margin = margin(r = 15))   # Add space to the right of y-axis title
  )+
  theme(text = element_text(size=20))

SIMD_plot2



################## 10. Health board of residence---------------------------------


### data set - health board data in SLiCCD

HB_group_data <- merged_data_total_births%>%
  group_by(healthboard_of_residence)%>%
  summarise(count_sliccd = sum(count_sliccd, na.rm = TRUE),
            count_denom = sum(count_denom, na.rm = TRUE))%>%
  ungroup()

HB_group_data$healthboard_of_residence<- as.factor(HB_group_data$healthboard_of_residence)
###removing unknown category 

HB_group_data <- HB_group_data%>%
  filter(!healthboard_of_residence == 'Unknown')

##### ADD TOTAL BIRTH PREVALENCE COLUMN TO DATASET 
HB_group_data <- HB_group_data%>%
  mutate(totalbirth_prevalence = count_sliccd/count_denom * 10000)

## grouping island boards into islands label 

HB_group_data<- HB_group_data%>%
  mutate(healthboard_of_residence = if_else(healthboard_of_residence %in% c('NHS Western Isles', 'NHS Orkney', 'NHS Shetland'),
                                            'NHS Island Boards',
                                            healthboard_of_residence))

HB_group_data<- HB_group_data%>%
  group_by(healthboard_of_residence)%>%
  summarise(count_sliccd = sum(count_sliccd), count_denom = sum(count_denom), .groups = 'drop')%>%
  mutate(totalbirth_prevalence = count_sliccd/ count_denom * 10000)%>%
  ungroup()

## crude CIs calculated for each health board

HB_group_data <- HB_group_data%>%
  rowwise()%>%
  mutate(ci = list(calculate_CI(count_sliccd, count_denom)),
         prevalence = ci$prevalence,
         ci_lower = ci$ci_lower,
         ci_upper = ci$ci_upper,
         births_prev_ci = paste0(round(prevalence,2), " (", round(ci_lower, 2), "," ,round(ci_upper, 2), ")")) %>%
  ungroup()


### re-structuring hb to a factor variable

HB_group_data$healthboard_of_residence<- as.factor(HB_group_data$healthboard_of_residence)
str(HB_group_data)


## Unadjusted HB model------------------------------------------------
HB_group_data$healthboard_of_residence = relevel(HB_group_data$healthboard_of_residence, ref = "NHS Ayrshire and Arran")#reference HB is ayrshire and arran

vglm1<-vglm(HB_group_data$count_sliccd ~ HB_group_data$healthboard_of_residence + offset(log(count_denom)), poissonff(bred= TRUE),
            data=HB_group_data)

##summary of model
vglm1summary <- summary(vglm1)
show.summary.vglm(vglm1summary)
vglm1summary <- summaryvglm(vglm1)

###### creating summary table with PRR and CIs --------------------------

##CIs
CI <- data.frame(confint(vglm1))

###create summary table

vglm1summary_tbl <- data.frame (coef(vglm1))
vglm1summary_tbl <- vglm1summary_tbl%>%
  mutate(Lower = CI$X2.5.., Upper = CI$X97.5.. )

## rename table columns 
colnames(vglm1summary_tbl)[1] <- "Coefficients_Est"

## add in coefficients names column 
HB_group_data$healthboard_of_residence

vglm1summary_tbl <- vglm1summary_tbl%>%
  mutate(Coefficients = c("Intercept" ,     "NHS Borders"     ,              "NHS Dumfries and Galloway"    , "NHS Fife"    ,                 
                          "NHS Forth Valley"  ,            "NHS Grampian"      ,            "NHS Greater Glasgow and Clyde", "NHS Highland"  ,  "NHS Island Boards"   ,          
                          "NHS Lanarkshire"    ,           "NHS Lothian"    ,  "NHS Tayside"))

#reorder columns in data
vglm1summary_tbl <- vglm1summary_tbl[c("Coefficients", "Coefficients_Est", "Lower", "Upper")]

##add in pRR for each coefficient and CI
vglm1summary_tbl <- vglm1summary_tbl%>%
  mutate(pRR = exp(coef(vglm1)), pRR_CI_lower = exp(Lower), PRR_CI_upper = exp(Upper))

### add in p values
vglm1summary_tbl<- vglm1summary_tbl%>%
  mutate("P value" = coef(summary(vglm1))[,'Pr(>|z|)'])

## rounding to 2sf

vglm1summary_tbl$pRR_CI_lower <- round(vglm1summary_tbl$pRR_CI_lower, 2)
vglm1summary_tbl$pRR_CI_upper <- round(vglm1summary_tbl$PRR_CI_upper, 2)
vglm1summary_tbl$pRR <-round(vglm1summary_tbl$pRR,2)
vglm1summary_tbl$`P value` <-round(vglm1summary_tbl$`P value`,3)
vglm1summary_tbl$Coefficients_Est <-round(vglm1summary_tbl$Coefficients_Est,2)
vglm1summary_tbl$Lower <-round(vglm1summary_tbl$Lower,2)
vglm1summary_tbl$Upper <-round(vglm1summary_tbl$Upper,2)

##final table

vglm1summary_tbl <- vglm1summary_tbl%>%
  mutate(estimate_CI = paste0(Coefficients_Est," (", Lower, ", ", Upper,")"))%>%
  mutate(pRR_CI = paste0(pRR," (", vglm1summary_tbl$pRR_CI_lower, ", ", pRR_CI_upper,")"))



####### plotting predicted TB prevalence in graph - by HB (Figure) --------------------

## For total birth prevalence 

# Predict values
HB_vglm1 <- data.frame(healthboard_of_residence = HB_group_data$healthboard_of_residence, count_sliccd = HB_group_data$count_sliccd, count_denom = HB_group_data$count_denom)
predicted_values <- predict(vglm1, newdata = HB_vglm1, type = "link", se.fit = TRUE)

# Extract predicted values and standard errors
predicted_means <- predicted_values$fit
standard_errors <- predicted_values$se.fit

# Calculate confidence intervals
z_value <- qnorm(0.975)  # 1.96 for a 95% confidence interval
lower_bound <- exp(predicted_means - z_value * standard_errors)
upper_bound <- exp(predicted_means + z_value * standard_errors)

# Display results
result_df_vglm1 <- data.frame(
  healthboard_of_residence = HB_vglm1$healthboard_of_residence,
  Predicted_Mean = (exp(predicted_means) / HB_group_data$count_denom) * 10000,
  Lower_Bound = (lower_bound / HB_group_data$count_denom) * 10000,
  Upper_Bound = (upper_bound / HB_group_data$count_denom) * 10000,
  birth_prev = (HB_group_data$count_sliccd / HB_group_data$count_denom) * 10000  # Include prevs in the result_df
)
print(result_df_vglm1)


##### final TB prevalence of DS model for HB of residence (Figure 3c) ------------------------

HB_plot1 <- ggplot(result_df_vglm1, aes(x = healthboard_of_residence, y = birth_prev)) +
  geom_bar(stat = "identity", fill = 'steel blue', alpha = 0.3)+
  geom_point(aes(y = birth_prev), shape = 19, size = 3, color = "steel blue")+
  geom_errorbar(aes(ymin = Lower_Bound, ymax = Upper_Bound), colour = "steel blue", width = 0.4)+
  labs(x = "Maternal NHS health board of residence", y = "Total birth prevalence per 10,000 total births", title = "c")+
  ylim(0,30)+
  coord_flip()+
  theme_minimal()+
  theme(
    axis.title.x = element_text(margin = margin(t = 15)),  # Add space above x-axis title
    axis.title.y = element_text(margin = margin(r = 15))   # Add space to the right of y-axis title
  )+
  theme(text = element_text(size=20))

HB_plot1


### 11. healthboard + time model -----------------------------------------

### creating time and HB dataset to use 
HB_group_data_bytime <- merged_data_total_births%>%
  group_by(time_period, healthboard_of_residence)%>%
  summarise(count_sliccd = sum(count_sliccd, na.rm = TRUE),
            count_denom = sum(count_denom, na.rm = TRUE))%>%
  ungroup()

###removing unknown category 
HB_group_data_bytime <- HB_group_data_bytime%>%
  filter(!healthboard_of_residence == 'Unknown')

##### ADD TOTAL BIRTH PREVALENCE COLUMN TO DATASET 
HB_group_data_bytime <- HB_group_data_bytime%>%
  mutate(totalbirth_prevalence = count_sliccd/count_denom * 10000)



### group island boards:

###### grouping island boards into islands label 
HB_group_data_bytime<- HB_group_data_bytime%>%
  mutate(healthboard_of_residence = if_else(healthboard_of_residence %in% c('NHS Western Isles', 'NHS Orkney', 'NHS Shetland'),
                                            'NHS Island Boards',
                                            healthboard_of_residence))

## full dataset for TB prev HB and time:
HB_group_data_bytime<- HB_group_data_bytime%>%
  group_by(time_period, healthboard_of_residence)%>%
  summarise(count_sliccd = sum(count_sliccd), count_denom = sum(count_denom), .groups = 'drop')%>%
  mutate(totalbirth_prevalence = count_sliccd/ count_denom * 10000)

########### plot total birth prev by year and health board 
ggplot(HB_group_data_bytime, aes(x = factor(time_period), y = totalbirth_prevalence, color = healthboard_of_residence, group = healthboard_of_residence))+
  geom_line(size=0.5)+
  geom_point(size=1)+
  labs(title = "Total birth prevalence of singleton pregnancies with DS, per NHS healthboard of residence over time (1)",
       x = "Year (of pregnancy end)",
       y = "Prevalence of DS (per 10,000 total births)",
       colour = "NHS healthboard of residence")+
  facet_wrap(healthboard_of_residence~., ncol = 1)+
  theme_minimal()

#### change time to 0-21:
HB_group_data_bytime <- transform(HB_group_data_bytime, time_period = time_period - min(time_period))

#### HB adjusted by (4 knots) time trend 
model_with_time <- glmmTMB(count_sliccd ~ healthboard_of_residence + rcs(time_period,knots2) +offset(log(count_denom)),
                           family = compois, data = HB_group_data_bytime)

## summary of the model 
summary (model_with_time)
model_time_summ<- summary(model_with_time)


## dataset of the summary
HB_model_with_time <- data.frame(model_time_summ$coefficients$cond)
setDT(HB_model_with_time, keep.rownames = TRUE)
HB_model_with_time$Estimate<- round(HB_model_with_time$Estimate, 5)
HB_model_with_time$Pr...z..<- round(HB_model_with_time$Pr...z.., 3)

###remove unneeded columns 

HB_model_with_time$Std..Error <- NULL
HB_model_with_time$z.value <- NULL

## CIs and rounded values
CI_model_with_time_summary <- data.frame(round(confint(model_with_time), 3))
CI_model_with_time_summary <- CI_model_with_time_summary%>%
  mutate(pRR = exp(Estimate), CI_lower = exp(X2.5..), CI_upper = exp(X97.5..))
CI_model_with_time_summary$pRR <- round(CI_model_with_time_summary$pRR, 3)
CI_model_with_time_summary$CI_lower <- round(CI_model_with_time_summary$CI_lower, 5)
CI_model_with_time_summary$CI_upper <- round(CI_model_with_time_summary$CI_upper, 5)
CI_model_with_time_summary <- CI_model_with_time_summary%>%
  mutate(estimate_CI = paste0(Estimate," (", X2.5.., ", ", X97.5..,")"))%>%
  mutate(pRR_CI = paste0(pRR," (", CI_lower, ", ", CI_upper,")"))

###add in CI column to summary table 

HB_model_with_time <- HB_model_with_time%>%
  mutate(Estimate_CI = CI_model_with_time_summary$estimate_CI)%>%
  mutate(pRR_CI = CI_model_with_time_summary$pRR_CI)

#### calculating predicted prev per year by health board ---------------------

# Predict values
HB_time <- data.frame(time_period = HB_group_data_bytime$time_period, healthboard_of_residence = HB_group_data_bytime$healthboard_of_residence, count_sliccd = HB_group_data_bytime$count_sliccd, count_denom = HB_group_data_bytime$count_denom)
predicted_values <- predict(model_with_time, newdata = HB_time, type = "link", se.fit = TRUE)

# Extract predicted values and standard errors
predicted_means <- predicted_values$fit
standard_errors <- predicted_values$se.fit

# Calculate confidence intervals
z_value <- qnorm(0.975)  # 1.96 for a 95% confidence interval
lower_bound <- exp(predicted_means - z_value * standard_errors)
upper_bound <- exp(predicted_means + z_value * standard_errors)

# Display results
result_df_vglm1 <- data.frame(
  time_period = HB_time$time_period,
  healthboard_of_residence = HB_time$healthboard_of_residence,
  Predicted_Mean = (exp(predicted_means) / HB_group_data_bytime$count_denom) * 10000,
  Lower_Bound = (lower_bound / HB_group_data_bytime$count_denom) * 10000,
  Upper_Bound = (upper_bound / HB_group_data_bytime$count_denom) * 10000,
  birth_prev = (HB_group_data_bytime$count_sliccd / HB_group_data_bytime$count_denom) * 10000  # Include prevs in the result_df
)
print(result_df_vglm1)

HB_group_data_bytime <- HB_group_data_bytime%>%
  mutate(predicted_mean = result_df_vglm1$Predicted_Mean, lower_bound = result_df_vglm1$Lower_Bound, upper_bound = result_df_vglm1$Upper_Bound)

#### plot HB data predicted by time (Figure 4c) ---------------------------------------------
HB_plot3 <- ggplot(HB_group_data_bytime, aes(x = time_period+2000, y = totalbirth_prevalence, group = healthboard_of_residence))+
  geom_line(aes(colour = healthboard_of_residence, y = predicted_mean))+
  geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound, fill = healthboard_of_residence), alpha = 0.3)+
  geom_line(aes(colour = healthboard_of_residence), alpha = 0.3, 
            size=0.5)+
  geom_point(aes(colour = healthboard_of_residence), size=3)+
  labs(x = "Year",
       y = "Total birth prevalence per 10,000 total births",
       title = "a",
       colour = "Maternal NHS health board of residence",
       fill = "Maternal NHS health board of residence")+
  facet_wrap(healthboard_of_residence~., ncol = 1)+
  theme_minimal()+
  theme(
    axis.title.x = element_text(margin = margin(t = 15)),  # Add space above x-axis title
    axis.title.y = element_text(margin = margin(r = 15))   # Add space to the right of y-axis title
  )+
  theme(text = element_text(size=20))

HB_plot3

### 12. Full model ----------------------------------------------------------

## model to include time (non-linear), MA, SIMD and healthboard 
## using original dataset with full breakdown

#### change time to 0-21:
merged_data_total_births <- transform(merged_data_total_births, time_period = time_period - min(time_period))

## remove unknown category for each variable 
final_model_data <- merged_data_total_births%>%
  group_by(time_period, age_category, maternal_simd_quintile, healthboard_of_residence)%>%
  summarise(count_sliccd = sum(count_sliccd, na.rm = TRUE),
            count_denom = sum(count_denom, na.rm = TRUE))%>%
  ungroup()

##removing unknown 
final_model_data <- final_model_data%>%
  filter(!age_category == 'Unknown')

final_model_data <- final_model_data%>%
  filter(!healthboard_of_residence == 'Unknown')

final_model_data <- final_model_data%>%
  filter(!maternal_simd_quintile == 'Unknown')

###### grouping island boards into islands label 
final_model_data<- final_model_data%>%
  mutate(healthboard_of_residence = if_else(healthboard_of_residence %in% c('NHS Western Isles', 'NHS Orkney', 'NHS Shetland'),
                                            'NHS Island Boards',
                                            healthboard_of_residence))

## final model dataset 
final_model_data<- final_model_data%>%
  group_by(time_period, age_category,maternal_simd_quintile ,healthboard_of_residence)%>%
  summarise(count_sliccd = sum(count_sliccd, na.rm = TRUE), count_denom = sum(count_denom, na.rm = TRUE), .groups = 'drop')%>%
  mutate(totalbirth_prevalence = count_sliccd/ count_denom * 10000)

### running full model with significant factors, 
full_model <- glmmTMB(count_sliccd ~ age_category + maternal_simd_quintile +healthboard_of_residence + rcs(time_period,knots2) +offset(log(count_denom)),
                      family = compois, data = final_model_data)

## summary of the model 
summary (full_model)

### table of model summary 
fullmodel_summ<- summary(full_model)
round(fullmodel_summ$coefficients$cond, 2)
full_model_tbl <- data.frame(fullmodel_summ$coefficients$cond)
setDT(full_model_tbl, keep.rownames = TRUE)
full_model_tbl$Estimate<- round(full_model_tbl$Estimate, 2)
full_model_tbl$P_value<- round(full_model_tbl$Pr...z.., 3)

###remove unneeded columns 
full_model_tbl$Std..Error <- NULL
full_model_tbl$z.value <- NULL

### pRR
full_model_tbl <- full_model_tbl%>%
  mutate(pRR = exp(Estimate))
round(confint(full_model), 4)
CI_model_with_time_summary <- data.frame(confint(full_model))
CI_model_with_time_summary <- CI_model_with_time_summary%>%
  mutate(pRR = exp(Estimate), CI_lower = exp(X2.5..), CI_upper = exp(X97.5..))
CI_model_with_time_summary$pRR <- round(CI_model_with_time_summary$pRR, 2)
CI_model_with_time_summary$CI_lower <- round(CI_model_with_time_summary$CI_lower, 2)
CI_model_with_time_summary$CI_upper <- round(CI_model_with_time_summary$CI_upper, 2)
CI_model_with_time_summary <- CI_model_with_time_summary%>%
  mutate(estimate_CI = paste0(Estimate," (", X2.5.., ", ", X97.5..,")"))%>%
  mutate(pRR_CI = paste0(pRR," (", CI_lower, ", ", CI_upper,")"))

###add in CI column to summary table
full_model_tbl <- full_model_tbl%>%
  mutate(Estimate_CI = CI_model_with_time_summary$estimate_CI)%>%
  mutate(pRR_CI = CI_model_with_time_summary$pRR_CI)



## ----------------------------------------------------------------------------
## --------- Live birth prevalence of DS --------------------------------------

### 1. reading in denominator data and prepping the sliccd numerator data #######

##### sliccd dataset prep -------------------------------------------------------
library(data.table)

##### dataset of live births only  ----------------------------------------------
sliccd_ds_extract <- sliccd_ds_extract%>%
  filter(pregnancy_end_type == "Livebirth")

## age categories in SLICCD data
sliccd_ds_extract <- sliccd_ds_extract %>%
  mutate(age_category = cut(maternal_age_at_conception,
                            breaks = c(-Inf, 20, 25, 30, 35, 40, Inf),
                            labels = c('<20', '20 - 24', '25-29', '30-34', '35-39', '40+'),
                            include.lowest = TRUE,
                            right = FALSE),
         age_category = ifelse(is.na(age_category), 'Unknown', as.character(age_category)))

sliccd_ds_extract <- sliccd_ds_extract%>%
  mutate(age_category = recode(age_category,
                               "<20" = "<25",
                               "20 - 24" = "<25"))

sliccd_ds_extract$age_category <- as.factor(sliccd_ds_extract$age_category)

## infant sex in sliccd 
# code the sex 1,2,3 into words 
sliccd_ds_extract <- sliccd_ds_extract %>%
  mutate(sex_2 = case_when(sex == 1 ~ 'Male',
                           sex == 2 ~ 'Female',
                           sex == 3 ~ 'Unknown',
                           TRUE ~ NA_character_))

sliccd_ds_extract$sex_2 <- as.factor(sliccd_ds_extract$sex_2)

#### sorting out structure and values in the dataset - characters/factors and 'unknown' groups 
str(sliccd_ds_extract$pregnancy_end_type)
sliccd_ds_extract$pregnancy_end_type <- as.factor(sliccd_ds_extract$pregnancy_end_type)

#### numerator data 
sliccd_numerator_data <- sliccd_ds_extract%>%
  group_by(time_period, age_category, maternal_simd_quintile, healthboard_of_residence, sex_2)%>%
  summarise(Count = n())

#### SIMD
sliccd_numerator_data$maternal_simd_quintile <- as.character(sliccd_numerator_data$maternal_simd_quintile)
sliccd_numerator_data <- sliccd_numerator_data%>%
  mutate( maternal_simd_quintile = replace_na(maternal_simd_quintile, "Unknown"))
sliccd_numerator_data$maternal_simd_quintile <- as.factor(sliccd_numerator_data$maternal_simd_quintile)

#### healthboard of residence 
sliccd_numerator_data <- sliccd_numerator_data%>%
  mutate(healthboard_of_residence = replace_na(healthboard_of_residence, "Unknown"))
sliccd_numerator_data$healthboard_of_residence<- as.factor(sliccd_numerator_data$healthboard_of_residence)
sliccd_numerator_data

##### read in denominator data for total births ---------------------------------
NRS_live_births <-read.csv("/PHI_conf/NIPT_evaluation_anon/Data/cohort_2/NRS_denom_data_livebirths.csv")

### recode SIMD
NRS_live_births <- NRS_live_births%>%
  mutate(maternal_simd_quintile = recode(maternal_simd_quintile,
                                         "1 (most deprived)" = "1",
                                         "5 (least deprived)" = "5"))

### re-code male/female
NRS_live_births <- NRS_live_births%>%
  mutate(sex_2 = recode(sex_2,
                        "M" = "Male",
                        "F" = "Female"))

### change to variables factors 
NRS_live_births$age_category <- as.factor(NRS_live_births$age_category)
NRS_live_births$healthboard_of_residence <- as.factor(NRS_live_births$healthboard_of_residence)
NRS_live_births$maternal_simd_quintile <- as.factor(NRS_live_births$maternal_simd_quintile)
NRS_live_births$sex_2 <- as.factor(NRS_live_births$sex_2)

###### merge the datasets
merged_data_live_births <- NRS_live_births%>%
  full_join(sliccd_numerator_data, by = c("time_period", "age_category", "maternal_simd_quintile", "healthboard_of_residence", "sex_2"))

#### rename the count columns 
merged_data_live_births <- merged_data_live_births%>%
  rename(count_denom = Count.x,
         count_sliccd = Count.y
  )

##### live births in both sliccd and denom data 
total_counts <- merged_data_live_births%>%
  summarise(total_denom = sum(count_denom, na.rm = TRUE),
            total_sliccd = sum(count_sliccd, na.rm = TRUE))

### 2. live birth prevalence ----------------------------------------------------
(total_counts$total_sliccd/ total_counts$total_denom * 10000 )

##### live birth prev per year in dataset --------------------------------------

### dataset of total counts per year 
total_counts_per_year <- merged_data_live_births%>%
  group_by(time_period)%>%
  summarise(total_denom = sum(count_denom, na.rm = TRUE),
            total_sliccd = sum(count_sliccd, na.rm = TRUE))

# function to calculate the prevalence and 95 CI in data using epitools
calculate_CI<- function(total_sliccd , total_denom){
  prevalence <- (total_sliccd / total_denom)* 10000
  ci <- pois.exact(total_sliccd, pt =  total_denom, conf.level = 0.95) * 10000
  ci_lower <- (ci$lower / total_denom)* 10000
  ci_upper <- (ci$upper / total_denom)* 10000
  list(prevalence = prevalence, ci_lower = ci$lower, ci_upper = ci$upper)
}

## dataset of LB prevalence per
prevalence_data_per_year <- total_counts_per_year%>%
  rowwise()%>%
  mutate(ci = list(calculate_CI(total_sliccd, total_denom)),
         prevalence = ci$prevalence,
         ci_lower = ci$ci_lower,
         ci_upper = ci$ci_upper,
         births_prev_ci = paste0(round(prevalence,2), " (", round(ci_lower, 2), "," ,round(ci_upper, 2), ")")) %>%
  ungroup()

total_counts%>%
  rowwise()%>%
  mutate(ci = list(calculate_CI(total_sliccd, total_denom)),
         prevalence = ci$prevalence,
         ci_lower = ci$ci_lower,
         ci_upper = ci$ci_upper,
         births_prev_ci = paste0(round(prevalence,2), " (", round(ci_lower, 2), "," ,round(ci_upper, 2), ")")) %>%
  ungroup()

### Modelling LB prevalence per year -------------------------------------------

#transforming year to 0-21
prevalence_data_per_year <- transform(prevalence_data_per_year, time_period = time_period - min(time_period))

## live birth prevalence is calculated using live births in sliccd and total live births in NRS to have live births per 10,000 live births 

#### a. assuming linear time trend ---------------------------------------------

livebirth_prev_linear1 <- glmmTMB(total_sliccd ~ time_period +offset(log(total_denom)),
                                  family = compois, data = prevalence_data_per_year)##total pregnancies / total NRS births 

summary(livebirth_prev_linear1)
round(confint(livebirth_prev_linear1), 2)
round(exp(confint(livebirth_prev_linear1)), 2)

#### plotting linear model ----------------------------------------------------

## For live birth prevalence 

# Predict values
livebirths_COMpoissonm1 <- data.frame(time_period = prevalence_data_per_year$time_period, total_denom = prevalence_data_per_year$total_denom, total_sliccd = prevalence_data_per_year$total_sliccd)
LBpredicted_values <- predict(livebirth_prev_linear1, newdata = livebirths_COMpoissonm1, type = "link", se.fit = TRUE)

# Extract predicted values and standard errors
LBpredicted_means <- LBpredicted_values$fit
LBstandard_errors <- LBpredicted_values$se.fit

# Calculate confidence intervals
z_value <- qnorm(0.975)  # 1.96 for a 95% confidence interval
LBlower_bound <- exp(LBpredicted_means - z_value * LBstandard_errors)
LBupper_bound <- exp(LBpredicted_means + z_value * LBstandard_errors)

# Display results
LBresult_df_COM1 <- data.frame(
  time_period = livebirths_COMpoissonm1$time_period,
  Predicted_Mean = (exp(LBpredicted_means) / prevalence_data_per_year$total_denom) * 10000,
  Lower_Bound = (LBlower_bound / prevalence_data_per_year$total_denom) * 10000,
  Upper_Bound = (LBupper_bound / prevalence_data_per_year$total_denom) * 10000,
  livebirth_prev = (prevalence_data_per_year$total_sliccd / prevalence_data_per_year$total_denom) * 10000  # Include prevs in the result_df
)
print(LBresult_df_COM1)

# Plotting linear model (Figure 5)
prevalence_data_per_year$livebirth_prev <- (prevalence_data_per_year$total_sliccd / prevalence_data_per_year$total_denom) * 10000

# Scatterplot of original data

scatterplot <- ggplot(prevalence_data_per_year, aes(x = time_period+2000, y = livebirth_prev)) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.6, colour = "steel blue")+
  geom_line(colour = "cornflower blue")+  
  geom_point( shape = 21, size = 3, fill = "steel blue")+   
  labs( x = "Year", y = "Live birth prevalence, per 10,000 live births")

scatterplot

# Add fitted values and confidence intervals as blue lines
lineplot <- scatterplot +
  geom_line(data = LBresult_df_COM1, aes(x = time_period+2000, y = Predicted_Mean), color = "blue") +
  geom_ribbon(data = LBresult_df_COM1, aes(x = time_period+2000, ymin = Lower_Bound, ymax = Upper_Bound), fill = "blue", alpha = 0.2) +
  
  # Customize plot labels and appearance
  labs(x = "Year",
       y = "Live birth prevalence per 10,000 live births") +
  scale_x_continuous(
    breaks = seq(2000,2021, by = 1),
    labels = seq (2000, 2021, by = 1))+
  theme_minimal()+
  ylim(0, 30)+
  theme(text = element_text(size=20))

#### b. assuming non-linear time trend -----------------------------------------

## 3 knots (RCS)

livebirth_prev_nonlin1 <- glmmTMB(total_sliccd ~ rcs(time_period,3) +offset(log(total_denom)),
                                  family = compois, data = prevalence_data_per_year)##total pregnancies / total NRS births 

summary(livebirth_prev_nonlin1)
round(confint(livebirth_prev_nonlin1), 2)
round(exp(confint(livebirth_prev_nonlin1)), 2)

### 4 knots 

livebirth_prev_nonlin1_3 <- glmmTMB(total_sliccd ~ rcs(time_period,4) +offset(log(total_denom)),
                                    family = compois, data = prevalence_data_per_year)##total pregnancies / total NRS births 

summary(livebirth_prev_nonlin1_3)
round(confint(livebirth_prev_nonlin1_3), 2)
round(exp(confint(livebirth_prev_nonlin1_3)), 2)

###3. Livebirth prevalence of DS and maternal age  -----------------------------
### dataset with maternal age group and LB prevalence per age group 

LB_MA_group_data <- merged_data_live_births%>%
  group_by(age_category)%>%
  summarise(count_sliccd = sum(count_sliccd, na.rm = TRUE),
            count_denom = sum(count_denom, na.rm = TRUE))%>%
  ungroup()

###removing unknown category 

LB_MA_group_data <- LB_MA_group_data%>%
  filter(!age_category == 'Unknown')


##### ADD live BIRTH PREVALENCE COLUMN TO DATASET 

LB_MA_group_data <- LB_MA_group_data%>%
  mutate(livebirth_prevalence = count_sliccd/count_denom * 10000)

### 95% CIs 

LB_MA_group_data <- LB_MA_group_data%>%
  rowwise()%>%
  mutate(ci = list(calculate_CI(count_sliccd, count_denom)),
         prevalence = ci$prevalence,
         ci_lower = ci$ci_lower,
         ci_upper = ci$ci_upper,
         births_prev_ci = paste0(round(prevalence,2), " (", round(ci_lower, 2), "," ,round(ci_upper, 2), ")")) %>%
  ungroup()

##### model of age and prevalence model - VGLM #########################
LB_vglm1<-vglm(LB_MA_group_data$count_sliccd ~ LB_MA_group_data$age_category, poissonff(bred=TRUE), offset=log(count_denom), 
            data=LB_MA_group_data)

LB_vglm1summary <- summary(LB_vglm1)
show.summary.vglm(LB_vglm1summary)
LB_vglm1summary <- summaryvglm(LB_vglm1)
CI <- data.frame(confint(LB_vglm1))

###create summary table
MA_LB_vglm1summary_tbl <- data.frame (coef(LB_vglm1))
MA_LB_vglm1summary_tbl <- MA_LB_vglm1summary_tbl%>%
  mutate(Lower = CI$X2.5.., Upper = CI$X97.5.. )

## rename table columns 
colnames(MA_LB_vglm1summary_tbl)[1] <- "Coefficients_Est"

## add in coefficients names column 

MA_LB_vglm1summary_tbl <- MA_LB_vglm1summary_tbl%>%
  mutate(Coefficients = c("Intercept" ,"25-29", "30-34", "35-39", "40+"))
#reorder columns in data
MA_LB_vglm1summary_tbl <- MA_LB_vglm1summary_tbl[c("Coefficients", "Coefficients_Est", "Lower", "Upper")]

##add in pRR for each coefficient and CI
MA_LB_vglm1summary_tbl <- MA_LB_vglm1summary_tbl%>%
  mutate(pRR = exp(coef(LB_vglm1)), pRR_CI_lower = exp(Lower), PRR_CI_upper = exp(Upper))
### add in p values

MA_LB_vglm1summary_tbl<- MA_LB_vglm1summary_tbl%>%
  mutate("P value" = coef(summary(LB_vglm1))[,'Pr(>|z|)'])

MA_LB_vglm1summary_tbl$pRR_CI_lower <- round(MA_LB_vglm1summary_tbl$pRR_CI_lower, 2)
MA_LB_vglm1summary_tbl$pRR_CI_upper <- round(MA_LB_vglm1summary_tbl$PRR_CI_upper, 2)

MA_LB_vglm1summary_tbl$pRR <-round(MA_LB_vglm1summary_tbl$pRR,2)

MA_LB_vglm1summary_tbl$`P value` <-round(MA_LB_vglm1summary_tbl$`P value`,3)
MA_LB_vglm1summary_tbl$Coefficients_Est <-round(MA_LB_vglm1summary_tbl$Coefficients_Est,2)
MA_LB_vglm1summary_tbl$Lower <-round(MA_LB_vglm1summary_tbl$Lower,2)
MA_LB_vglm1summary_tbl$Upper <-round(MA_LB_vglm1summary_tbl$Upper,2)

MA_LB_vglm1summary_tbl <- MA_LB_vglm1summary_tbl%>%
  mutate(estimate_CI = paste0(Coefficients_Est," (", Lower, ", ", Upper,")"))%>%
  mutate(pRR_CI = paste0(pRR," (", MA_LB_vglm1summary_tbl$pRR_CI_lower, ", ", pRR_CI_upper,")"))


##########plotting MA predicted prevalence (Figure 3d)####################

## For birth prevalence 

# Predict values
age_vglm1 <- data.frame(age_category = LB_MA_group_data$age_category, count_sliccd = LB_MA_group_data$count_sliccd, count_denom = LB_MA_group_data$count_denom)
predicted_values <- predict(LB_vglm1, newdata = age_vglm1, type = "link", se.fit = TRUE)

# Extract predicted values and standard errors
predicted_means <- predicted_values$fit
standard_errors <- predicted_values$se.fit

# Calculate confidence intervals
z_value <- qnorm(0.975)  # 1.96 for a 95% confidence interval
lower_bound <- exp(predicted_means - z_value * standard_errors)
upper_bound <- exp(predicted_means + z_value * standard_errors)

# Display results
LB_MA_result_df_vglm1 <- data.frame(
  age_category = age_vglm1$age_category,
  Predicted_Mean = (exp(predicted_means) / LB_MA_group_data$count_denom) * 10000,
  Lower_Bound = (lower_bound / LB_MA_group_data$count_denom) * 10000,
  Upper_Bound = (upper_bound / LB_MA_group_data$count_denom) * 10000,
  birth_prev = (LB_MA_group_data$count_sliccd / LB_MA_group_data$count_denom) * 10000  # Include prevs in the result_df
)
print(LB_MA_result_df_vglm1)


######## plot this data --------------------------------------------------------

ggplot(LB_MA_result_df_vglm1, aes(x = age_category, y = birth_prev)) +
  geom_line(aes(y = Predicted_Mean, group = 1), colour = "purple")+
  geom_errorbar(aes(ymin = Lower_Bound, ymax = Upper_Bound), colour = "purple", alpha = 0.5, width = 0.4)+
  geom_point(aes(y = birth_prev), shape = 19, size = 2, alpha = 0.7, color = "black")+
  labs(x = "Maternal age group (yrs)", y = "Live birth prevalence per 10,000 live births")+
  theme_minimal()+
  theme(text = element_text(size=20))

## Final plot

MA_plot3 <- ggplot(LB_MA_result_df_vglm1, aes(x=age_category, y = birth_prev))+
  geom_line(aes(y = birth_prev, group = 1), colour = "purple")+
  geom_errorbar(aes(ymin = Lower_Bound, ymax = Upper_Bound), width = 0.6, colour = "purple")+
  geom_point(aes(y = birth_prev), shape = 19, size = 3, color = "purple")+
  labs(x = "Maternal age group (yrs)", y = "Live birth prevalence per 10,000 live births", title = "d")+
  ylim(0,120)+
  theme_minimal()+
  theme(
    axis.title.x = element_text(margin = margin(t = 15)),  # Add space above x-axis title
    axis.title.y = element_text(margin = margin(r = 15))   # Add space to the right of y-axis title
  )+
  theme(text = element_text(size=20))

MA_plot3

##### Maternal age and time ----------------------------------------------------

######## creating time and MA dataset to use ------------------------------

MA_group_data_bytime <- merged_data_live_births%>%
  group_by(time_period, age_category)%>%
  summarise(count_sliccd = sum(count_sliccd, na.rm = TRUE),
            count_denom = sum(count_denom, na.rm = TRUE))%>%
  ungroup()

### removing unknown category 

MA_group_data_bytime <- MA_group_data_bytime%>%
  filter(!age_category == 'Unknown')

##### ADD TOTAL BIRTH PREVALENCE COLUMN TO DATASET 

MA_group_data_bytime <- MA_group_data_bytime%>%
  mutate(totalbirth_prevalence = count_sliccd/count_denom * 10000)

#### change time to 0-21:

MA_group_data_bytime <- transform(MA_group_data_bytime, time_period = time_period - min(time_period))


#### maternal age adjusted by linear time trend model

model_with_time <- glmmTMB(count_sliccd ~ age_category + time_period +offset(log(count_denom)),
                           family = compois, data = MA_group_data_bytime) 

summary (model_with_time)
model_time_summ<- summary(model_with_time)

MA_model_with_time <- data.frame(model_time_summ$coefficients$cond)

setDT(MA_model_with_time, keep.rownames = TRUE)

MA_model_with_time$Estimate<- round(MA_model_with_time$Estimate, 2)

MA_model_with_time$Pr...z..<- round(MA_model_with_time$Pr...z.., 3)

###remove unneeded columns 

MA_model_with_time$Std..Error <- NULL
MA_model_with_time$z.value <- NULL

CI_model_with_time_summary <- data.frame(round(confint(model_with_time), 2))

CI_model_with_time_summary <- CI_model_with_time_summary%>%
  mutate(pRR = exp(Estimate), CI_lower = exp(X2.5..), CI_upper = exp(X97.5..))

CI_model_with_time_summary$pRR <- round(CI_model_with_time_summary$pRR, 2)
CI_model_with_time_summary$CI_lower <- round(CI_model_with_time_summary$CI_lower, 2)
CI_model_with_time_summary$CI_upper <- round(CI_model_with_time_summary$CI_upper, 2)

CI_model_with_time_summary <- CI_model_with_time_summary%>%
  mutate(estimate_CI = paste0(Estimate," (", X2.5.., ", ", X97.5..,")"))%>%
  mutate(pRR_CI = paste0(pRR," (", CI_lower, ", ", CI_upper,")"))

###add in CI column to summary table 

LB_MA_model_with_time<- MA_model_with_time%>%
  mutate(Estimate_CI = CI_model_with_time_summary$estimate_CI)%>%
  mutate(pRR_CI = CI_model_with_time_summary$pRR_CI)


####### plotting maternal age and time ---------------------------------------------

# Predict values
MA_time <- data.frame(time_period = MA_group_data_bytime$time_period, age_category = MA_group_data_bytime$age_category, count_sliccd = MA_group_data_bytime$count_sliccd, count_denom = MA_group_data_bytime$count_denom)
predicted_values <- predict(model_with_time, newdata = MA_time, type = "link", se.fit = TRUE)

# Extract predicted values and standard errors
predicted_means <- predicted_values$fit
standard_errors <- predicted_values$se.fit

# Calculate confidence intervals
z_value <- qnorm(0.975)  # 1.96 for a 95% confidence interval
lower_bound <- exp(predicted_means - z_value * standard_errors)
upper_bound <- exp(predicted_means + z_value * standard_errors)

# Display results
MA_LB_TIME_result_df_vglm1 <- data.frame(
  time_period = MA_time$time_period,
  age_category = MA_time$age_category,
  Predicted_Mean = (exp(predicted_means) / MA_group_data_bytime$count_denom) * 10000,
  Lower_Bound = (lower_bound / MA_group_data_bytime$count_denom) * 10000,
  Upper_Bound = (upper_bound / MA_group_data_bytime$count_denom) * 10000,
  birth_prev = (MA_group_data_bytime$count_sliccd / MA_group_data_bytime$count_denom) * 10000  # Include prevs in the result_df
)

MA_group_data_bytime <- MA_group_data_bytime%>%
  mutate(predicted_mean = MA_LB_TIME_result_df_vglm1$Predicted_Mean, lower_bound = MA_LB_TIME_result_df_vglm1$Lower_Bound, upper_bound = MA_LB_TIME_result_df_vglm1$Upper_Bound)

### plotting predicted prevalence per year + MA

## Final plot for manuscript (Figure 4d)

ma_plot4 <- ggplot(MA_group_data_bytime, aes(x = time_period+2000, y = totalbirth_prevalence, group = age_category))+
  geom_line(aes(colour = age_category, y = predicted_mean))+
  geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound, fill = age_category), alpha = 0.3)+
  geom_line(aes(colour = age_category), alpha = 0.3, 
            size=0.5)+
  geom_point(aes(colour = age_category), size=3)+
  ylim(0,150)+
  theme_minimal()+ 
  labs(x = "Year",
       y = "Live birth prevalence per 10,000 live births",
       title = "c",
       colour = "Maternal age group",
       fill = "Maternal age group")+
  theme(
    axis.title.x = element_text(margin = margin(t = 15)),  # Add space above x-axis title
    axis.title.y = element_text(margin = margin(r = 15))   # Add space to the right of y-axis title
  )+ 
  theme(text = element_text(size=20))

ma_plot4

### 4. SIMD and live birth prevalence ------------------------------------------

### dataset with SIMD and prevalence per age group 

SIMD_group_data <- merged_data_live_births%>%
  group_by(maternal_simd_quintile)%>%
  summarise(count_sliccd = sum(count_sliccd, na.rm = TRUE),
            count_denom = sum(count_denom, na.rm = TRUE))%>%
  ungroup()

###removing unknown category 

SIMD_group_data <- SIMD_group_data%>%
  filter(!maternal_simd_quintile == 'Unknown')

##### ADD live BIRTH PREVALENCE COLUMN TO DATASET 

SIMD_group_data <- SIMD_group_data%>%
  mutate(livebirth_prevalence = count_sliccd/count_denom * 10000)

SIMD_group_data <- SIMD_group_data%>%
  rowwise()%>%
  mutate(ci = list(calculate_CI(count_sliccd, count_denom)),
         prevalence = ci$prevalence,
         ci_lower = ci$ci_lower,
         ci_upper = ci$ci_upper,
         births_prev_ci = paste0(round(prevalence,2), " (", round(ci_lower, 2), "," ,round(ci_upper, 2), ")")) %>%
  ungroup()


##### Model of SIMD and LB prevalence --------------------------------

vglm1<-vglm(SIMD_group_data$count_sliccd ~ SIMD_group_data$maternal_simd_quintile, poissonff(bred=TRUE), offset=log(count_denom), 
            data=SIMD_group_data)
vglm1summary <- summary(vglm1)
show.summary.vglm(vglm1summary)
vglm1summary <- summaryvglm(vglm1)
CI <- data.frame(confint(vglm1))

###create summary table
vglm1summary_tbl <- data.frame (coef(vglm1))

vglm1summary_tbl <- vglm1summary_tbl%>%
  mutate(Lower = CI$X2.5.., Upper = CI$X97.5.. )
## rename table columns 
colnames(vglm1summary_tbl)[1] <- "Coefficients_Est"

## add in coefficients names column 
vglm1summary_tbl <- vglm1summary_tbl%>%
  mutate(Coefficients = c("Intercept" ,"2", "3", "4", "5"))

#reorder columns in data
vglm1summary_tbl <- vglm1summary_tbl[c("Coefficients", "Coefficients_Est", "Lower", "Upper")]

##add in pRR for each coefficient and CI
vglm1summary_tbl <- vglm1summary_tbl%>%
  mutate(pRR = exp(coef(vglm1)), pRR_CI_lower = exp(Lower), PRR_CI_upper = exp(Upper))
### add in p values

vglm1summary_tbl<- vglm1summary_tbl%>%
  mutate("P value" = coef(summary(vglm1))[,'Pr(>|z|)'])

vglm1summary_tbl$pRR_CI_lower <- round(vglm1summary_tbl$pRR_CI_lower, 2)
vglm1summary_tbl$pRR_CI_upper <- round(vglm1summary_tbl$PRR_CI_upper, 2)
vglm1summary_tbl$pRR <-round(vglm1summary_tbl$pRR,2)
vglm1summary_tbl$`P value` <-round(vglm1summary_tbl$`P value`,3)
vglm1summary_tbl$Coefficients_Est <-round(vglm1summary_tbl$Coefficients_Est,2)
vglm1summary_tbl$Lower <-round(vglm1summary_tbl$Lower,2)
vglm1summary_tbl$Upper <-round(vglm1summary_tbl$Upper,2)

LB_SIMD_vglm1summary_tbl <- vglm1summary_tbl%>%
  mutate(estimate_CI = paste0(Coefficients_Est," (", Lower, ", ", Upper,")"))%>%
  mutate(pRR_CI = paste0(pRR," (", vglm1summary_tbl$pRR_CI_lower, ", ", pRR_CI_upper,")"))


##### plotting SIMD predicted LB prevalence (Figure 3e) -------------------------------------

## For birth prevalence 

# Predict values
LBSIMD_vglm1 <- data.frame(maternal_simd_quintile = SIMD_group_data$maternal_simd_quintile, count_sliccd = SIMD_group_data$count_sliccd, count_denom = SIMD_group_data$count_denom)
predicted_values <- predict(vglm1, newdata = LBSIMD_vglm1, type = "link", se.fit = TRUE)

# Extract predicted values and standard errors
predicted_means <- predicted_values$fit
standard_errors <- predicted_values$se.fit

# Calculate confidence intervals
z_value <- qnorm(0.975)  # 1.96 for a 95% confidence interval
lower_bound <- exp(predicted_means - z_value * standard_errors)
upper_bound <- exp(predicted_means + z_value * standard_errors)

# Display results
LB_SIMD_result_df_vglm1 <- data.frame(
  maternal_simd_quintile = LBSIMD_vglm1$maternal_simd_quintile,
  Predicted_Mean = (exp(predicted_means) / SIMD_group_data$count_denom) * 10000,
  Lower_Bound = (lower_bound / SIMD_group_data$count_denom) * 10000,
  Upper_Bound = (upper_bound / SIMD_group_data$count_denom) * 10000,
  birth_prev = (SIMD_group_data$count_sliccd / SIMD_group_data$count_denom) * 10000  # Include prevs in the result_df
)
print(LB_SIMD_result_df_vglm1)

#### 


########## plot this data --------------------------------------------------------

ggplot(LB_SIMD_result_df_vglm1, aes(x = maternal_simd_quintile, y = birth_prev)) +
  geom_line(aes(y = Predicted_Mean, group = 1), colour = "purple")+
  geom_errorbar(aes(ymin = Lower_Bound, ymax = Upper_Bound), colour = "purple", alpha = 0.5, width = 0.4)+
  geom_point(aes(y = birth_prev), shape = 19, size = 2, color = "black")+ labs(
    x = "Maternal SIMD group (1 - most deprived, 5 - least deprived)", y = "Live birth prevalence (per 10,000 live births)")+
  theme_minimal()+
  theme(text = element_text(size=20))+
  ylim(0,20)


SIMD_plot3 <- ggplot(LB_SIMD_result_df_vglm1, aes(x=maternal_simd_quintile, y = birth_prev))+
  geom_line(group =1, colour = "purple")+
  geom_errorbar(aes(ymin = Lower_Bound, ymax = Upper_Bound), width = 0.6, colour = "purple")+
  geom_point(aes(y = birth_prev), shape = 19, size = 3, color = "purple")+
  labs(x = "Maternal SIMD group (1 - most deprived, 5 - least deprived)", y = "Live birth prevalence per 10,000 live births", title = "e")+
  ylim(0,30)+
  theme_minimal()+
  theme(
    axis.title.x = element_text(margin = margin(t = 15)),  # Add space above x-axis title
    axis.title.y = element_text(margin = margin(r = 15))   # Add space to the right of y-axis title
  )+
  theme(text = element_text(size=20))

SIMD_plot3

##### SIMD and time ------------------------------------------------------

### creating time and SIMD dataset

SIMD_group_data_bytime <- merged_data_live_births%>%
  group_by(time_period, maternal_simd_quintile)%>%
  summarise(count_sliccd = sum(count_sliccd, na.rm = TRUE),
            count_denom = sum(count_denom, na.rm = TRUE))%>%
  ungroup()

###removing unknown category 

SIMD_group_data_bytime <- SIMD_group_data_bytime%>%
  filter(!maternal_simd_quintile == 'Unknown')

##### ADD TOTAL BIRTH PREVALENCE COLUMN TO DATASET 

SIMD_group_data_bytime <- SIMD_group_data_bytime%>%
  mutate(totalbirth_prevalence = count_sliccd/count_denom * 10000)

#### change time to 0-21:

SIMD_group_data_bytime <- transform(SIMD_group_data_bytime, time_period = time_period - min(time_period))

#### SIMD adjusted by linear time trend 

model_with_time <- glmmTMB(count_sliccd ~ maternal_simd_quintile + time_period +offset(log(count_denom)),
                           family = compois, data = SIMD_group_data_bytime) 

summary (model_with_time)
model_time_summ<- summary(model_with_time)

SIMD_model_with_time <- data.frame(model_time_summ$coefficients$cond)

setDT(SIMD_model_with_time, keep.rownames = TRUE)

SIMD_model_with_time$Estimate<- round(SIMD_model_with_time$Estimate, 2)

SIMD_model_with_time$Pr...z..<- round(SIMD_model_with_time$Pr...z.., 3)

###remove unneeded columns 

SIMD_model_with_time$Std..Error <- NULL
SIMD_model_with_time$z.value <- NULL

CI_model_with_time_summary <- data.frame(round(confint(model_with_time), 2))

CI_model_with_time_summary <- CI_model_with_time_summary%>%
  mutate(pRR = exp(Estimate), CI_lower = exp(X2.5..), CI_upper = exp(X97.5..))

CI_model_with_time_summary$pRR <- round(CI_model_with_time_summary$pRR, 2)
CI_model_with_time_summary$CI_lower <- round(CI_model_with_time_summary$CI_lower, 2)
CI_model_with_time_summary$CI_upper <- round(CI_model_with_time_summary$CI_upper, 2)

CI_model_with_time_summary <- CI_model_with_time_summary%>%
  mutate(estimate_CI = paste0(Estimate," (", X2.5.., ", ", X97.5..,")"))%>%
  mutate(pRR_CI = paste0(pRR," (", CI_lower, ", ", CI_upper,")"))

###add in CI column to summary table 

LB_SIMD_model_with_time<- SIMD_model_with_time%>%
  mutate(Estimate_CI = CI_model_with_time_summary$estimate_CI)%>%
  mutate(pRR_CI = CI_model_with_time_summary$pRR_CI)

####### plotting SIMD with time ---------------------------------------------------

# Predict values
SIMD_time <- data.frame(time_period = SIMD_group_data_bytime$time_period, maternal_simd_quintile = SIMD_group_data_bytime$maternal_simd_quintile, count_sliccd = SIMD_group_data_bytime$count_sliccd, count_denom = SIMD_group_data_bytime$count_denom)
predicted_values <- predict(model_with_time, newdata = SIMD_time, type = "link", se.fit = TRUE)

# Extract predicted values and standard errors
predicted_means <- predicted_values$fit
standard_errors <- predicted_values$se.fit

# Calculate confidence intervals
z_value <- qnorm(0.975)  # 1.96 for a 95% confidence interval
lower_bound <- exp(predicted_means - z_value * standard_errors)
upper_bound <- exp(predicted_means + z_value * standard_errors)

# Display results
LB_SIMD_result_df_vglm1 <- data.frame(
  time_period = SIMD_time$time_period,
  maternal_simd_quintile = SIMD_time$maternal_simd_quintile,
  Predicted_Mean = (exp(predicted_means) / SIMD_group_data_bytime$count_denom) * 10000,
  Lower_Bound = (lower_bound / SIMD_group_data_bytime$count_denom) * 10000,
  Upper_Bound = (upper_bound / SIMD_group_data_bytime$count_denom) * 10000,
  birth_prev = (SIMD_group_data_bytime$count_sliccd / SIMD_group_data_bytime$count_denom) * 10000  # Include prevs in the result_df
)
print(LB_SIMD_result_df_vglm1)

LB_SIMD_group_data_bytime <- SIMD_group_data_bytime%>%
  mutate(predicted_mean = LB_SIMD_result_df_vglm1$Predicted_Mean, lower_bound = LB_SIMD_result_df_vglm1$Lower_Bound, upper_bound = LB_SIMD_result_df_vglm1$Upper_Bound)

#### plotting predicted prevalence per year + SIMD (Figure 4d)------------------------------------

ggplot(LB_SIMD_group_data_bytime, aes(x = time_period+2000, y = totalbirth_prevalence, group = maternal_simd_quintile))+
  geom_line(aes(colour = maternal_simd_quintile), alpha = 0.3, 
            size=0.5)+
  geom_line(aes(colour = maternal_simd_quintile, y = predicted_mean))+
  geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound, fill = maternal_simd_quintile), alpha = 0.3)+
  geom_point(aes(colour = maternal_simd_quintile), size=3)+
  labs(x = "Year of birth",
       y = "Live birth prevalence per 10,000 live births",
       colour = "Maternal SIMD group",
       fill = "Maternal SIMD group")+
  facet_wrap(maternal_simd_quintile~., ncol=1)+
  ylim(0,20)+
  theme_minimal()+
  theme(text = element_text(size=20))


simd_plot4 <- ggplot(LB_SIMD_group_data_bytime, aes(x = time_period+2000, y = totalbirth_prevalence, group = maternal_simd_quintile))+
  geom_line(aes(colour = maternal_simd_quintile), alpha = 0.3, 
            size=0.5)+
  geom_line(aes(colour = maternal_simd_quintile, y = predicted_mean))+
  geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound, fill = maternal_simd_quintile), alpha = 0.3)+
  geom_point(aes(colour = maternal_simd_quintile), size=3)+
  ylim(0,40)+
  labs(x = "Year",
       y = "Live birth prevalence per 10,000 live births",
       title = "d",
       colour = "Maternal SIMD group",
       fill = "Maternal SIMD group")+
  theme_minimal()+
  theme(
    axis.title.x = element_text(margin = margin(t = 15)),  # Add space above x-axis title
    axis.title.y = element_text(margin = margin(r = 15))   # Add space to the right of y-axis title
  )+
  theme(text = element_text(size=20))

simd_plot4

## 5. HEALTHBOARD and livebirth prevalence of DS  ------------------------------

### dataset with HB and LB prevalence 

HB_group_data <- merged_data_live_births%>%
  group_by(healthboard_of_residence)%>%
  summarise(count_sliccd = sum(count_sliccd, na.rm = TRUE),
            count_denom = sum(count_denom, na.rm = TRUE))%>%
  ungroup()

###removing unknown category 

HB_group_data <- HB_group_data%>%
  filter(!healthboard_of_residence == 'Unknown')


##### ADD live BIRTH PREVALENCE COLUMN TO DATASET 

HB_group_data <- HB_group_data%>%
  mutate(livebirth_prevalence = count_sliccd/count_denom * 10000)

HB_group_data <- HB_group_data%>%
  rowwise()%>%
  mutate(ci = list(calculate_CI(count_sliccd, count_denom)),
         prevalence = ci$prevalence,
         ci_lower = ci$ci_lower,
         ci_upper = ci$ci_upper,
         births_prev_ci = paste0(round(prevalence,2), " (", round(ci_lower, 2), "," ,round(ci_upper, 2), ")")) %>%
  ungroup()

##### plot prevalence by HB

ggplot(HB_group_data, aes(x=healthboard_of_residence, y = livebirth_prevalence))+
  geom_bar(stat="identity", fill = "purple", alpha= 0.5)+
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, colour = "purple")+
  labs(title = "The live birth prevalence of babies with DS in Scotland by maternal healthboard of residence", x = "NHS healthboard of residence", y = "Live birth prevalence (per 10,000 live births)")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.4), text = element_text(size=15))

## grouping smaller health boards

HB_group_data<- HB_group_data%>%
  mutate(healthboard_of_residence = if_else(healthboard_of_residence %in% c('NHS Western Isles', 'NHS Orkney', 'NHS Shetland'),
                                            'NHS Island Boards',
                                            healthboard_of_residence))

HB_group_data<- HB_group_data%>%
  group_by(healthboard_of_residence)%>%
  summarise(count_sliccd = sum(count_sliccd), count_denom = sum(count_denom), .groups = 'drop')%>%
  mutate(totalbirth_prevalence = count_sliccd/ count_denom * 10000)

HB_group_data <- HB_group_data%>%
  rowwise()%>%
  mutate(ci = list(calculate_CI(count_sliccd, count_denom)),
         prevalence = ci$prevalence,
         ci_lower = ci$ci_lower,
         ci_upper = ci$ci_upper,
         births_prev_ci = paste0(round(prevalence,2), " (", round(ci_lower, 2), "," ,round(ci_upper, 2), ")")) %>%
  ungroup()

##### healthboard and LB prevalence model -----------------------------------

vglm1<-vglm(HB_group_data$count_sliccd ~ HB_group_data$healthboard_of_residence, poissonff(bred=TRUE), offset=log(count_denom), 
            data=HB_group_data)
vglm1summary <- summary(vglm1)
show.summary.vglm(vglm1summary)
vglm1summary <- summaryvglm(vglm1)
CI <- data.frame(confint(vglm1))

###create summary table

vglm1summary_tbl <- data.frame (coef(vglm1))
vglm1summary_tbl <- vglm1summary_tbl%>%
  mutate(Lower = CI$X2.5.., Upper = CI$X97.5.. )
## rename table columns 
colnames(vglm1summary_tbl)[1] <- "Coefficients_Est"

## add in coefficients names column 

vglm1summary_tbl <- vglm1summary_tbl%>%
  mutate(Coefficients = c("Intercept" ,     "NHS Borders"     ,              "NHS Dumfries and Galloway"    , "NHS Fife"    ,                 
                          "NHS Forth Valley"  ,            "NHS Grampian"      ,            "NHS Greater Glasgow and Clyde", "NHS Highland"  ,               
                          "NHS Island Boards"     ,        "NHS Lanarkshire"    ,           "NHS Lothian"    ,               "NHS Tayside"))
#reorder columns in data
vglm1summary_tbl <- vglm1summary_tbl[c("Coefficients", "Coefficients_Est", "Lower", "Upper")]

##add in pRR for each coefficient and CI
vglm1summary_tbl <- vglm1summary_tbl%>%
  mutate(pRR = exp(coef(vglm1)), pRR_CI_lower = exp(Lower), PRR_CI_upper = exp(Upper))
### add in p values

vglm1summary_tbl<- vglm1summary_tbl%>%
  mutate("P value" = coef(summary(vglm1))[,'Pr(>|z|)'])

vglm1summary_tbl$pRR_CI_lower <- round(vglm1summary_tbl$pRR_CI_lower, 2)
vglm1summary_tbl$pRR_CI_upper <- round(vglm1summary_tbl$PRR_CI_upper, 2)

vglm1summary_tbl$pRR <-round(vglm1summary_tbl$pRR,2)

vglm1summary_tbl$`P value` <-round(vglm1summary_tbl$`P value`,3)
vglm1summary_tbl$Coefficients_Est <-round(vglm1summary_tbl$Coefficients_Est,2)
vglm1summary_tbl$Lower <-round(vglm1summary_tbl$Lower,2)
vglm1summary_tbl$Upper <-round(vglm1summary_tbl$Upper,2)

HB_LB_vglm1summary_tbl <- vglm1summary_tbl%>%
  mutate(estimate_CI = paste0(Coefficients_Est," (", Lower, ", ", Upper,")"))%>%
  mutate(pRR_CI = paste0(pRR," (", vglm1summary_tbl$pRR_CI_lower, ", ", pRR_CI_upper,")"))

########## plotting MA predicted LB prevalence -------------------------------------

## For birth prevalence 

# Predict values
HB_vglm1 <- data.frame(healthboard_of_residence = HB_group_data$healthboard_of_residence, count_sliccd = HB_group_data$count_sliccd, count_denom = HB_group_data$count_denom)
predicted_values <- predict(vglm1, newdata = HB_vglm1, type = "link", se.fit = TRUE)

# Extract predicted values and standard errors
predicted_means <- predicted_values$fit
standard_errors <- predicted_values$se.fit

# Calculate confidence intervals
z_value <- qnorm(0.975)  # 1.96 for a 95% confidence interval
lower_bound <- exp(predicted_means - z_value * standard_errors)
upper_bound <- exp(predicted_means + z_value * standard_errors)

# Display results
HB_LB_result_df_vglm1 <- data.frame(
  healthboard_of_residence = HB_vglm1$healthboard_of_residence,
  Predicted_Mean = (exp(predicted_means) / HB_group_data$count_denom) * 10000,
  Lower_Bound = (lower_bound / HB_group_data$count_denom) * 10000,
  Upper_Bound = (upper_bound / HB_group_data$count_denom) * 10000,
  birth_prev = (HB_group_data$count_sliccd / HB_group_data$count_denom) * 10000  # Include prevs in the result_df
)
print(HB_LB_result_df_vglm1)

# final plot for manuscript (Figure 3e)

HB_PLOT2<- ggplot(HB_LB_result_df_vglm1, aes(x = healthboard_of_residence, y = birth_prev)) +
  geom_bar(stat = "identity", fill = "purple", alpha = 0.3)+
  coord_flip()+
  geom_errorbar(aes(ymin = Lower_Bound, ymax = Upper_Bound), colour = "purple", alpha = 0.7, width = 0.4)+
  geom_point(aes(y = birth_prev), shape = 19, size = 3, color = "purple")+
  labs(x = "Maternal NHS health board of residence", y = "Live birth prevalence per 10,000 live births", title = "f")+
  theme_minimal()+
  theme(
    axis.title.x = element_text(margin = margin(t = 15)),  # Add space above x-axis title
    axis.title.y = element_text(margin = margin(r = 15))   # Add space to the right of y-axis title
  )+
  theme(text = element_text(size=20))+
  ylim(0,30)

HB_PLOT2


#### healthboard and time---------------------------------------------

####### creating time and hb dataset to use -----------------------------

HB_group_data_bytime <- merged_data_live_births%>%
  group_by(time_period, healthboard_of_residence)%>%
  summarise(count_sliccd = sum(count_sliccd, na.rm = TRUE),
            count_denom = sum(count_denom, na.rm = TRUE))%>%
  ungroup()

###removing unknown category 

HB_group_data_bytime <- HB_group_data_bytime%>%
  filter(!healthboard_of_residence == 'Unknown')


HB_group_data_bytime<- HB_group_data_bytime%>%
  mutate(healthboard_of_residence = if_else(healthboard_of_residence %in% c('NHS Western Isles', 'NHS Orkney', 'NHS Shetland'),
                                            'NHS Island Boards',
                                            healthboard_of_residence))

HB_group_data_bytime<- HB_group_data_bytime%>%
  group_by(time_period, healthboard_of_residence)%>%
  summarise(count_sliccd = sum(count_sliccd), count_denom = sum(count_denom), .groups = 'drop')%>%
  mutate(totalbirth_prevalence = count_sliccd/ count_denom * 10000)

##### ADD live BIRTH PREVALENCE COLUMN TO DATASET 

HB_group_data_bytime <- HB_group_data_bytime%>%
  mutate(totalbirth_prevalence = count_sliccd/count_denom * 10000)

#### change time to 0-21:

HB_group_data_bytime <- transform(HB_group_data_bytime, time_period = time_period - min(time_period))

#### maternal age adjusted by linear time trend 

HB_model_with_time <- glmmTMB(count_sliccd ~ healthboard_of_residence + time_period +offset(log(count_denom)),
                           family = compois, data = HB_group_data_bytime) 

summary (HB_model_with_time)
HB_model_time_summ<- summary(HB_model_with_time)

LB_HB_model_with_time <- data.frame(HB_model_time_summ$coefficients$cond)

setDT(LB_HB_model_with_time, keep.rownames = TRUE)

LB_HB_model_with_time$Estimate<- round(LB_HB_model_with_time$Estimate, 2)

LB_HB_model_with_time$Pr...z..<- round(LB_HB_model_with_time$Pr...z.., 3)

###remove unneeded columns 

LB_HB_model_with_time$Std..Error <- NULL
LB_HB_model_with_time$z.value <- NULL

CI_model_with_time_summary <- data.frame(round(confint(HB_model_with_time), 2))

CI_model_with_time_summary <- CI_model_with_time_summary%>%
  mutate(pRR = exp(Estimate), CI_lower = exp(X2.5..), CI_upper = exp(X97.5..))

CI_model_with_time_summary$pRR <- round(CI_model_with_time_summary$pRR, 2)
CI_model_with_time_summary$CI_lower <- round(CI_model_with_time_summary$CI_lower, 2)
CI_model_with_time_summary$CI_upper <- round(CI_model_with_time_summary$CI_upper, 2)

CI_model_with_time_summary <- CI_model_with_time_summary%>%
  mutate(estimate_CI = paste0(Estimate," (", X2.5.., ", ", X97.5..,")"))%>%
  mutate(pRR_CI = paste0(pRR," (", CI_lower, ", ", CI_upper,")"))

###add in CI column to summary table 

LB_HB_model_with_time<- LB_HB_model_with_time%>%
  mutate(Estimate_CI = CI_model_with_time_summary$estimate_CI)%>%
  mutate(pRR_CI = CI_model_with_time_summary$pRR_CI)

#### plotting HB AND time -----------------------------------------------------


# Predict values
HB_time <- data.frame(time_period = HB_group_data_bytime$time_period, healthboard_of_residence = HB_group_data_bytime$healthboard_of_residence, count_sliccd = HB_group_data_bytime$count_sliccd, count_denom = HB_group_data_bytime$count_denom)
predicted_values <- predict(HB_model_with_time, newdata = HB_time, type = "link", se.fit = TRUE)

# Extract predicted values and standard errors
predicted_means <- predicted_values$fit
standard_errors <- predicted_values$se.fit

# Calculate confidence intervals
z_value <- qnorm(0.975)  # 1.96 for a 95% confidence interval
lower_bound <- exp(predicted_means - z_value * standard_errors)
upper_bound <- exp(predicted_means + z_value * standard_errors)

# Display results
LB_HB_result_df_vglm1 <- data.frame(
  time_period = HB_time$time_period,
  healthboard_of_residence = HB_time$healthboard_of_residence,
  Predicted_Mean = (exp(predicted_means) / HB_group_data_bytime$count_denom) * 10000,
  Lower_Bound = (lower_bound / HB_group_data_bytime$count_denom) * 10000,
  Upper_Bound = (upper_bound / HB_group_data_bytime$count_denom) * 10000,
  birth_prev = (HB_group_data_bytime$count_sliccd / HB_group_data_bytime$count_denom) * 10000  # Include prevs in the result_df
)
print(LB_HB_result_df_vglm1)

HB_group_data_bytime <- HB_group_data_bytime%>%
  mutate(predicted_mean = LB_HB_result_df_vglm1$Predicted_Mean, lower_bound = LB_HB_result_df_vglm1$Lower_Bound, upper_bound = LB_HB_result_df_vglm1$Upper_Bound)

### plotting predicted prevalence per year + SIMD

ggplot(HB_group_data_bytime, aes(x = time_period+2000, y = totalbirth_prevalence, group = healthboard_of_residence))+
  geom_line(aes(colour = healthboard_of_residence), alpha = 0.3, 
            size=0.5)+
  geom_line(aes(colour = healthboard_of_residence, y = predicted_mean))+
  geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound, fill = healthboard_of_residence), alpha = 0.3)+
  geom_point(aes(colour = healthboard_of_residence), size=3)+
  labs(x = "Year of birth",
       y = "Live birth prevalence per 10,000 live births",
       colour = "Maternal health board of residence",
       fill = "Maternal health board of residence")+
  facet_grid(healthboard_of_residence ~  .)+
  ylim(0,30)+
  theme_minimal()+
  theme(text = element_text(size=20))+
  theme(strip.text.y = element_blank() )


## Final plot for manuscript 

hb_PLOT_4 <- ggplot(HB_group_data_bytime, aes(x = time_period+2000, y = totalbirth_prevalence, group = healthboard_of_residence))+
  geom_line(aes(colour = healthboard_of_residence), alpha = 0.3, 
            size=0.5)+
  geom_line(aes(colour = healthboard_of_residence, y = predicted_mean))+
  geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound, fill = healthboard_of_residence), alpha = 0.3)+
  geom_point(aes(colour = healthboard_of_residence), size=3)+
  labs(x = "Year",
       y = "Live birth prevalence per 10,000 live births",
       title = "b",
       colour = "Maternal NHS health board of residence",
       fill = "Maternal NHS health board of residence")+
  facet_wrap(healthboard_of_residence~., ncol = 1)+
  ylim(0,50)+
  theme_minimal()+
  theme(
    axis.title.x = element_text(margin = margin(t = 15)),  # Add space above x-axis title
    axis.title.y = element_text(margin = margin(r = 15))   # Add space to the right of y-axis title
  )+theme(text = element_text(size=20))

hb_PLOT_4


## 6. Sex and live birth prevalence ############################################

### dataset with sex and prevalence per age group 

SEX_group_data <- merged_data_live_births%>%
  group_by(sex_2)%>%
  summarise(count_sliccd = sum(count_sliccd, na.rm = TRUE),
            count_denom = sum(count_denom, na.rm = TRUE))%>%
  ungroup()

###removing unknown category 

SEX_group_data <- SEX_group_data%>%
  filter(!sex_2 == '#N/A')

##### ADD live BIRTH PREVALENCE COLUMN TO DATASET 

SEX_group_data <- SEX_group_data%>%
  mutate(livebirth_prevalence = count_sliccd/count_denom * 10000)

SEX_group_data <- SEX_group_data%>%
  rowwise()%>%
  mutate(ci = list(calculate_CI(count_sliccd, count_denom)),
         prevalence = ci$prevalence,
         ci_lower = ci$ci_lower,
         ci_upper = ci$ci_upper,
         births_prev_ci = paste0(round(prevalence,2), " (", round(ci_lower, 2), "," ,round(ci_upper, 2), ")")) %>%
  ungroup()

##### model of sex and LB prevalence model ----------------------------------------


vglm1<-vglm(SEX_group_data$count_sliccd ~ SEX_group_data$sex_2, poissonff(bred=TRUE), offset=log(count_denom), 
            data=SEX_group_data)

vglm1summary <- summary(vglm1)

show.summary.vglm(vglm1summary)


vglm1summary <- summaryvglm(vglm1)

CI <- data.frame(confint(vglm1))

###create summary table

vglm1summary_tbl <- data.frame (coef(vglm1))

vglm1summary_tbl <- vglm1summary_tbl%>%
  mutate(Lower = CI$X2.5.., Upper = CI$X97.5.. )
## rename table columns 
colnames(vglm1summary_tbl)[1] <- "Coefficients_Est"

## add in coefficients names column 

vglm1summary_tbl <- vglm1summary_tbl%>%
  mutate(Coefficients = c("Intercept" ,"Male"))
#reorder columns in data
vglm1summary_tbl <- vglm1summary_tbl[c("Coefficients", "Coefficients_Est", "Lower", "Upper")]

##add in pRR for each coefficient and CI
vglm1summary_tbl <- vglm1summary_tbl%>%
  mutate(pRR = exp(coef(vglm1)), pRR_CI_lower = exp(Lower), PRR_CI_upper = exp(Upper))
### add in p values

vglm1summary_tbl<- vglm1summary_tbl%>%
  mutate("P value" = coef(summary(vglm1))[,'Pr(>|z|)'])

vglm1summary_tbl$pRR_CI_lower <- round(vglm1summary_tbl$pRR_CI_lower, 2)
vglm1summary_tbl$pRR_CI_upper <- round(vglm1summary_tbl$PRR_CI_upper, 2)

vglm1summary_tbl$pRR <-round(vglm1summary_tbl$pRR,2)

vglm1summary_tbl$`P value` <-round(vglm1summary_tbl$`P value`,3)
vglm1summary_tbl$Coefficients_Est <-round(vglm1summary_tbl$Coefficients_Est,2)
vglm1summary_tbl$Lower <-round(vglm1summary_tbl$Lower,2)
vglm1summary_tbl$Upper <-round(vglm1summary_tbl$Upper,2)

LB_SEX_vglm1summary_tbl <- vglm1summary_tbl%>%
  mutate(estimate_CI = paste0(Coefficients_Est," (", Lower, ", ", Upper,")"))%>%
  mutate(pRR_CI = paste0(pRR," (", vglm1summary_tbl$pRR_CI_lower, ", ", pRR_CI_upper,")"))


##### plotting sex and LB  predicted prevalence -----------------------------

## For birth prevalence 

# Predict values
sex_vglm1 <- data.frame(sex_2 = SEX_group_data$sex_2, count_sliccd = SEX_group_data$count_sliccd, count_denom = SEX_group_data$count_denom)
predicted_values <- predict(vglm1, newdata = sex_vglm1, type = "link", se.fit = TRUE)

# Extract predicted values and standard errors
predicted_means <- predicted_values$fit
standard_errors <- predicted_values$se.fit

# Calculate confidence intervals
z_value <- qnorm(0.975)  # 1.96 for a 95% confidence interval
lower_bound <- exp(predicted_means - z_value * standard_errors)
upper_bound <- exp(predicted_means + z_value * standard_errors)

# Display results
sex_LB_result_df_vglm1 <- data.frame(
  sex_2 = sex_vglm1$sex_2,
  Predicted_Mean = (exp(predicted_means) / SEX_group_data$count_denom) * 10000,
  Lower_Bound = (lower_bound / SEX_group_data$count_denom) * 10000,
  Upper_Bound = (upper_bound / SEX_group_data$count_denom) * 10000,
  birth_prev = (SEX_group_data$count_sliccd / SEX_group_data$count_denom) * 10000  # Include prevs in the result_df
)
print(sex_LB_result_df_vglm1)

#### infant sex and time --------------------------------------------------------

### creating time and sex dataset

SEX_group_data_bytime <- merged_data_live_births%>%
  group_by(time_period, sex_2)%>%
  summarise(count_sliccd = sum(count_sliccd, na.rm = TRUE),
            count_denom = sum(count_denom, na.rm = TRUE))%>%
  ungroup()

###removing unknown category 

SEX_group_data_bytime <- SEX_group_data_bytime%>%
  filter(!sex_2 == '#N/A')

##### ADD live BIRTH PREVALENCE COLUMN TO DATASET 

SEX_group_data_bytime <- SEX_group_data_bytime%>%
  mutate(totalbirth_prevalence = count_sliccd/count_denom * 10000)

#### change time to 0-21:

SEX_group_data_bytime <- transform(SEX_group_data_bytime, time_period = time_period - min(time_period))

#### infant sex adjusted by linear time trend model 

model_with_time <- glmmTMB(count_sliccd ~ sex_2 + time_period +offset(log(count_denom)),
                           family = compois, data = SEX_group_data_bytime) 

summary (model_with_time)
model_time_summ<- summary(model_with_time)

sex_model_with_time <- data.frame(model_time_summ$coefficients$cond)

setDT(sex_model_with_time, keep.rownames = TRUE)

sex_model_with_time$Estimate<- round(sex_model_with_time$Estimate, 2)

sex_model_with_time$Pr...z..<- round(sex_model_with_time$Pr...z.., 3)

###remove unneeded columns 

sex_model_with_time$Std..Error <- NULL
sex_model_with_time$z.value <- NULL

CI_model_with_time_summary <- data.frame(round(confint(model_with_time), 2))

CI_model_with_time_summary <- CI_model_with_time_summary%>%
  mutate(pRR = exp(Estimate), CI_lower = exp(X2.5..), CI_upper = exp(X97.5..))

CI_model_with_time_summary$pRR <- round(CI_model_with_time_summary$pRR, 2)
CI_model_with_time_summary$CI_lower <- round(CI_model_with_time_summary$CI_lower, 2)
CI_model_with_time_summary$CI_upper <- round(CI_model_with_time_summary$CI_upper, 2)

CI_model_with_time_summary <- CI_model_with_time_summary%>%
  mutate(estimate_CI = paste0(Estimate," (", X2.5.., ", ", X97.5..,")"))%>%
  mutate(pRR_CI = paste0(pRR," (", CI_lower, ", ", CI_upper,")"))

###add in CI column to summary table 

LB_SEX_model_with_time<- sex_model_with_time%>%
  mutate(Estimate_CI = CI_model_with_time_summary$estimate_CI)%>%
  mutate(pRR_CI = CI_model_with_time_summary$pRR_CI)


##### plot LB prev of DS and infant sex over time ----------------------------

#### plotting SEX with time

# Predict values
sex_time <- data.frame(time_period = SEX_group_data_bytime$time_period, sex_2 = SEX_group_data_bytime$sex_2, count_sliccd = SEX_group_data_bytime$count_sliccd, count_denom = SEX_group_data_bytime$count_denom)
predicted_values <- predict(model_with_time, newdata = sex_time, type = "link", se.fit = TRUE)

# Extract predicted values and standard errors
predicted_means <- predicted_values$fit
standard_errors <- predicted_values$se.fit

# Calculate confidence intervals
z_value <- qnorm(0.975)  # 1.96 for a 95% confidence interval
lower_bound <- exp(predicted_means - z_value * standard_errors)
upper_bound <- exp(predicted_means + z_value * standard_errors)

# Display results
sex_LB_result_df_vglm1 <- data.frame(
  time_period = sex_time$time_period,
  sex_2 = sex_time$sex_2,
  Predicted_Mean = (exp(predicted_means) / SEX_group_data_bytime$count_denom) * 10000,
  Lower_Bound = (lower_bound / SEX_group_data_bytime$count_denom) * 10000,
  Upper_Bound = (upper_bound / SEX_group_data_bytime$count_denom) * 10000,
  birth_prev = (SEX_group_data_bytime$count_sliccd / SEX_group_data_bytime$count_denom) * 10000  # Include prevs in the result_df
)
print(sex_LB_result_df_vglm1)
SEX_group_data_bytime <- SEX_group_data_bytime%>%
  mutate(predicted_mean = sex_LB_result_df_vglm1$Predicted_Mean, lower_bound = sex_LB_result_df_vglm1$Lower_Bound, upper_bound = sex_LB_result_df_vglm1$Upper_Bound)



## 7. final model - live birth prevalence -----------------------------------

## including maternal age and healthboard 

##### time trasnforming 

merged_data_live_births <- transform(merged_data_live_births, time_period = time_period - min(time_period))

### final model dataset for livebirths 

final_model_data <- merged_data_live_births%>%
  group_by(time_period, age_category, healthboard_of_residence)%>%
  summarise(count_sliccd = sum(count_sliccd, na.rm = TRUE),
            count_denom = sum(count_denom, na.rm = TRUE))%>%
  ungroup()

##removing unknown 
final_model_data <- final_model_data%>%
  filter(!age_category == 'Unknown')

final_model_data <- final_model_data%>%
  filter(!healthboard_of_residence == 'Unknown')

###### grouping island boards into islands label 

final_model_data<- final_model_data%>%
  mutate(healthboard_of_residence = if_else(healthboard_of_residence %in% c('NHS Western Isles', 'NHS Orkney', 'NHS Shetland'),
                                            'NHS Island Boards',
                                            healthboard_of_residence))


LB_final_model_data<- final_model_data%>%
  group_by(time_period, age_category, healthboard_of_residence)%>%
  summarise(count_sliccd = sum(count_sliccd), count_denom = sum(count_denom), .groups = 'drop')%>%
  mutate(totalbirth_prevalence = count_sliccd/ count_denom * 10000)


#### final model

LB_final_model <- glmmTMB(count_sliccd ~ age_category+healthboard_of_residence +offset(log(count_denom)),
                       family = compois, data = LB_final_model_data)

summary (LB_final_model)


FINALMODEL2<- summary(LB_final_model)

final_model_LB <- data.frame(FINALMODEL2$coefficients$cond)

setDT(final_model_LB, keep.rownames = TRUE)

final_model_LB$Estimate<- round(final_model_LB$Estimate, 3)

final_model_LB$Pr...z..<- round(final_model_LB$Pr...z.., 3)

###remove unneeded columns 

final_model_LB$Std..Error <- NULL
final_model_LB$z.value <- NULL

CI_model_with_time_summary <- data.frame(round(confint(LB_final_model), 2))

CI_model_with_time_summary <- CI_model_with_time_summary%>%
  mutate(pRR = exp(Estimate), CI_lower = exp(X2.5..), CI_upper = exp(X97.5..))

CI_model_with_time_summary$pRR <- round(CI_model_with_time_summary$pRR, 2)
CI_model_with_time_summary$CI_lower <- round(CI_model_with_time_summary$CI_lower, 2)
CI_model_with_time_summary$CI_upper <- round(CI_model_with_time_summary$CI_upper, 2)

CI_model_with_time_summary <- CI_model_with_time_summary%>%
  mutate(estimate_CI = paste0(Estimate," (", X2.5.., ", ", X97.5..,")"))%>%
  mutate(pRR_CI = paste0(pRR," (", CI_lower, ", ", CI_upper,")"))

###add in CI column to summary table 

final_model_LB <- final_model_LB%>%
  mutate(Estimate_CI = CI_model_with_time_summary$estimate_CI)%>%
  mutate(pRR_CI = CI_model_with_time_summary$pRR_CI)




