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

###################### SET UP SCRIPT #################################

############ RUN THIS BEFORE SLICCD analyses #########################


######## 1. Load packages ###########################


library(readxl)      ##for reading in excel files
library(writexl)     ##for writing to excel files 
library(dplyr)       ##for working with 'tidy' data
library(stringr)     ##for working with strings
library(janitor)     ##for rounding 0.5 upwards and clean_names function
library(magrittr)    ##for %<>% operator
library(ggplot2)
library(epiDisplay)
library(reshape2)
library(tidyr)
library(tidyverse)


######### 2. Read in the data ########################


SLiCCD_cohort_2 <- readRDS("/PHI_conf/NIPT_evaluation_anon/Data/cohort_2/sliccd_cohort2_2000-2021_anon_with_eurocat_groups.rds")
#updated SLICCD file with eurocat groups and month of birth 
## filepaths for storing data and outputs 

###################################################
############## 3- define filepath where data is held ############


#define filepath where data held
filepath1 <- "/PHI_conf/NIPT_evaluation_anon/Data/"

#e.g. define a working directly
wd_phd <- "/PHI_conf/NIPT_evaluation_anon/Data/wd_phd/"

objective1 <- "/PHI_conf/NIPT_evaluation_anon/Data/wd_phd/Objective1/"

objective2 <- "/PHI_conf/NIPT_evaluation_anon/Data/wd_phd/Objective2/"

################## 4- format the dataset - extract cohort #####################

########### order dataset by year of preganncy end:###########

sliccd <- SLiCCD_cohort_2 %>%
  arrange(time_period)

######### extract only those with singleton pregnancy and DS #######

sliccd_ds_extract <- sliccd %>%
  filter(ds_singleton_cohort == 1)




