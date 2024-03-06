## Mount Sinai URM PID project
## 1/07/2024
## Examining PID/IEI dx in all children in Mount Sinai system to see if there is a racial disparity in rate of diagnosis 
## This analysis uses all data sent over from MSDW encompassing every record at MSH and w/ CHD pts as controls. 

setwd("~/Documents/Research Projects/URM_PID_project")
library(tidyverse)
library(readr)
library(writexl)
library(readxl)
library(data.table)
library(rstatix)
library(jtools)
library(pubh)
library(sjlabelled)
library(sjPlot)
library(gtsummary)
library(eeptools)
library(tidycensus)
library(multcomp)


# First read in main data file 
rawV2 <- read_excel(path = "./data_files/IEI_CHD_data_v2_raw.xlsx") %>% as.data.frame() %>% setnames(., "Patient", "MRN")

#are there pt duplicates? YES
#remove duplicates
wo_dupV2 <- rawV2 %>% dplyr::distinct(MRN, .keep_all = T) # 156 duplicates in this dataset

# What samples were excluded 
col_cln_V2 <- wo_dupV2 %>% select(., -matches(c("Dates_of_every_inpt_visit_before_spec_ICD10_dx", "Dates_of_every_ED_visit_before_spec_ICD10_dx")))
excluded <- anti_join(wo_dup, col_cln_V2, by = "MRN")

# Which records are duplicated (there are some pts w/ both IEI and CHD)
dupsV2 <- setdiff(rawV2, wo_dupV2) #156 pts have dual dx in the  dataset. 89 of them in IEI cohort

#read in global dx data
raw_dxsV2 <- read_excel(path = "./data_files/IEI_CHD_data_v2_dxs.xlsx") %>% as.data.frame()

# dx's
dx_names_exclud <- semi_join(raw_dxsV2, excluded, by = "MRN") %>% filter(Cohort == "IEI")
dx_table <- table(dx_names_exclud$MRN, dx_names_exclud$Diagnosis_name) %>% as.data.frame() #most common dxs w/ MRN attached
dx_table2 <- table(dx_names_exclud$Diagnosis_name) %>% as.data.frame() #most common dxs

###--------------------Actual Analysis----------------------###

#Will start w/ analysis using the V2 dataset that already has the HIV and other protected classes removed. 

# split up cohorts
IEI <- wo_dupV2 %>% filter(Cohort == "IEI")
CHD <- wo_dupV2 %>% filter(Cohort == "CHD")
IEI_2 <- anti_join(IEI, CHD, by = "MRN") # no change so all dups removed
CHD_2 <- anti_join(CHD, IEI, by = "MRN") # no change so all dups removed

#remove pts who did not identify a race, ethnicity, or sex
unique(IEI_2$Race)
table(IEI$Race)
IEI_2 <- IEI_2 %>% filter(!Race %in% c("Prefer not to say", "PATIENT DECLINED", "UNKNOWN")) #removed 499
IEI_2 <- IEI_2 %>% filter(!Ethnicity %in% c("PATIENT DECLINED")) #removed 1
IEI_2 <- IEI_2 %>% filter(!Gender %in% c("Unknown", "Indeterminate")) #removed 3 
unique(IEI_2$Ethnicity)
table(IEI_2$Ethnicity)
unique(IEI_2$Gender)
table(IEI_2$Gender)
unique(IEI_2$Race_ethnicity_combined)
table(IEI_2$Race_ethnicity_combined)

CHD_2 <- CHD_2 %>% filter(!Race %in% c("Prefer not to say", "PATIENT DECLINED", "UNKNOWN")) 
CHD_2 <- CHD_2 %>% filter(!Ethnicity %in% c("PATIENT DECLINED")) 
CHD_2 <- CHD_2 %>% filter(!Gender %in% c("Unknown")) 
unique(CHD_2$Ethnicity)
table(CHD_2$Ethnicity)
unique(CHD_2$Gender)
table(CHD_2$Gender)
unique(CHD_2$Race_ethnicity_combined)
table(CHD_2$Race_ethnicity_combined)

#aggregate races (IEI)
agg_race_IEI <- IEI_2 %>% 
  mutate(agg_race = case_when(Race %in% c("AFRICAN-AMERICAN", "JAMAICAN", "Black or African-American", "BARBADIAN", "UGANDAN", "HAITIAN", "WEST INDIAN", "TRINIDADIAN", "DOMINICA ISLANDER", "OTHER: NORTH AFRICAN") ~ "black", 
                              Race %in% c("WHITE", "Hispanic/Latino") ~ "white",
                              Race %in% c("CHINESE", "JAPANESE", "BHUTANESE", "Asian", "TAIWANESE", "NEPALESE", "BURMESE", "Other Asian", "ASIAN INDIAN", "BANGLADESHI", "PAKISTANI", "KOREAN", "THAI", "LAOTIAN", "VIETNAMESE") ~ "asian",
                              Race %in% c("Asian (Pacific Islander)", "FILIPINO") ~ "pacific_islander",
                              Race %in% c("Native American", "AMERICAN INDIAN OR ALASKAN") ~ "native_american",
                              TRUE ~ "OTHER"))
table(agg_race_IEI$agg_race)

#aggregate ethnicity (IEI) 
agg_race_eth_IEI <- agg_race_IEI %>% 
  mutate(agg_ethnicity = case_when(Ethnicity %in% c("NON-HISPANIC") ~ "Non_hispanic", 
                                   Ethnicity %in% c("CENTRAL AMERICAN", "GUATEMALAN", "MEXICAN AMERICAN INDIAN", "PERUVIAN", "COLOMBIAN", "HONDURAN", "MEXICANO", "PUERTO RICAN", "COSTA RICAN", "LATIN AMERICAN", "NICARAGUAN", "SALVADORAN", "DOMINICAN", "MEXICAN", "SOUTH AMERICAN", "ECUADORIAN", "MEXICAN AMERICAN", "PANAMANIAN", "SPANIARD") ~ "hispanic",
                                   Ethnicity %in% c("UNKNOWN") ~ "unknown",
  ))
table(agg_race_eth_IEI$agg_ethnicity)

#aggregate races (CHD)
agg_race_CHD <- CHD_2 %>% 
  mutate(agg_race = case_when(Race %in% c("AFRICAN-AMERICAN", "JAMAICAN", "Black or African-American", "BARBADIAN", "UGANDAN", "HAITIAN", "WEST INDIAN", "TRINIDADIAN", "DOMINICA ISLANDER", "OTHER: NORTH AFRICAN", "BLACK OR AFRICAN-AMERICAN", "GUINEAN", "OTHER: WEST AFRICAN", "IVORY COASTIAN", "ETHIOPIAN", "GHANAIAN", "CONGOLESE", "CONGOLESE", "NIGERIAN", "OTHER: EAST AFRICAN", "KENYAN", "MALIAN") ~ "black", 
                              Race %in% c("WHITE", "Hispanic/Latino") ~ "white",
                              Race %in% c("CHINESE", "JAPANESE", "BHUTANESE", "Asian", "TAIWANESE", "NEPALESE", "BURMESE", "Other Asian", "ASIAN INDIAN", "BANGLADESHI", "PAKISTANI", "KOREAN", "THAI", "LAOTIAN", "VIETNAMESE", "OKINAWAN") ~ "asian",
                              Race %in% c("Asian (Pacific Islander)", "FILIPINO", "GUAMANIAN", "NATIVE HAWAIIAN", "Pacific Islander", "PAPUA NEW GUINEAN", "POLYNESIAN", "INDONESIAN", "OTHER PACIFIC ISLANDER") ~ "pacific_islander",
                              Race %in% c("Native American", "AMERICAN INDIAN OR ALASKAN") ~ "native_american",
                              TRUE ~ "OTHER"))
table(agg_race_CHD$agg_race)

#aggregate ethnicity (CHD) 
agg_race_eth_CHD <- agg_race_CHD %>% 
  mutate(agg_ethnicity = case_when(Ethnicity %in% c("NON-HISPANIC") ~ "Non_hispanic", 
                                   Ethnicity %in% c("CENTRAL AMERICAN", "GUATEMALAN", "MEXICAN AMERICAN INDIAN", "PERUVIAN", "COLOMBIAN", "HONDURAN", "MEXICANO", "PUERTO RICAN", "COSTA RICAN", "LATIN AMERICAN", "NICARAGUAN", "SALVADORAN", "DOMINICAN", "MEXICAN", "SOUTH AMERICAN", "ECUADORIAN", "MEXICAN AMERICAN", "PANAMANIAN", "SPANIARD", "ANDALUSIAN", "CRIOLLO", "CENTRAL AMERICAN INDIAN", "CHILEAN", "SOUTH AMERICAN INDIAN", "SPANISH BASQUE", "BOLIVIAN", "ASTURIAN", "ARGENTINEAN", "VENEZUELAN", "PARAGUAYAN", "LA RAZA", "URUGUAYAN", "CATALONIAN", "CANAL ZONE", "CUBAN", "CASTILLIAN") ~ "hispanic",
                                   Ethnicity %in% c("UNKNOWN") ~ "unknown",
  ))
table(agg_race_eth_CHD$agg_ethnicity)

# make column for combined race and ethnicity for white and black pts
table(agg_race_eth_IEI$agg_race, agg_race_eth_IEI$agg_ethnicity)

comb_race_eth_IEI <- agg_race_eth_IEI %>% 
  mutate(comb_race_eth = case_when(agg_race == "white" & agg_ethnicity == "hispanic" ~ "white_hispanic",
                                   agg_race == "white" & agg_ethnicity == "Non_hispanic" ~ "white_non_hispanic",
                                   agg_race == "white" & agg_ethnicity == "unknown" ~ "white_unknown_eth",
                                   agg_race == "black" & agg_ethnicity == "hispanic" ~ "black_hispanic", 
                                   agg_race == "black" & agg_ethnicity == "Non_hispanic" ~ "black_non_hispanic",
                                   agg_race == "black" & agg_ethnicity == "unknown" ~ "black_unknown_eth",
                                   agg_race == "asian" ~ "asian", 
                                   agg_race == "native_american" ~ "native_american",
                                   agg_race == "OTHER" ~ "OTHER",
                                   agg_race == "pacific_islander" ~"pacific_islander")
  )
table(comb_race_eth_IEI$comb_race_eth)

comb_race_eth_CHD <- agg_race_eth_CHD %>% 
  mutate(comb_race_eth = case_when(agg_race == "white" & agg_ethnicity == "hispanic" ~ "white_hispanic",
                                   agg_race == "white" & agg_ethnicity == "Non_hispanic" ~ "white_non_hispanic",
                                   agg_race == "white" & agg_ethnicity == "unknown" ~ "white_unknown_eth",
                                   agg_race == "black" & agg_ethnicity == "hispanic" ~ "black_hispanic", 
                                   agg_race == "black" & agg_ethnicity == "Non_hispanic" ~ "black_non_hispanic",
                                   agg_race == "black" & agg_ethnicity == "unknown" ~ "black_unknown_eth",
                                   agg_race == "asian" ~ "asian", 
                                   agg_race == "native_american" ~ "native_american",
                                   agg_race == "OTHER" ~ "OTHER",
                                   agg_race == "pacific_islander" ~"pacific_islander")
  )
table(comb_race_eth_CHD$comb_race_eth)


# check pt age dist
final_IEI_norm %>%
  ggplot(aes(x=Current_age)) +
  geom_histogram(binwidth=1, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  ggtitle("IEI cohort ages") +
  theme(plot.title = element_text(size=15))

#calculate age of each pt at dx
comb_race_eth_IEI$DOB <- comb_race_eth_IEI$DOB %>% as.Date()
comb_race_eth_IEI$Initial_qualifying_dx_date <- comb_race_eth_IEI$Initial_qualifying_dx_date %>% as.Date()

dx_age_IEI <- comb_race_eth_IEI %>% 
  mutate(dx_age = as.period(interval(start = DOB, end = Initial_qualifying_dx_date))$year)

#pt MRN: 7430968 has wrong DOB. should be: 07-27-1921. This excludes him from this study and will remove
dx_age_IEI <- dx_age_IEI %>% filter(MRN != "7430968")
# 3pts have dx age of 22yrs old, but will keep this as this is likely just a small rounding error from calculating dates

comb_race_eth_CHD$DOB <- comb_race_eth_CHD$DOB %>% as.Date()
comb_race_eth_CHD$Initial_qualifying_dx_date <- comb_race_eth_CHD$Initial_qualifying_dx_date %>% as.Date()

dx_age_CHD <- comb_race_eth_CHD %>% 
  mutate(dx_age = as.period(interval(start = DOB, end = Initial_qualifying_dx_date))$year)
# pt MRN: 7175340 has wrong initial dx date as this occurred before pt was born and EMR has no record of this encounter. Will remove. 
dx_age_CHD <- dx_age_CHD %>% filter(MRN != "7175340")
# 1pt has dx age of 22yrs old, but will keep this as this is likely just a small rounding error from calculating dates



#combine insurances down to medicaid, private, medicare, other, unknown
unique(dx_age_IEI$Insurance)
table(dx_age_IEI$Insurance)
dx_age_insur_IEI <- dx_age_IEI %>% 
  mutate(agg_insur = case_when(Insurance %in% c("AETNA", "EMPIRE-32BJ","EMPIRE-MA", "HORIZON", "BCBSOOS", "MULTIPLAN", "CIGNA", "LOCAL1199", "EMBLEM", "MAGNA", "UHCOHP", "OPTUM", "EMPIRE", "OSCAR", "FIDELIS") ~ "private",
                               Insurance %in% c("HF-MCAID", "METROPLUS", "AFFINITY-MCAID", "METROPLUS-MCAID", "OTHER-MCAID", "EMPIRE-MCAID", "HORIZON-MCAID", "FIDELIS", "NEIGHBOR-MCAID", "FIDELIS-MCAID", "EMBLEM-MCAID", "MEDICAID", "UHCOHP-MCAID") ~ "medicaid",
                               Insurance %in% c("MEDICARE") ~ "medicare",
                               Insurance %in% c("NULL") ~ "unknown",
                               TRUE ~ "other"
  ))
table(dx_age_insur_IEI$agg_insur)

table(dx_age_CHD$Insurance)
dx_age_insur_CHD <- dx_age_CHD %>% 
  mutate(agg_insur = case_when(Insurance %in% c("AETNA", "EMPIRE-32BJ", "BEECH", "AFFINITY", "EMPIRE-MA", "HORIZON", "BCBSOOS", "MULTIPLAN", "CIGNA", "LOCAL1199", "EMBLEM", "MAGNA", "UHCOHP", "OPTUM", "EMPIRE", "OSCAR", "FIDELIS") ~ "private",
                               Insurance %in% c("HF-MCAID", "METROPLUS", "AFFINITY-MCAID", "METROPLUS-MCAID", "OTHER-MCAID", "EMPIRE-MCAID", "HORIZON-MCAID", "FIDELIS", "NEIGHBOR-MCAID", "FIDELIS-MCAID", "EMBLEM-MCAID", "MEDICAID", "UHCOHP-MCAID", "BEACON-MCAID", "WELLCARE-MCAID") ~ "medicaid",
                               Insurance %in% c("MEDICARE", "AMERIGROUP-MA") ~ "medicare",
                               Insurance %in% c("NULL") ~ "unknown",
                               TRUE ~ "other"
  ))
table(dx_age_insur_CHD$agg_insur)

# remove pts who have confounding dxs or don't have true IEI
#first check all ICD10 codes in qualifying dx category 
table(dx_age_insur_IEI$Initial_qualifying_dx_icd10_code) %>% as.data.frame() %>%
  write_xlsx(., "~/Documents/Research Projects/URM_PID_project/MSDW_analy_output/most_com_IEI_dxs.xlsx")
# They all appear correct. 


# calculate socioeconomic status score based on zip codes
census_api_key("da1605649b5e8c84551b28804daa772723c2c1b4")
cen_incom <- get_acs(geography = "zcta", 
                     variables = c(medincome = "B19013_001"), #this is the code for median household income in 12mo
                     year = 2021) %>% select(c("GEOID", "estimate")) %>% rename(Zip_code = GEOID) %>% rename(medincome = estimate)

#link median income based on zip code to our pts
comb_data <- bind_rows(dx_age_insur_CHD, dx_age_insur_IEI)
comb_incom <- left_join(comb_data, cen_incom, by = "Zip_code")
summary(comb_incom$medincome)

#assign score of 1-4 for income level based on medincome for each pt
t <- comb_incom %>% filter(!is.na(medincome)) %>% select(c("medincome", "MRN")) 
quant <- quantile(t$medincome)

q <- (quantile(t$medincome))
t$income_score <- cut(t$medincome, q, include.lowest=TRUE, labels=paste(1:4)) 
comb_incom <- left_join(comb_incom, t, by = "MRN") %>% select(-medincome.y) %>% rename(medincome = medincome.x)

# calculate difference in time to dx after first encounter w/ diagnosing provider
comb_incom$Date_of_first_enc_qualifying_dx_provider <- comb_incom$Date_of_first_enc_qualifying_dx_provider %>% as.Date()

comb_incom <- comb_incom %>%
  mutate(ttd_mo = (interval(ymd(Date_of_first_enc_qualifying_dx_provider), ymd(Initial_qualifying_dx_date))) %/% months(1))
summary(comb_incom$ttd_mo)


##--------Generate final cohorts--------##
final_IEI <- comb_incom %>% filter(Cohort == "IEI") 
final_CHD <- comb_incom %>% filter(Cohort == "CHD")
providers <- table(final_IEI$Initial_qualifying_dx_provider) %>% as.data.frame()

# One of the top IEI dxing providers, NT, is giving unconfirmed PANDAS and other unspecified IEI dx's to pts
# w/o them meeting diagnostic criteria. Will remove all her pts from this analysis. 
nin <- final_IEI %>% filter(Initial_qualifying_dx_provider == "REDACTED")
nin_tab <- table(nin$Initial_qualifying_dx_icd10_code) %>% as.data.frame()

green <- final_IEI %>% filter(Initial_qualifying_dx_provider == "REDACTED")
green_tab <- table(green$Initial_qualifying_dx_icd10_code) %>% as.data.frame()
# One of the top IEI dxing providers, AG, is coding IBD as an unspecified immunodeficiency. 
# Will remove all his pts from this analysis.

redacted <- final_IEI %>% filter(Initial_qualifying_dx_provider == "AC")
redacted_tab <- table(redacted$Initial_qualifying_dx_icd10_code) %>% as.data.frame()
# This provider makes real IEI dx's

redacted <- final_IEI %>% filter(Initial_qualifying_dx_provider == "SB")
redacted_tab <- table(redacted$Initial_qualifying_dx_icd10_code) %>% as.data.frame()
# Dx's OK. 

redacted <- final_IEI %>% filter(Initial_qualifying_dx_provider == "TF")
redacted_tab <- table(redacted$Initial_qualifying_dx_icd10_code) %>% as.data.frame()
# remove these pts as she works w/ NT and has also provided spurious immune deficiency dx's. 15 pts removed

# removing false providers
final_IEI <- final_IEI %>% filter(!Initial_qualifying_dx_provider %in% c("REDACTED")) 

# regenerate final df of both cohorts
final_comb <- bind_rows(final_IEI, final_CHD)


#remove native americans as N is too low
final_comb <- final_comb %>% filter(agg_race != "native_american") 
final_IEI <- final_comb %>% filter(Cohort == "IEI") 
final_CHD <- final_comb %>% filter(Cohort == "CHD") 
saveRDS(final_comb, file = "~/Documents/Research Projects/URM_PID_project/MSDW_analy_output/final_comb.rds")

final_comb <- readRDS("MSDW_analy_output/final_comb.rds")

## 1/20/2024 UPDATE: searching to see if there are transplant pts, CF pts who do not have IEI in the IEI cohort
# Z94 is global ICD10 transplant code. E84.9 is CF code. 
raw_dxsV2 <- read_excel(path = "./data_files/IEI_CHD_data_v2_dxs.xlsx") %>% as.data.frame()

z94_remove <- raw_dxsV2 %>% 
  filter(Cohort == "IEI") %>%
  filter(Icd10_code %in% c("Z94.0", "Z94.4", "Z94.84", "Z94.83", "Z94.82", "Z94.2", "Z94.1", "E84.9")) %>%
  dplyr::distinct(MRN, .keep_all = T) %>%
  anti_join(final_IEI, ., by = "MRN") 
# 151 pts w/ transplant and CF ICD10 codes removed from IEI cohort

# make final cohorts again (1/30/2024)
final_comb <- bind_rows(z94_remove, final_CHD)

final_IEI <- final_comb %>% filter(Cohort == "IEI") 
final_CHD <- final_comb %>% filter(Cohort == "CHD") 
saveRDS(final_comb, file = "~/Documents/Research Projects/URM_PID_project/MSDW_analy_output/final_comb.rds")


## ---- calculate percent of each racial group w/ a dx of an IEI---##
#This measure is most telling of racial disparity 
# This calculation performed as of 1/14/2024 uses totals from all pts <=21yrs old current at MSSM EPIC database. 
# This is not completely accurate as the total number of pediatric pts in each racial category will have changed over time 
# so it will not be completely accurate for older pts in out cohort, but this is a minority of our cohort. 

# Note: the discrepancy in hispanic vs. non-hispanic designation and total of each race differs b/c of individuals w/ no reported and unknown ethnicity

# Totals
all_white <- 204560
all_black <- 84081 
all_native_amer <- 1370
all_pac_island <- 11589
all_asian <- 29756
all_hispanic <- 69235 
all_non_hispanic <- 227456 
black_hispanic <- 4984
black_not_hispanic <- 49418
white_hispanic <- 6861
white_not_hospanic <- 105521

final_IEI %>%
  group_by(agg_race) %>%
  summarise(percent_of_pts=(n()/nrow(final_IEI))*100) %>%
  write_xlsx(., "~/Documents/Research Projects/URM_PID_project/MSDW_analy_output/IEI_per_race.xlsx")

final_CHD %>%
  group_by(agg_race) %>%
  summarise(percent_of_pts=(n()/nrow(final_CHD))*100) %>%
  write_xlsx(., "~/Documents/Research Projects/URM_PID_project/MSDW_analy_output/CHD_per_race.xlsx")


##------ calculate percent of IEI/CHD pts per ethnicity --------##
final_IEI %>%
  group_by(agg_ethnicity) %>%
  summarise(percent_of_pts=(n()/nrow(final_IEI))*100) %>%
  write_xlsx(., "~/Documents/Research Projects/URM_PID_project/MSDW_analy_output/IEI_per_eth.xlsx")

final_CHD %>%
  group_by(agg_ethnicity) %>%
  summarise(percent_of_pts=(n()/nrow(final_CHD))*100) %>%
  write_xlsx(., "~/Documents/Research Projects/URM_PID_project/MSDW_analy_output/CHD_per_eth.xlsx")


##------ calculate percent of IEI/CHD pts per combined race and ethnicity  --------##
final_IEI %>%
  group_by(comb_race_eth) %>%
  summarise(percent_of_pts=(n()/nrow(final_IEI))*100) %>%
  write_xlsx(., "~/Documents/Research Projects/URM_PID_project/MSDW_analy_output/IEI_per_comb_race_eth.xlsx")

final_CHD %>%
  group_by(comb_race_eth) %>%
  summarise(percent_of_pts=(n()/nrow(final_CHD))*100) %>%
  write_xlsx(., "~/Documents/Research Projects/URM_PID_project/MSDW_analy_output/CHD_per_comb_race_eth.xlsx")

#calculate percent of each racial group w/ a dx of an IEI/CHD
#This measure is most telling of racial disparity 
final_IEI %>%
  group_by(agg_race) %>%
  summarise(N=(n())) %>%
  mutate(percent_of_group_w_IEI_dx = case_when(agg_race == "white" ~ (N/all_white)*100,
                                               agg_race == "black" ~ (N/all_black)*100,
                                               agg_race == "asian" ~ (N/all_asian)*100,
                                               agg_race == "native_american" ~ (N/all_native_amer)*100,
                                               agg_race == "pacific_islander" ~ (N/all_pac_island)*100)) %>%
  write_xlsx(., "~/Documents/Research Projects/URM_PID_project/MSDW_analy_output/IEI_per_race_normalized.xlsx")

final_CHD %>%
  group_by(agg_race) %>%
  summarise(N=(n())) %>%
  mutate(percent_of_group_w_CHD_dx = case_when(agg_race == "white" ~ (N/all_white)*100,
                                               agg_race == "black" ~ (N/all_black)*100,
                                               agg_race == "asian" ~ (N/all_asian)*100,
                                               agg_race == "native_american" ~ (N/all_native_amer)*100,
                                               agg_race == "pacific_islander" ~ (N/all_pac_island)*100)) %>%
  write_xlsx(., "~/Documents/Research Projects/URM_PID_project/MSDW_analy_output/CHD_per_race_normalized.xlsx")

# now repeat above w/ ethnicity
final_IEI %>%
  group_by(agg_ethnicity) %>%
  summarise(N=(n())) %>%
  mutate(percent_of_group_w_IEI_dx = case_when(agg_ethnicity == "Non_hispanic" ~ (N/all_non_hispanic)*100,
                                               agg_ethnicity == "hispanic" ~ (N/all_hispanic)*100)) %>%
  write_xlsx(., "~/Documents/Research Projects/URM_PID_project/MSDW_analy_output/IEI_per_eth_normalized.xlsx")  

final_CHD %>%
  group_by(agg_ethnicity) %>%
  summarise(N=(n())) %>%
  mutate(percent_of_group_w_CHD_dx = case_when(agg_ethnicity == "Non_hispanic" ~ (N/all_non_hispanic)*100,
                                               agg_ethnicity == "hispanic" ~ (N/all_hispanic)*100)) %>%
  write_xlsx(., "~/Documents/Research Projects/URM_PID_project/MSDW_analy_output/CHD_per_eth_normalized.xlsx")  

# now repeat above w/ combined race/ethnicity
final_IEI %>%
  group_by(comb_race_eth) %>%
  summarise(N=(n())) %>%
  mutate(percent_of_group_w_IEI_dx = case_when(comb_race_eth == "white_non_hispanic" ~ (N/white_not_hospanic)*100,
                                               comb_race_eth == "white_hispanic" ~ (N/white_hispanic)*100,
                                               comb_race_eth == "black_non_hispanic" ~ (N/black_not_hispanic)*100,
                                               comb_race_eth == "black_hispanic" ~ (N/black_hispanic)*100,
                                               comb_race_eth == "asian" ~ (N/all_asian)*100,
                                               comb_race_eth == "native_american" ~ (N/all_native_amer)*100,
                                               comb_race_eth == "pacific_islander" ~ (N/all_pac_island)*100)) %>%
  write_xlsx(., "~/Documents/Research Projects/URM_PID_project/MSDW_analy_output/IEI_per_race_eth_normalized.xlsx")

final_CHD %>%
  group_by(comb_race_eth) %>%
  summarise(N=(n())) %>%
  mutate(percent_of_group_w_CHD_dx = case_when(comb_race_eth == "white_non_hispanic" ~ (N/white_not_hospanic)*100,
                                               comb_race_eth == "white_hispanic" ~ (N/white_hispanic)*100,
                                               comb_race_eth == "black_non_hispanic" ~ (N/black_not_hispanic)*100,
                                               comb_race_eth == "black_hispanic" ~ (N/black_hispanic)*100,
                                               comb_race_eth == "asian" ~ (N/all_asian)*100,
                                               comb_race_eth == "native_american" ~ (N/all_native_amer)*100,
                                               comb_race_eth == "pacific_islander" ~ (N/all_pac_island)*100)) %>%
  write_xlsx(., "~/Documents/Research Projects/URM_PID_project/MSDW_analy_output/CHD_per_race_eth_normalized.xlsx")


# calculate difference in gender and IEI/CHD dx
final_IEI %>%
  group_by(Gender) %>%
  summarise(percent_of_pts=(n()/nrow(final_IEI))*100) %>%
  write_xlsx(., "~/Documents/Research Projects/URM_PID_project/MSDW_analy_output/IEI_by_sex.xlsx")

final_CHD %>%
  group_by(Gender) %>%
  summarise(percent_of_pts=(n()/nrow(final_CHD))*100) %>%
  write_xlsx(., "~/Documents/Research Projects/URM_PID_project/MSDW_analy_output/CHD_by_sex.xlsx")


##---------------------------------------------------------##



# chi square test 
chisq.test(final_comb$Cohort, final_comb$agg_race)
chisq.test(final_comb$Cohort, final_comb$agg_ethnicity)
chisq.test(final_comb$Cohort, final_comb$Gender)
chisq.test(final_comb$Cohort, final_comb$comb_race_eth)
chisq.test(final_comb$Cohort, final_comb$income_score)
wo_medicare <- final_comb %>% filter(agg_insur != "medicare")
chisq.test(wo_medicare$Cohort, wo_medicare$agg_insur)

#-------------Logistic regression of race/eth normalized dx rates-------------
final_comb_norm <- final_comb 
final_comb_norm$white_non_hisp <- factor(final_comb_norm$white_non_hisp)
final_comb_norm$white_non_hisp <- relevel(final_comb_norm$white_non_hisp, ref = "no")
final_comb_norm$income_score <- relevel(final_comb_norm$income_score, ref = "4")
final_comb_norm$agg_insur <- factor(final_comb_norm$agg_insur)
final_comb_norm$agg_insur <- relevel(final_comb_norm$agg_insur, ref = "private")
final_comb_norm$comb_race_eth <- factor(final_comb_norm$comb_race_eth)
final_comb_norm$comb_race_eth <- relevel(final_comb_norm$comb_race_eth, ref = "white_non_hispanic")
final_comb_norm$Cohort <- factor(final_comb_norm$Cohort)
final_comb_norm$Cohort <- relevel(final_comb_norm$Cohort, ref = "CHD")

final_IEI_norm <- final_comb_norm %>% filter(Cohort == "IEI")
final_CHD_norm <- final_comb_norm %>% filter(Cohort == "CHD")

# Actually need to use cohort as outcome and run GLM w/ CHD included. 
comb_logit <- glm(Cohort ~ dx_age + ttd_mo + agg_insur + income_score + Gender + comb_race_eth, data = final_comb_norm, family = "binomial")
summary(comb_logit) 
# Odds ratios w/ 95% CI
exp(cbind(OR = coef(comb_logit), confint(comb_logit)))
saveRDS(final_comb_norm, file = "~/Documents/Research Projects/URM_PID_project/MSDW_analy_output/final_comb_norm.rds")
#-------------------------------------------------------------------------  


#--------Sub variable analysis in IEI cohort based on race/eth---------
# Compute the analysis of variance
res.aov <- aov(Num_inpatient_enc_after_qualifying_dx ~ comb_race_eth, data = final_IEI_norm)
# Summary of the analysis
summary(res.aov)
TukeyHSD(res.aov)
plot(res.aov, 2)
# Extract the residuals
aov_residuals <- residuals(object = res.aov )
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals )
## data is not normally distributed; need non-parametric test
res.kwt <- kruskal.test(Num_inpatient_enc_after_qualifying_dx ~ comb_race_eth, data = final_IEI_norm)
# multiple comparisons
pairwise.wilcox.test(final_IEI_norm$Num_inpatient_enc_after_qualifying_dx, final_IEI_norm$comb_race_eth,
                     p.adjust.method = "BH")
summary(glht(res.aov, linfct = mcp(comb_race_eth = "Tukey")))  

#calculate median numbers, IQR, and sd for Num_inpatient_enc_after_qualifying_dx
stderror <- function(x) sd(x)/sqrt(length(x))

final_IEI_norm %>%
  group_by(comb_race_eth) %>%
  summarise(
    count = n(),
    mean = mean(Num_inpatient_enc_after_qualifying_dx, na.rm = TRUE),
    sd = sd(Num_inpatient_enc_after_qualifying_dx, na.rm = TRUE),
    IQR = IQR(Num_inpatient_enc_after_qualifying_dx, na.rm = TRUE), 
    SEM = stderror(Num_inpatient_enc_after_qualifying_dx)
  )

#repeat above for ER visits after dx
## data is not normally distributed; need non-parametric test
res.kwt1 <- kruskal.test(Num_ed_enc_after_qualifying_dx ~ comb_race_eth, data = final_IEI_norm)
# multiple comparisons
pairwise.wilcox.test(final_IEI_norm$Num_ed_enc_after_qualifying_dx, final_IEI_norm$comb_race_eth,
                     p.adjust.method = "BH")

final_IEI_norm %>%
  group_by(comb_race_eth) %>%
  summarise(
    count = n(),
    mean = mean(Num_ed_enc_after_qualifying_dx, na.rm = TRUE),
    sd = sd(Num_ed_enc_after_qualifying_dx, na.rm = TRUE),
    IQR = IQR(Num_ed_enc_after_qualifying_dx, na.rm = TRUE),
    SEM = stderror(Num_ed_enc_after_qualifying_dx)
  )

#repeat above for inpt visits before dx
## data is not normally distributed; need non-parametric test
res.kwt2 <- kruskal.test(Num_inpatient_enc_before_qualifying_dx ~ comb_race_eth, data = final_IEI_norm)
# multiple comparisons
pairwise.wilcox.test(final_IEI_norm$Num_inpatient_enc_before_qualifying_dx, final_IEI_norm$comb_race_eth,
                     p.adjust.method = "BH")

final_IEI_norm %>%
  group_by(comb_race_eth) %>%
  summarise(
    count = n(),
    mean = mean(Num_inpatient_enc_before_qualifying_dx, na.rm = TRUE),
    sd = sd(Num_inpatient_enc_before_qualifying_dx, na.rm = TRUE),
    IQR = IQR(Num_inpatient_enc_before_qualifying_dx, na.rm = TRUE),
    SEM = stderror(Num_inpatient_enc_before_qualifying_dx)
  )

#repeat above for ED visits before dx
## data is not normally distributed; need non-parametric test
res.kwt3 <- kruskal.test(Num_ed_enc_before_qualifying_dx ~ comb_race_eth, data = final_IEI_norm)
# multiple comparisons
pairwise.wilcox.test(final_IEI_norm$Num_ed_enc_before_qualifying_dx, final_IEI_norm$comb_race_eth,
                     p.adjust.method = "BH")

final_IEI_norm %>%
  group_by(comb_race_eth) %>%
  summarise(
    count = n(),
    mean = mean(Num_ed_enc_before_qualifying_dx, na.rm = TRUE),
    sd = sd(Num_ed_enc_before_qualifying_dx, na.rm = TRUE),
    IQR = IQR(Num_ed_enc_before_qualifying_dx, na.rm = TRUE),
    SEM = stderror(Num_ed_enc_before_qualifying_dx)
  )



# Calculate mean, SEM and stats for diagnosis age
res.aov5 <- aov(dx_age ~ comb_race_eth, data = final_IEI_norm)
# Summary of the analysis
summary(res.aov5)
TukeyHSD(res.aov5)
plot(res.aov5, 2)
# Extract the residuals
aov_residuals5 <- residuals(object = res.aov5 )
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals5 )
# multiple comparisons
pairwise.wilcox.test(final_IEI_norm$dx_age, final_IEI_norm$comb_race_eth,
                     p.adjust.method = "BH")

#calculate mean and SEM numbers for dx age
final_IEI_norm %>%
  group_by(comb_race_eth) %>%
  summarise(
    count = n(),
    mean = mean(dx_age, na.rm = TRUE),
    sd = sd(dx_age, na.rm = TRUE),
    IQR = IQR(dx_age, na.rm = TRUE), 
    SEM = stderror(dx_age)
  )

# Calculate mean, SEM and stats for time to diagnosis
res.aov6 <- aov(ttd_mo ~ comb_race_eth, data = final_IEI_norm)
# Summary of the analysis
summary(res.aov6)
TukeyHSD(res.aov6)
plot(res.aov6, 2)
# Extract the residuals
aov_residuals6 <- residuals(object = res.aov6 )
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals6 )
# multiple comparisons
pairwise.wilcox.test(final_IEI_norm$dx_age, final_IEI_norm$comb_race_eth,
                     p.adjust.method = "BH")

#calculate mean and SEM numbers for dx age
final_IEI_norm %>%
  group_by(comb_race_eth) %>%
  summarise(
    count = n(),
    mean = mean(ttd_mo, na.rm = TRUE),
    sd = sd(ttd_mo, na.rm = TRUE),
    IQR = IQR(ttd_mo, na.rm = TRUE), 
    SEM = stderror(ttd_mo)
  )


# Statistics and frequency of Income score and Insurance across race/eth in IEI cohort
# chi square test 
test1 <- final_IEI_norm %>% filter(comb_race_eth %in% c("OTHER", "white_non_hispanic", "white_unknown_eth", "black_non_hispanic", "black_unknown_eth", "asian"))
chisq.test(test1$comb_race_eth, test1$income_score)
wo_medicare <- final_IEI_norm %>% filter(agg_insur != "medicare")
chisq.test(wo_medicare$comb_race_eth, wo_medicare$agg_insur)

#Demographic table w/ ethnicity and race seperate for IEI only and by race/eth
final_IEI_norm$comb_race_eth <- as.factor(final_IEI_norm$comb_race_eth)
final_IEI_norm %>%
  dplyr::select(c(comb_race_eth, agg_insur, income_score)) %>%
  mutate(
    comb_race_eth = relevel(comb_race_eth, ref = "white_non_hispanic")
  ) %>%
  copy_labels(final_IEI_norm) %>%
  cross_tbl(by = "comb_race_eth") %>%
  pubh::theme_pubh(2)
#-----------------------------------------------------------------------

# make table
theme_set(sjPlot::theme_sjplot2(base_size = 10))
theme_update(legend.position = "top")
options('huxtable.knit_print_df' = FALSE)
options('huxtable.autoformat_number_format' = list(numeric = "%5.2f"))
knitr::opts_chunk$set(comment = NA)

#Demographic table w/ ethnicity and race seperate
final_comb$Cohort <- as.factor(final_comb$Cohort)
final_comb %>%
  dplyr::select(c(Cohort, agg_race, agg_ethnicity, Gender, agg_insur, income_score)) %>%
  mutate(
    Cohort = relevel(Cohort, ref = "IEI")
  ) %>%
  copy_labels(final_comb) %>%
  cross_tbl(by = "Cohort") %>%
  pubh::theme_pubh(2)

#Demographic table w/ ethnicity and combined
final_comb %>%
  dplyr::select(c(Cohort, comb_race_eth, Gender, agg_insur, income_score)) %>%
  mutate(
    Cohort = relevel(Cohort, ref = "IEI")
  ) %>%
  copy_labels(final_comb) %>%
  pubh::cross_tbl(by = "Cohort") %>%
  pubh::theme_pubh(2) %>%
  write_xlsx(., "~/Documents/Research Projects/URM_PID_project/MSDW_analy_output/cohort_demographics.xlsx")
#-----------------------------------------------------------------------


##-------------Plot of age at diagnosis by race/ethnicity----------------##
ggplot(final_IEI_norm, aes(x = dx_age, colour = as.factor(agg_race), label = as.factor(agg_race))) +
  geom_density() +
  theme_classic() 
# not terribly informative


##-------------Plot time to diagnosis by race/ethnicity----------------##
pp <- ggplot(final_IEI_norm, aes(x = ttd_mo, colour = as.factor(comb_race_eth), label = as.factor(comb_race_eth))) +
  geom_density(size = 1) +
  theme_classic() + xlim(0,25) +
  scale_color_manual(values = c("#955251", "#B565A7", "#009B77", "#DD4124", "#EFC050", "#5B5EA6", "#9B2335", "#0a0000", "#858585"))

pp + theme(legend.text = element_text(size = 9),
           legend.title = element_text(size = 12),
           legend.position = c(0.7, 0.5)) +labs(x = "Time to Diagnosis (mo)", y = "Density",
                                                colour = "Race/Ethnicity") + theme(axis.text = element_text(size = 13))


#-------------------------------------------------------------------------
# check numbers of apts before and after dx distribution
final_IEI_norm %>%
  ggplot( aes(x=Num_ed_enc_before_qualifying_dx, fill = agg_race)) +
  geom_histogram( binwidth=1, alpha=0.9) +
  ggtitle("IEI ED encounters b/f dx") +
  theme(
    plot.title = element_text(size=15)
  )

final_IEI_norm %>%
  ggplot( aes(x=Num_ed_enc_after_qualifying_dx, fill = agg_race)) +
  geom_histogram( binwidth=1, alpha=0.9) +
  ggtitle("IEI ED encounters after dx") +
  theme(
    plot.title = element_text(size=15)
  )

final_IEI_norm %>%
  ggplot( aes(x=Num_inpatient_enc_before_qualifying_dx, fill = agg_race)) +
  geom_histogram( binwidth=1, alpha=0.9) +
  ggtitle("IEI inpt encounters b/f dx") +
  theme(
    plot.title = element_text(size=15)
  )

final_IEI_norm %>%
  ggplot( aes(x=Num_inpatient_enc_after_qualifying_dx, fill = agg_race)) +
  geom_histogram( binwidth=1, alpha=0.9) +
  ggtitle("IEI inpt encounters after dx") +
  theme(
    plot.title = element_text(size=15)
  )

#calculate mean and SEM numbers for current age
final_IEI_norm %>%
  group_by(Cohort) %>%
  summarise(
    count = n(),
    mean = mean(Current_age, na.rm = TRUE),
    sd = sd(Current_age, na.rm = TRUE),
    IQR = IQR(Current_age, na.rm = TRUE), 
    SEM = stderror(Current_age)
  )
