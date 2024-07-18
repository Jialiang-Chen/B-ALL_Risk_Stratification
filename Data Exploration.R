setwd("/Users/chenjialiang/Desktop/MSc AI for Medicine&Medical Research/Research Project")

library(openxlsx)
library(dplyr)
library(survival)
library(survminer)
library(contsurvplot)
library(glmnet)

data_folder = "all_phase2_target_2018_pub"

sample <- read.delim(paste(data_folder, "data_clinical_sample.txt", sep = "/"), skip = 4)
patient <- read.delim(paste(data_folder, "data_clinical_patient.txt", sep = "/"), skip = 4)

#write.xlsx(sample, "sample.xlsx")
#write.xlsx(patient, "patient.xlsx")

# sample count per patient
patient_sample_cnt <- sample |>
  group_by(PATIENT_ID) |>
  summarize(SAMPLE_NUM = n())

hist(patient_sample_cnt$SAMPLE_NUM, labels = TRUE, ylim = c(0, 2000))
count(patient_sample_cnt, SAMPLE_NUM)

# sample data cleaning
sample$MOLECULAR_SUBTYPE[sample$MOLECULAR_SUBTYPE == 'None of above'] <- 'None of the above'
# sample$ETV6_RUNX1_FUSION_STATUS[sample$ETV6_RUNX1_FUSION_STATUS == ''] <- 'Unknown'
# sample$TRISOMY_4_10[sample$TRISOMY_4_10 == ''] <- 'Unknown'
# sample$MLL_STATUS[sample$MLL_STATUS == ''] <- 'Unknown'
# sample$BCR_ABL1_STATUS[sample$BCR_ABL1_STATUS == ''] <- 'Unknown'
# sample$TCF3_PBX1_STATUS[sample$TCF3_PBX1_STATUS == ''] <- 'Unknown'

# sample data check if patients have same results for all samples
sample_ <- select(sample, -SAMPLE_ID)
dup <- sample_[duplicated(sample_), ]
#write.xlsx(dup, "duplicates.xlsx")

# check_same_data <- function(sample_) {
#   df_same <- apply(df[, -1], 1, function(x) all(x == x[1]))
#   return(df_same)
# }
# 
# same_data_by_id <- aggregate(sample_[, -1], by = list(sample_$PATIENT_ID), FUN = check_same_data)


# patient data cleaning
patient_removed <- patient[which(patient$SEX == '' & patient$RACE == '' & patient$ETHNICITY == ''), ] # 6 patients are removed
patient <- patient[-which(patient$SEX == '' & patient$RACE == '' & patient$ETHNICITY == ''), ] # 1545 patients left
patient_removed_os <- patient[which(patient$OS_STATUS == '' | (is.na(patient$OS_MONTHS) & is.na(patient$OS_DAYS))), ] # 10 patients are removed
patient <- patient[-which(patient$OS_STATUS == '' | (is.na(patient$OS_MONTHS) & is.na(patient$OS_DAYS))), ] # 1535 patients left
patient$OS_MONTHS <- ifelse(is.na(patient$OS_MONTHS) & !is.na(patient$OS_DAYS), ceiling(patient$OS_DAYS/30), patient$OS_MONTHS) # convert os days to months if months na and days not na

# merge clinical info with patient data
patient_sample <- unique(sample[, c('PATIENT_ID'
                                    , 'BONE_MARROW_SITE_OF_RELAPSE'
                                    , 'CNS_SITE_OF_RELAPSE'
                                    , 'TESTES_SITE_OF_RELAPSE'
                                    , 'OTHER_SITE_OF_RELAPSE'
                                    , 'ETV6_RUNX1_FUSION_STATUS'
                                    , 'TRISOMY_4_10'
                                    , 'MLL_STATUS'
                                    , 'BCR_ABL1_STATUS'
                                    , 'TCF3_PBX1_STATUS'
                                    , 'DNA_INDEX'
                                    , 'CELL_OF_ORIGIN'
                                    , 'MOLECULAR_SUBTYPE'
                                    , 'ONCOTREE_CODE'
                                    #, 'ANALYSIS_COHORT' # not 1-1 relation
                                    , 'CANCER_TYPE'
                                    , 'CANCER_TYPE_DETAILED'
                                    #, 'TMB_NONSYNONYMOUS' # not 1-1 relation
                                    )]) 

data_combined <- patient |>
  mutate(OS_STATUS_ = as.numeric(substr(OS_STATUS,1,1))) |>
  inner_join(patient_sample, by = "PATIENT_ID") |>
  mutate(RELAPSE = if_else(BONE_MARROW_SITE_OF_RELAPSE == "Yes" | CNS_SITE_OF_RELAPSE == "Yes" | TESTES_SITE_OF_RELAPSE == "Yes" | OTHER_SITE_OF_RELAPSE == "Yes", "Yes"
                           , if_else(BONE_MARROW_SITE_OF_RELAPSE == "" & CNS_SITE_OF_RELAPSE == "" & TESTES_SITE_OF_RELAPSE == "" & OTHER_SITE_OF_RELAPSE == "", "Unknown", "No"))) |>
  mutate(AGE_GROUP = if_else(AGE < 10, "<10", ">=10")) |>
  mutate(MRD_STATUS = if_else(MRD_PERCENT_DAY_29 <= 0.01, "Negative", if_else(is.na(MRD_PERCENT_DAY_29), NA, "Positive"))) |>
  mutate(CNS_STATUS = substr(CNS_STATUS, 1, 5 )) |>
  mutate(DNA_INDEX_LEVEL = if_else(DNA_INDEX < 1, "<1", if_else(DNA_INDEX >= 1.16 & DNA_INDEX <= 1.6, "1.16-1.6", if_else(DNA_INDEX > 1.6, ">1.6", if_else(is.na(DNA_INDEX), NA, "1-1.16")))))

data_combined$MOLECULAR_SUBTYPE[data_combined$MOLECULAR_SUBTYPE == ''] <- 'Unknown'


# keep useful features only
feature_list <- c('PATIENT_ID'
                  , 'SEX'
                  , 'ETHNICITY'
                  , 'RACE'
                  , 'AGE'
                  , 'AGE_GROUP'
                  , 'WBC'
                  , 'CNS_STATUS'
                  , 'TESTICULAR_INVOLVEMENT'
                  , 'MRD_PERCENT_DAY_29'
                  , 'MRD_STATUS'
                  #, 'BM_DAY_8'
                  #, 'BM_DAY_15'
                  , 'BM_DAY_29'
                  #, 'BM_DAY_43'
                  , 'KARYOTYPE'
                  , 'CONGENITAL_ABNORMALITY'
                  , 'PROTOCOL'
                  , 'ALTERNATE_THERAPY'
                  , 'ALTERNATE_THERAPY_OTHER'
                  , 'BONE_MARROW_SITE_OF_RELAPSE'
                  , 'CNS_SITE_OF_RELAPSE'
                  , 'TESTES_SITE_OF_RELAPSE'
                  , 'OTHER_SITE_OF_RELAPSE'
                  , 'RELAPSE'
                  , 'MOLECULAR_SUBTYPE'
                  , 'BCR_ABL1_STATUS'
                  , 'ETV6_RUNX1_FUSION_STATUS'
                  , 'TRISOMY_4_10'
                  , 'MLL_STATUS'
                  , 'TCF3_PBX1_STATUS'
                  , 'DNA_INDEX'
                  , 'DNA_INDEX_LEVEL'
                  , 'OS_MONTHS'
                  , 'OS_STATUS_'
)

data_combined <- data_combined[ , feature_list]

# remove patient if important features missing
data_combined <- data_combined |>
  mutate(MISSING_COUNT = rowSums(is.na(data_combined) | data_combined == ""))

data_combined <- data_combined |>
  filter(MISSING_COUNT < 10
         , !is.na(DNA_INDEX)
         , !is.na(MRD_PERCENT_DAY_29)
         , !is.na(BM_DAY_29)
         , CNS_STATUS != '.'
         , !(BCR_ABL1_STATUS == 'Unknown' & ETV6_RUNX1_FUSION_STATUS == 'Unknown' & TRISOMY_4_10 == 'Unknown' & MLL_STATUS == 'Unknown' & TCF3_PBX1_STATUS == 'Unknown'))


write.xlsx(data_combined, "data_combined.xlsx")

# patient basic clinical info 
count(data_combined, SEX)
count(data_combined, AGE_GROUP)
count(data_combined, RACE)
count(data_combined, ETHNICITY)
count(data_combined, CNS_STATUS)
count(data_combined, TESTICULAR_INVOLVEMENT)


hist(data_combined$AGE, labels = TRUE, ylim = c(0, 500))

hist(data_combined$WBC, labels = TRUE, ylim = c(0, 1500))

# check MRD availability 
sum(is.na(data_combined$MRD_PERCENT_DAY_8))
sum(is.na(data_combined$MRD_PERCENT_DAY_29))
sum(is.na(data_combined$MRD_PERCENT_DAY_43))

hist(data_combined$MRD_PERCENT_DAY_29, labels = TRUE, ylim = c(0, 2000))

sum(data_combined$MRD_PERCENT_DAY_29 <= 0.01) #MRD negative
sum(data_combined$MRD_PERCENT_DAY_29 > 0.01) #MRD positive

count(data_combined, MRD_STATUS)

# check BM availability 
sum(is.na(data_combined$BM_DAY_8))
sum(is.na(data_combined$BM_DAY_15))
sum(is.na(data_combined$BM_DAY_29))
sum(is.na(data_combined$BM_DAY_43))

hist(data_combined$BM_DAY_29, labels = TRUE, ylim = c(0, 2000))

# check relapse
count(data_combined, RELAPSE)
count(data_combined, BONE_MARROW_SITE_OF_RELAPSE)
count(data_combined, CNS_SITE_OF_RELAPSE)
count(data_combined, TESTES_SITE_OF_RELAPSE)
count(data_combined, OTHER_SITE_OF_RELAPSE)

# check protocal(treatment)
count(data_combined, PROTOCOL)

# check congenital abnormality
count(data_combined, CONGENITAL_ABNORMALITY)

# check first event
count(data_combined, FIRST_EVENT)

# molecular subtype count
count(data_combined, MOLECULAR_SUBTYPE)

# cell of origin count
#count(data_combined, CELL_OF_ORIGIN)

# BCR-ABL1 count
count(data_combined, BCR_ABL1_STATUS)

# ETV6-RUNX1 count
count(data_combined, ETV6_RUNX1_FUSION_STATUS)

# TRISOMY of chromosome 4 & 10 count
count(data_combined, TRISOMY_4_10)

# MLL count
count(data_combined, MLL_STATUS)

# TCF3-PBX1 count
count(data_combined, TCF3_PBX1_STATUS)

# DNA index
boxplot(data_combined$DNA_INDEX)
summary(data_combined$DNA_INDEX)

# survival curve - all
cox_fit <- survfit(Surv(OS_MONTHS, OS_STATUS_) ~ 1, data = data_combined)

ggsurvplot(cox_fit, conf.int=TRUE, pval=TRUE, risk.table=TRUE
           , risk.table.height=.3
           , xlab = "Time (Month)"
           , title="Survival Curve - All")

# survival curve - sex
cox_fit <- survfit(Surv(OS_MONTHS, OS_STATUS_) ~ SEX, data = data_combined)

ggsurvplot(cox_fit, conf.int=TRUE, pval=TRUE, risk.table=TRUE
           , risk.table.height=.3
           , xlab = "Time (Month)"
           , title="Survival Curve - Sex")

# survival curve - age group
cox_fit <- survfit(Surv(OS_MONTHS, OS_STATUS_) ~ AGE_GROUP, data = data_combined)

ggsurvplot(cox_fit, conf.int=TRUE, pval=TRUE, risk.table=TRUE
           , risk.table.height=.3
           , xlab = "Time (Month)"
           , title="Survival Curve - Age Group")

# survival curve - race
cox_fit <- survfit(Surv(OS_MONTHS, OS_STATUS_) ~ RACE, data = data_combined)

ggsurvplot(cox_fit, conf.int=TRUE, pval=TRUE, risk.table=TRUE
           , risk.table.height=.3
           , xlab = "Time (Month)"
           , title="Survival Curve - Race")

# survival curve - ethnicity
cox_fit <- survfit(Surv(OS_MONTHS, OS_STATUS_) ~ ETHNICITY, data = data_combined)

ggsurvplot(cox_fit, conf.int=TRUE, pval=TRUE, risk.table=TRUE
           , risk.table.height=.3
           , xlab = "Time (Month)"
           , title="Survival Curve - Ethnicity")

# survival curve - CNS Status
cox_fit <- survfit(Surv(OS_MONTHS, OS_STATUS_) ~ CNS_STATUS, data = data_combined)

ggsurvplot(cox_fit, conf.int=TRUE, pval=TRUE, risk.table=TRUE
           , risk.table.height=.3
           , xlab = "Time (Month)"
           , title="Survival Curve - CNS Status")

# survival curve - Testicular involvement
cox_fit <- survfit(Surv(OS_MONTHS, OS_STATUS_) ~ TESTICULAR_INVOLVEMENT , data = data_combined)

ggsurvplot(cox_fit, conf.int=TRUE, pval=TRUE, risk.table=TRUE
           , risk.table.height=.3
           , xlab = "Time (Month)"
           , title="Survival Curve - Testicular Involvement")

# survival curve - DNA index level
cox_fit <- survfit(Surv(OS_MONTHS, OS_STATUS_) ~ DNA_INDEX_LEVEL, data = data_combined)

ggsurvplot(cox_fit, conf.int=TRUE, pval=TRUE, risk.table=TRUE
           , risk.table.height=.3
           , xlab = "Time (Month)"
           , title="Survival Curve - DNA Index Level")

# survival curve - molecular subtype
cox_fit <- survfit(Surv(OS_MONTHS, OS_STATUS_) ~ MOLECULAR_SUBTYPE, data = data_combined)

ggsurvplot(cox_fit, conf.int=TRUE, pval=TRUE, risk.table=TRUE
           , risk.table.height=.3
           , xlab = "Time (Month)"
           , title="Survival Curve - Molecular Subtype")

# survival curve - BCR-ABL1
cox_fit <- survfit(Surv(OS_MONTHS, OS_STATUS_) ~ BCR_ABL1_STATUS, data = data_combined)

ggsurvplot(cox_fit, conf.int=TRUE, pval=TRUE, risk.table=TRUE
           , risk.table.height=.3
           , xlab = "Time (Month)"
           , title="Survival Curve - BCR-ABL1 Status")

# survival curve - MRD Status
cox_fit <- survfit(Surv(OS_MONTHS, OS_STATUS_) ~ MRD_STATUS, data = data_combined)

ggsurvplot(cox_fit, conf.int=TRUE, pval=TRUE, risk.table=TRUE
           , risk.table.height=.3
           , xlab = "Time (Month)"
           , title="Survival Curve - MRD Status")

# survival curve - congenital abnormality Status
cox_fit <- survfit(Surv(OS_MONTHS, OS_STATUS_) ~ CONGENITAL_ABNORMALITY, data = data_combined)

ggsurvplot(cox_fit, conf.int=TRUE, pval=TRUE, risk.table=TRUE
           , risk.table.height=.3
           , xlab = "Time (Month)"
           , title="Survival Curve - Congenital Abnormality Status")

# survival curve - relapse
cox_fit <- survfit(Surv(OS_MONTHS, OS_STATUS_) ~ RELAPSE, data = data_combined)

ggsurvplot(cox_fit, conf.int=TRUE, pval=TRUE, risk.table=TRUE
           , risk.table.height=.3
           , xlab = "Time (Month)"
           , title="Survival Curve - Relapse Status")

# survival curve - protocol
cox_fit <- survfit(Surv(OS_MONTHS, OS_STATUS_) ~ PROTOCOL, data = data_combined)

ggsurvplot(cox_fit, conf.int=TRUE, pval=TRUE, risk.table=TRUE
           , risk.table.height=.3
           , xlab = "Time (Month)"
           , title="Survival Curve - Protocol")

# cox regression
# cox_res <- coxph(Surv(OS_MONTHS, OS_STATUS_) ~ SEX + AGE + BCR_ABL1_STATUS + ETV6_RUNX1_FUSION_STATUS + TRISOMY_4_10 + MLL_STATUS + TCF3_PBX1_STATUS
#                  , data = data_combined, x = TRUE)

# cox_res <- coxph(Surv(OS_MONTHS, OS_STATUS_) ~ SEX + AGE + WBC + BCR_ABL1_STATUS + HYPODIPLOIDY
#                  , data = data_combined, x = TRUE)

# cox_res <- coxph(Surv(OS_MONTHS, OS_STATUS_) ~ SEX + AGE + MOLECULAR_SUBTYPE
#                  , data = data_combined, x = TRUE)

cox_res <- coxph(Surv(OS_MONTHS, OS_STATUS_) ~ SEX + AGE_GROUP + ETHNICITY+ WBC + CNS_STATUS + BCR_ABL1_STATUS 
                 + ETV6_RUNX1_FUSION_STATUS + TRISOMY_4_10 + MLL_STATUS + TCF3_PBX1_STATUS + DNA_INDEX 
                 + MRD_STATUS + BM_DAY_29 + RELAPSE + PROTOCOL + CONGENITAL_ABNORMALITY
                  , data = data_combined, x = TRUE)

# cox_res <- coxph(Surv(OS_MONTHS, OS_STATUS_) ~ MOLECULAR_SUBTYPE
#                   , data = data_combined, x = TRUE)

cox_res
summary(cox_res)

# check statistics
cox.zph(cox_res)

# check data pattern overtime
ggcoxzph(cox.zph(cox_res), point.size = 1)

# hazard ratio
ggforest(cox_res, data = data_combined)

# some continuous features
# survival curve - age
plot_surv_area(time = "OS_MONTHS"
               , status = "OS_STATUS_"
               , variable = "AGE"
               , data = data_combined
               , model = cox_res
               , title = "Survival Curve - Age"
               , xlab = "Time (Month)"
               , start_color = "lightblue"
               , end_color = "orchid2")

# survival curve - WBC
plot_surv_area(time = "OS_MONTHS"
               , status = "OS_STATUS_"
               , variable = "WBC"
               , data = data_combined
               , model = cox_res
               , title = "Survival Curve - WBC"
               , xlab = "Time (Month)"
               , start_color = "lightblue"
               , end_color = "orchid2")

# survival curve - BM DAY 29
plot_surv_area(time = "OS_MONTHS"
               , status = "OS_STATUS_"
               , variable = "BM_DAY_29"
               , data = data_combined
               , model = cox_res
               , title = "Survival Curve - BM_DAY_29"
               , xlab = "Time (Month)"
               , start_color = "lightblue"
               , end_color = "orchid2")

# survival curve - DNA INDEX
plot_surv_area(time = "OS_MONTHS"
               , status = "OS_STATUS_"
               , variable = "DNA_INDEX"
               , data = data_combined
               , model = cox_res
               , title = "Survival Curve - DNA Index"
               , xlab = "Time (Month)"
               , start_color = "lightblue"
               , end_color = "orchid2")

# survival curve - MRD_PERCENT_DAY_29
plot_surv_area(time = "OS_MONTHS"
               , status = "OS_STATUS_"
               , variable = "MRD_PERCENT_DAY_29"
               , data = data_combined
               , model = cox_res
               , title = "Survival Curve - MRD_PERCENT_DAY_29"
               , xlab = "Time (Month)"
               , start_color = "lightblue"
               , end_color = "orchid2")



