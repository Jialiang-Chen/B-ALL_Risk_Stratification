setwd("/Users/chenjialiang/Desktop/MSc AI for Medicine&Medical Research/Research Project")
set.seed(22)

library(openxlsx)
library(dplyr)
library(tidyr)
library(GenVisR)
library(data.table)

data_folder = "all_phase2_target_2018_pub"

# clinical data preparation --------------------------------------------------------

# read sample and patient data
sample <- read.delim(paste(data_folder, "data_clinical_sample.txt", sep = "/"), skip = 4)
patient <- read.delim(paste(data_folder, "data_clinical_patient.txt", sep = "/"), skip = 4)

# sample data cleaning
sample$MOLECULAR_SUBTYPE[sample$MOLECULAR_SUBTYPE == 'None of above'] <- 'None of the above'

# patient data cleaning
patient_removed <- patient[which(patient$SEX == '' & patient$RACE == '' & patient$ETHNICITY == ''), ] # 6 patients are removed
patient <- patient[-which(patient$SEX == '' & patient$RACE == '' & patient$ETHNICITY == ''), ] # 1545 patients left
patient_removed_os <- patient[which(patient$OS_STATUS == '' | (is.na(patient$OS_MONTHS) & is.na(patient$OS_DAYS))), ] # 10 patients are removed
patient <- patient[-which(patient$OS_STATUS == '' | (is.na(patient$OS_MONTHS) & is.na(patient$OS_DAYS))), ] # 1535 patients left
patient$OS_MONTHS <- ifelse(is.na(patient$OS_MONTHS) & !is.na(patient$OS_DAYS), ceiling(patient$OS_DAYS/30), patient$OS_MONTHS) # convert os days to months if months na and days not na

# merge clinical info in sample data with patient data
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

# keep useful features for important feature missing check
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

# write.xlsx(data_combined, "data_combined.xlsx")

# keep useful features only
feature_list <- c('PATIENT_ID'
                  , 'SEX'
                  , 'ETHNICITY'
                  , 'AGE_GROUP'
                  , 'WBC'
                  , 'CNS_STATUS'
                  , 'MRD_STATUS'
                  , 'BM_DAY_29'
                  , 'PROTOCOL'
                  , 'RELAPSE'
                  , 'MOLECULAR_SUBTYPE'
                  , 'BCR_ABL1_STATUS'
                  , 'ETV6_RUNX1_FUSION_STATUS'
                  , 'TRISOMY_4_10'
                  , 'MLL_STATUS'
                  , 'TCF3_PBX1_STATUS'
                  , 'DNA_INDEX_LEVEL'
                  , 'OS_MONTHS'
                  , 'OS_STATUS_'
)

clinical <- data_combined[ , feature_list]


# mutation plot -----------------------------------------------------------
mutation0 <- read.delim(paste0(data_folder, "/data_mutations.txt"))
mutation0$Tumor_Sample_Barcode <- substr(mutation0$Tumor_Sample_Barcode, 1, 16)
mutation0 <- filter(mutation0, !Tumor_Sample_Barcode %in% c("TARGET-10-PAPEWB", "TARGET-10-PANWHW", "TARGET-10-PAPCUR", "TARGET-10-PAPIGX"))

patient_list <- data.frame(data_combined$PATIENT_ID)
colnames(patient_list) <- "Tumor_Sample_Barcode"
mutation0 <- left_join(patient_list, mutation0, by = "Tumor_Sample_Barcode")
#write.table(mutation0, file = "mutation_plot.txt", sep = "\t", row.names = TRUE, col.names = NA)

#mutation <- MutationAnnotationFormat("mutation_plot.txt", version = 2.4, verbose = FALSE)
mutation <- select(mutation0, Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification)
setnames(mutation, c("sample", "gene", "mutation"))

classes <- c("Nonsense_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del", "In_Frame_Del", "Translation_Start_Site", "Nonstop_Mutation", "5’Flank", "3’UTR", "Targeted_Region")

myHierarchy <- data.table("mutation" = classes,
                          color=c("#FF0000", "#00A08A", "#F2AD00", "#F98400", "#5BBCD6", "#046C9A", "#D69C4E", "#000000", "#446455"))


clinical_plot <- clinical |>
  select(PATIENT_ID, AGE_GROUP, CNS_STATUS, WBC, MRD_STATUS, BM_DAY_29, PROTOCOL, RELAPSE, MOLECULAR_SUBTYPE
         , BCR_ABL1_STATUS, ETV6_RUNX1_FUSION_STATUS, TRISOMY_4_10, MLL_STATUS
         , TCF3_PBX1_STATUS, DNA_INDEX_LEVEL, OS_STATUS_) |>
  mutate(AGE_GROUP = paste0("AGE_GROUP:", AGE_GROUP)) |>
  mutate(CNS_STATUS = paste0("CNS_STATUS:", CNS_STATUS)) |>
  mutate(WBC = if_else(WBC<50, "WBC:<50", "WBC:>=50")) |>
  mutate(MRD_STATUS = paste0("MRD_STATUS:", MRD_STATUS)) |>
  mutate(BM_DAY_29 = if_else(BM_DAY_29<5, "BM_DAY_29:<5", if_else(BM_DAY_29<50, "BM_DAY_29:5-50", "BM_DAY_29>50"))) |>
  mutate(PROTOCOL = paste0("PROTOCOL:", PROTOCOL)) |>
  mutate(RELAPSE = paste0("RELAPSE:", RELAPSE)) |>
  mutate(MOLECULAR_SUBTYPE = paste0("MOLECULAR_SUBTYPE:", MOLECULAR_SUBTYPE)) |>
  mutate(BCR_ABL1_STATUS = paste0("BCR_ABL1_STATUS:", BCR_ABL1_STATUS)) |>
  mutate(ETV6_RUNX1_FUSION_STATUS = paste0("ETV6_RUNX1_FUSION_STATUS:", ETV6_RUNX1_FUSION_STATUS)) |>
  mutate(TRISOMY_4_10 = paste0("TRISOMY_4_10:", TRISOMY_4_10)) |>
  mutate(MLL_STATUS = paste0("MLL_STATUS:", MLL_STATUS)) |>
  mutate(TCF3_PBX1_STATUS = paste0("TCF3_PBX1_STATUS:", TCF3_PBX1_STATUS)) |>
  mutate(DNA_INDEX_LEVEL = paste0("DNA_INDEX_LEVEL:", DNA_INDEX_LEVEL)) |>
  mutate(OS_STATUS_ = if_else(OS_STATUS_ == 0, "0:LIVING", "1:DECEASED")) 
  
clinical_plot = rename(clinical_plot, sample = PATIENT_ID, OS_STATUS = OS_STATUS_)

myClinicalColors <- c("1:DECEASED"="#F2AD00"
                      , "AGE_GROUP:>=10"="#F8BBD0"
                      , "CNS_STATUS:CNS 3"="#C5CAE9"
                      , "CNS_STATUS:CNS 2"="#7CB342"
                      , "WBC:>=50"="#FFD54F"
                      , "MRD_STATUS:Positive"="#5BBCD6"
                      , "BM_DAY_29:5-50"="#9FE2BF"
                      , "BM_DAY_29:>50"="#FA8072"
                      , "PROTOCOL:9906"="#DFFF00"
                      , "PROTOCOL:AALL0232"="#CCC591"
                      , "PROTOCOL:AALL0331"="#F48FB1"
                      , "PROTOCOL:AALL0434"="#E6EE9C"
                      , "RELAPSE:Yes"="#CCCCFF"
                      , "MOLECULAR_SUBTYPE:BCR-ABL1"="#E6EE9C"
                      , "MOLECULAR_SUBTYPE:ETV6-RUNX1"="#EF9A9A"
                      , "MOLECULAR_SUBTYPE:Hypodiploid"="#80CBC4"
                      , "MOLECULAR_SUBTYPE:MLL-Rearranged"="#558B2F"
                      , "BCR_ABL1_STATUS:Positive"="#FFEE58"
                      , "ETV6_RUNX1_FUSION_STATUS:Positive"="#85D4E3"
                      , "TRISOMY_4_10:Positive"="#EF9A9A"
                      , "MLL_STATUS:Positive"="#E1BEE7"
                      , "TCF3_PBX1_STATUS:Positive"="#FFF59D"
                      , "DNA_INDEX_LEVEL:<1"="#EF5350"
                      , "DNA_INDEX_LEVEL:1.16-1.6"="#AED581"
                      , "DNA_INDEX_LEVEL:>1.6"="#80DEEA")


clinicalData <- Clinical(inputData = clinical_plot, palette = myClinicalColors, legendColumns = 2)

#plotData <- Waterfall(mutation, geneMax = 50, clinical = clinicalData)
plotData <- Waterfall(mutation, geneMax = 50, mutationHierarchy = myHierarchy, clinical = clinicalData
                      , sampleNames = F, labelSize = 1)

pdf(file="Mutation Plot.pdf", height=24, width=24)
drawPlot(plotData)
dev.off()
