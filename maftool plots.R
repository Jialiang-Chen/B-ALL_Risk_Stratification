library(maftools)
library(dplyr)
library(tidyr)

setwd("/Users/chenjialiang/Desktop/MSc AI for Medicine&Medical Research/Research Project")

patient_list <- clinical$PATIENT_ID

mutation <- read.delim(paste(data_folder, "data_mutations.txt", sep = "/"))

mutation <- mutation |>
  mutate(PATIENT_ID = substr(Tumor_Sample_Barcode, 1, 16)) |>
  mutate(SAMPLE_ORDER = as.numeric(substr(Tumor_Sample_Barcode, 18, nchar(Tumor_Sample_Barcode)))) |>
  group_by(PATIENT_ID) |>
  filter(SAMPLE_ORDER == min(SAMPLE_ORDER)) |>
  select(-SAMPLE_ORDER) 

mutation <- mutation[mutation$PATIENT_ID %in% patient_list, ]

maf_object <- read.maf(maf = mutation)

plotmafSummary(maf = maf_object, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)

oncoplot(maf = maf_object, top = 10)
