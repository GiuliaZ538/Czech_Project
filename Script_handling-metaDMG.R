#############Script for handling output file from metaDMG ############

#Packages
library(tidyverse)
library(scales)
library(vegan)
library(rioja)
library(readxl)
library(ggplot2)
library(dplyr)

#Loading files
output_data <- read_csv("~/path-file-location/output.csv")
context_data <- read.csv ("~/path-file-location/context_data.csv", sep = ";")
depth_data <- read.csv ("~/path-file-location/depth_data.csv", sep = ";")

#Changing sample name to short format 
df <- output_data
df$sample[df$sample == "long_name_sample"] <- "short_name"

#Merge context_data and depth_data with dataframe(adding new column)
df$new <- depth_data$Depth_cm[match(df$sample, depth_data$Sample_ID)]
names(df)[names(df) == 'new'] <- 'Depth'

df$new <- context_data$Date[match(df$sample, context_data$Sample_ID)]
names(df)[names(df) == 'new'] <- 'Date'

##Adding factor for filtering (to keep updated based on the variable used)
output1_data <- df

output2_data <- output1_data %>%   # overwriting our data frame
  mutate(lambda_LR_new =   # creating our new column
           ifelse(lambda_LR > 10, "LR>10", "LR<10"))

##Adding grouping factor on Kingdoms
output3_data <- output2_data %>% 
  mutate(Kingdom =   # creating our new column
           case_when(grepl("Viridiplantae", tax_path) ~ "Viridiplantae",
                     grepl("2:Bacteria:",tax_path) ~ "Bacteria",
                     grepl("Metazoa",tax_path) ~ "Metazoa",
                     grepl("Fungi", tax_path) ~ "Fungi",
                     grepl("Archaea", tax_path) ~ "Archaea"))

#set NAs
output3_data[is.na(output3_data)] <- "Control"
output3_data$Kingdom[output3_data$Kingdom == "Control"] <- "Indet"

#Subset dataframe to show only genus, filtering reads
subset_genus <- output3_data %>% filter(D_max > 0.000, N_reads > 20, grepl("\\bgenus\\b", tax_rank))

#Other subsets for parsing reads, for ex. plants, animals, filtering D_max, lambda_LR.
subset1 <- output3_data %>% filter(D_max > 0.05, lambda_LR > 10, grepl("Viridiplantae", tax_path))
subset2 <- output3_data %>% filter(D_max > 0.02, N_reads > 100, lambda_LR > 5, grepl("Metazoa", tax_path))
subset3 <- output3_data %>% filter (N_reads > 50, grepl("Bacteria", Kingdom))
