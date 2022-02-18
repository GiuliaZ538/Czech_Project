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

#set all NAs as zeros
df$Date[is.na(df$Date)] <- "Control"


