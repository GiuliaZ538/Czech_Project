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

df$new <- context_data$Material[match(df$sample, context_data$Sample_ID)]
names(df)[names(df) == 'new'] <- 'Material'

##Adding factor that correspond to the threshold applied (to keep updated based on the variable used)
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

#Subset dataframe to show only genus level, filtering reads
subset_genus <- output3_data %>% filter(D_max > 0.000, N_reads > 20, grepl("\\bgenus\\b", tax_rank))

#Other subsets for filtering the data, for ex. plants, animals, filtering D_max, lambda_LR...
subset1 <- output3_data %>% filter(D_max > 0.05, lambda_LR > 10, grepl("Viridiplantae", tax_path))
subset2 <- output3_data %>% filter(D_max > 0.02, N_reads > 100, lambda_LR > 5, grepl("Metazoa", tax_path))
subset3 <- output3_data %>% filter (N_reads > 50, grepl("Bacteria", Kingdom))

##Overview of the data D-max vs lambda_LR
# counting number of unique samples for number of different symbols
f <- length(unique(subset1$sample))
#plotting  D_max vs lambda_LR
p1 <- ggplot(subset1, aes(y=D_max, x=lambda_LR)) + geom_point(aes(col=tax_name, size=N_reads,shape=sample)) 
+ scale_shape_manual(values=1:f) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust =1)) 
+ scale_x_continuous(breaks =seq(0,150, by=10), limits=c(0,150)) + scale_y_continuous(breaks =seq(0,0.6, by=0.05))
+ ggtitle("Damage vs LR_output3") +
xlab("likelihood ratio (lambda_LR) for fitting the beta binomial 'DNA damage' model") + ylab("D_max")

##Plotting Damage_vs_depth
p2 <- ggplot(data=subset1, aes(x=D_max, y=depth))
  + geom_point(aes(colour = factor(tax_name), size = N_reads))
  + scale_x_continuous(breaks =seq(0,0.6, by = 0.05), limits=c(0,0.6))
  + scale_y_reverse(breaks=rev(subset1$depth), limits = c(460,0)) #limits of the depth_data
 + ggtitle("Damage vs Depth") +
  xlab("D-max") + ylab("Depth")

##Plotting distribution mean_L vs damage
p3 <- ggplot(subset1, aes(y=D_max, x=mean_L)) + geom_point(aes(col=sample, size=N_reads)) 
+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust =1)) 
+ scale_y_continuous(breaks =seq(0,0.6, by=0.05))
+ ggtitle("Damage vs mean length") +
  xlab("Mean length (BP)") + ylab("D-max")

##Checking presence Taxa per layer and their damage factor (using factor=lambda_LR >10)
p5 <- ggplot(subset1, aes(y=Material, x=tax_name)) + geom_point(aes(x=tax_name, y=Material, col=lambda_LR_new, size=N_reads)) + geom_line(aes(group = tax_name))
p5 + theme_minimal() + ggtitle("Presence of taxa and their damage") + xlab ("Taxa") + ylab("Material")

###Visualization of datasets with box plot/violin plot 
#Plotting damage per sample
p1d <- ggplot(subset1, aes(x=D_max, y=sample, fill=sample)) + geom_violin(width=1.4) + geom_boxplot(width=0.1, color="grey", alpha=0.2) + +theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust =1)) + ggtitle("Damage distribution per sample") +
  xlab("D_max") + ylab("Depth") + scale_y_discrete(limits=rev)

#Plotting damage per taxa
p2d <-ggplot(dm2, aes(x=D_max, y=tax_name)) + geom_boxplot() + theme_minimal() + ggtitle("Damage per taxa") + ylab("Taxa") + xlab("D-max")

#Damage vs Material (factor=Kingdom)
p3d <- subset1 %>%
  mutate(Material = fct_relevel(Material,
                                "laminated-clay-sand","clayey-gyttja", "impure-sand", "Control")) %>%
  ggplot(aes(x=D_max, y=Material)) + geom_boxplot(aes(x=D_max, y=Material,fill=Kingdom))+ geom_jitter()
p3d + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#changing factor with the threshold applied (fill=lambda_LR_new) it is possible to distinguish damaged-not damaged reads.

##Plotting mean_L per sample
p1m <- ggplot(subset1, aes(x=mean_L, y=sample, fill=sample)) + geom_violin(width=1.4) + geom_boxplot(width=0.1, color="grey", alpha=0.2) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust =1)) + ggtitle("Meanlength distribution per sample") +
  xlab("Mean_L") + ylab("Depth") + scale_y_discrete(limits=rev)

#Mean length per taxa 
p2m <-ggplot(dm2, aes(x=mean_L, y=tax_name)) + geom_boxplot() + theme_minimal() + ggtitle("Mean read lengths per taxa") + ylab("Taxa") + xlab("Mean Length (Bp)")

#Mean length vs Material (factor=Kingdom)
p3m <- subset1 %>%
  mutate(Material = fct_relevel(Material,
                                "laminated-clay-sand","clayey-gyttja", "impure-sand", "Control")) %>%
  ggplot(aes(x=mean_L, y=Material)) + geom_boxplot(aes(x=mean_L, y=Material,fill=Kingdom)) + geom_jitter()
p3m + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#changing factor with the threshold applied (fill=lambda_LR_new) it is possible to distinguish damaged-not damaged reads.

##Pollen diagram for plant data
#Wide table with N_reads as variable
data_wide1_subset1 <- dcast(subset1, depth ~ tax_name, value.var="N_reads") #N_reads each depth(=each sample)
data_wide2_subset1 <- dcast(subset1, Material ~ tax_name, value.var="N_reads", sum) #sum of reads each material type (includes more samples)

n <- ncol(data_wide1_subset1)
dw1 <- data_wide1_subset1[,2:n]
dw1[is.na(dw1)] <-0 

dw1_perc <- prop.table(data.matrix(dw1), margin=1)*100 # makes proportion table, margin 1 (per each depth)
dw2_perc <- prop.table(data.matrix(dw1), margin=2)*100 # makes proportion table, margin 2 (per each taxa)
rowSums(dw1_perc) # should give 100 for each row
colSums(dw2_perc) # should give 100 for each column

#Strat.plot() function using package "rioja" 
y.scale <- c(50,150,200,250,325,460) #setting y axis by depths
p.col <- c("forestgreen") # setting color bars

#Plot with bars for each row
pol.plot1 <- strat.plot(dw1_perc, yvar=y.scale, y.tks=y.scale, y.rev=TRUE, plot.line=FALSE, plot.poly=TRUE, plot.bar=TRUE, col.bar="black", col.poly=p.col, col.poly.line="black", scale.percent=TRUE, wa.orde="topleft", xSpace=0.001, x.pc.lab=5, x.pc.omit0=TRUE, las=2)

#Filled plot with bars
pol.plot2 <- strat.plot(dw1_perc, yvar=y.scale, y.tks=y.scale, y.rev=TRUE, plot.line=FALSE, plot.poly=FALSE, plot.bar=TRUE, col.bar=p.col, lwd.bar=10, scale.percent=TRUE, wa.orde="topleft", xSpace=0.01, x.pc.lab=5, x.pc.omit0=5, las=2)

#Exaggerated curve to visualise low abundance taxa (only for fillen plot with bars)
# Specify which pollen types
ex <- c(rep(TRUE, times=3), rep(FALSE, times=2), TRUE, rep(FALSE, times=2), rep(TRUE, times=8), FALSE, rep(TRUE, times=2), FALSE, FALSE, TRUE)

# Filled plot with exaggeration curve
pol.plot <- strat.plot(dw1_perc, yvar=y.scale, y.tks=y.scale, y.rev=TRUE, plot.line=FALSE, plot.poly=TRUE,plot.bar=TRUE, col.bar="black", col.poly=p.col, col.poly.line="black", scale.percent=TRUE, wa.orde="topleft", xSpace=0.01, x.pc.lab=TRUE, x.pc.omit0=TRUE,las=2, exag=ex, col.exag="auto", exag.alpha=0.50, exag.mult=10)

