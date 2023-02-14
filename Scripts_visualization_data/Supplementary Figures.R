###Supplementary Figures 

library(reshape2)
library(tidyverse)
library(data.table)
library(dplyr)
library(ggplot2)
library(forcats)
library (readr)
library(tidytext)

#MawSpecies_BARPLOT_meanlength
metaDMG_aMaw_all <- read_csv("/path_to_metaDMG_output/VM_metaDMG-aMaw.csv")
subset100_species <- metaDMG_aMaw_all %>% filter(N_reads>100, grepl("\\bspecies\\b", tax_rank))

# Magic-palette
subset100_species$sample <- as.factor(subset100_species$sample)
length(levels(subset100_species$sample))    # that's 11 levels 

# Create Color Palette
magic.palette <- c("#9A32CD", "#00CD00", "#050505", "#CD3333", "#008B00", "#1C86EE", "#00EEEE", "#EEB422", "#A3A3A3", "#00008B", "#EE6A50")    # defining 7 colours
names(magic.palette) <- levels(subset100_species$sample)                                                    

#MawSpecies_Damage
# counting number of unique samples for number of different symbols
f <- length(unique(subset100_species$sample))
#plotting  D_max vs lambda_LR
p1 <- ggplot(subset100_species, aes(y=Bayesian_D_max, x=Bayesian_z)) + geom_point(aes(size=N_reads,col=sample)) +
  scale_colour_manual(values = magic.palette,name = "Sample") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust =1)) +
  scale_x_continuous(breaks =seq(-2,10, by=2), limits=c(-2,10))+ 
  ggtitle("Damage vs BayesianZ") + xlab("Bayesian Z") + ylab("D_max")
p1 + theme_minimal()
#+ scale_x_continuous(breaks =seq(0,150, by=10), limits=c(0,150)) + scale_y_continuous(breaks =seq(0,0.6, by=0.05))

#MawSpecies_Barplot_Dmax
Depth2 <- read.csv ("/path_to_metadata/VM_Sediment_context.csv", sep=";") 

subset100_species$new <- Depth2$depth_cm[match(subset100_species$sample, Depth2$sample_ID)]
names(subset100_species)[names(subset100_species) == 'new'] <- 'Depth'

write.csv(subset100_species, file = "subset100_species.csv")
subset100_species$Depth <- as.factor(subset100_species$Depth)

#Mean_Length
p6d <- subset100_species %>%
  mutate(Depth = fct_relevel(Depth,
                             "236","222","196","185","151","134","116","107","64", "26", "18")) %>%
  ggplot(aes(x=mean_L, y=Depth)) + geom_boxplot(aes(x=mean_L, y=Depth, colour=sample))+ scale_colour_manual(values = magic.palette,name = "Sample") + 
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#+ scale_x_continuous(breaks =seq(0,0.4, by = 0.10), limits=c(0,0.4))
p6d 

#D-max
p7d <- subset100_species %>%
  mutate(Depth = fct_relevel(Depth,
                             "236","222","196","185","151","134","116","107","64", "26", "18")) %>%
  ggplot(aes(x=Bayesian_D_max, y=Depth)) + geom_boxplot(aes(x=Bayesian_D_max, y=Depth, colour=sample))+ scale_colour_manual(values = magic.palette,name = "Sample") + 
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p7d 

grid.arrange(p6d, p7d, 
             ncol = 2, nrow = 1)