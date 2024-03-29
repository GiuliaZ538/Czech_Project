library(dplyr)
library(readr)
library(ggplot2)
library(tidyverse)
library(scales)
library(vegan)
library(rioja)
library(readxl)
library(yaml)
library(data.table)
library(gridExtra)
library(forcats)
library(tidytext)

#####Supplementary Figures 2
##Rarefaction with Vegan package
#make wide table from meaDMG output
sample <- metaDMGoutput$sample
tax_id <- metaDMGoutput$tax_name
reads<- metaDMGoutput$N_reads

raw <- data.frame(tax_id, sample, reads)

data_wide <- dcast(raw, tax_id ~ sample, value.var = "reads", fun.aggregate = sum)
head(data_wide) 

#Parsing reads for median of reads or median of sequencing depth
k <- ncol(data_wide)
dd = data_wide[,2:k]
rownames(dd) <- data_wide$tax_id

colSums(dd)
col <- median.default(colSums(dd))/2 #median sequencing depth
row <- median.default(rowSums(dd))/2 #median reads

drop <- colSums(dd) > col
dd1 <- dd[c(drop)]
dd1a <- as.data.frame(t(dd1))
#drop4 <- colSums(dd1a) > row ##drop the median reads each row
#dd2 <- dd1a[c(drop4)]
#dd2a <- as.data.frame(t(dd2))

#Data_t1=as.data.frame(t(dd100a)) #transpose columns-rows for rarecurve
S <- specnumber(dd1a) # observed number of species
raremax <- min(rowSums(dd1a)) #rarefy the sample counts to he minimum sample count
raremax
Srare <- rarefy(dd1a, raremax)

plot(S, Srare, xlab = "Observed No. of Genuses", ylab = "Rarefied No. of Genuses")
abline(0, 1)
out <- rarecurve(dd22a, step = 20, sample = raremax, col = "blue", cex = 0.6)

##Scatter plot Discarded
#subset
subset_plant <- MetaCzechBayesian %>% filter(D_max > 0.0000, N_reads > 100,  grepl("Viridiplant",tax_path), grepl("VM", sample), grepl("\\bgenus\\b", tax_rank))
subset_animal <- MetaCzechBayesian %>% filter(D_max > 0.0000, N_reads > 100,  grepl("Metazoa",tax_path), grepl("VM", sample), grepl("\\bgenus\\b", tax_rank))

#creating subset for low abundance/high abundance taxa for better visualization
subset_plant1 <- subset_plant %>% filter(N_reads < 200) 
subset_plant2 <- subset_plant %>% filter(N_reads > 200) 

subset_plant_discarded1 <- subset_plant1 %>% filter(Bayesian_D_max < 0.05, Bayesian_z < 2 )
subset_plant_discarded2 <- subset_plant2 %>% filter(Bayesian_D_max < 0.05, Bayesian_z < 2 )
subset_animal_discarded <- subset_plant %>% filter(Bayesian_D_max < 0.05, Bayesian_z < 2 )

# Create Color Palette for animal taxa
subset_animal_discarded$tax_name <- as.factor(subset_animal_discarded$tax_name) 
length(levels(subset_animal_discarded$tax_name))  
magic.palette <- c(Bison ="#8B2323", Bos = "#EE7600", Bubalus = "#458B00", Canis = "#ABABAB", Capra= "#00008B", Dendroctonus = "#8B4726", Homo = "#79CDCD", Odocoileus = "#B2DFEE", Ovis = "#6CA6CD", Rangifer = "#8A2BE2", Rhagoletis = "#EEC591",Spodoptera = "#C71585", Sus = "#EEA2AD")
names(magic.palette) <- levels(subset_animal_discarded$tax_name) 

# Create shape selection for samples
subset_animal_discarded$sample <- as.factor(subset_animal_discarded$sample) 
length(levels(subset_animal_discarded$sample))
magic.shape <- c(6,15,16,17,18,4,7,8,3,10)
names(magic.shape) <- levels(subset_animal_discarded$sample) 

#Create Color Palette for plant taxa
library(RColorBrewer)
display.brewer.all()
library(RColorBrewer)
n <- 42
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

subset_plant_discarded$tax_name <- as.factor(subset_plant_discarded$tax_name) 
length(levels(subset_plant_discarded$tax_name))    
magic.palette2 <- col_vector
names(magic.palette2) <- levels(subset_plant_discarded$tax_name) 

# Create shape selection for samples
subset_plant_discarded$sample <- as.factor(subset_plant_discarded$sample) 
length(levels(subset_plant_discarded$sample))  
magic.shape2 <- c(6,15,16,17,18,4,7,8,3,10)
names(magic.shape2) <- levels(subset_plant_discarded$sample) 

#Scatter plot D_max vs Bayesian_z
p1 <- ggplot(subset_animal_discarded, aes(y=Bayesian_D_max, x=Bayesian_z)) + 
  geom_point(aes(col=tax_name, size=N_reads, shape=sample), show.legend = FALSE)+
  scale_size_continuous(range = c(3, 8))+
  scale_shape_manual(values = magic.shape,name = "Sample") + 
  scale_colour_manual (values = magic.palette,name = "Tax name") +
  scale_y_continuous (breaks = seq(0,0.150, by=0.025), limits = c(0,0.150))+
  guides(color = guide_legend(override.aes = list(size = 5)))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust =1)) + theme_minimal()
p1

p1a <- ggplot(subset_plant_discarded, aes(y=Bayesian_D_max, x=Bayesian_z)) + 
  geom_point(aes(col=tax_name, size=N_reads, shape=sample), show.legend = FALSE)+
  scale_size_continuous(range = c(3, 8))+
  scale_shape_manual(values = magic.shape2,name = "Sample") + 
  scale_colour_manual (values = magic.palette2,name = "Tax name") +
  scale_y_continuous (breaks = seq(0,0.125, by=0.025), limits = c(0,0.125))+
  guides(color = guide_legend(override.aes = list(size = 5)))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust =1)) + theme_minimal()
p1a

p1b <- ggplot(VM_Virid100_DISCARDED2below, aes(y=Bayesian_D_max, x=Bayesian_z)) + 
  geom_point(aes(col=tax_name, size=N_reads, shape=sample), show.legend = FALSE)+
  scale_size_continuous(range = c(3, 8))+
  scale_shape_manual(values = magic.shape2,name = "Sample") + 
  scale_colour_manual (values = magic.palette2,name = "Tax name") +
  scale_y_continuous (breaks = seq(0,0.125, by=0.025), limits = c(0,0.125))+
  guides(color = guide_legend(override.aes = list(size = 5)))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust =1)) + theme_minimal()
p1b

grid.arrange(p1, p1a, p1b, ncol = 3, nrow = 1)

####Supplementary Figures 3
##Read Length Distribution
# grepping read lengths from lca output of metaDMG and sorting (extracted reads that taxa), counting unique lengths on server, example:
grep 'nameGenus' nameafile.lca.txt | cut -f9 -d: | sort | uniq -c | sort -k2 -n > Genus.readLength.txt

#download data and input
Genus1 <- read_table("Sus.readLength.txt", col_names = FALSE)
Genus2 <- read_table("Ovis.readLength.txt", col_names = FALSE)
Genus3 <- read_table("Bos.readLength.txt", col_names = FALSE)
Genus4 <- read_table("Homo.readLength.txt", col_names = FALSE)
...

#Remove higher than bp 80
Genus1 <- Genus1[-c(52), ]

##One plot
d1 <- rbind(Genus1, Genus2, Genus3, Genus4,..)
colnames(d1) <- c("Number_of_reads", "Read_length")

#Multiple plots (Facets)
pdf("Readlengths.pdf", height = 6,width = 7.5)

p1 <- ggplot(d1, aes(Read_length, Number_of_reads)) + geom_point(stat="identity", col ="Dark green", alpha = 0.7)
p1 + xlab("Read length (basepairs Bp)")+ 
  ylab("Number of reads") + theme_bw()+ ggtitle("Readlength") + facet_wrap(~Taxa, scales="free_y") + 
+ theme(plot.title=element_text(size=12, hjust=0.5, vjust=0.5, face='bold'))
dev.off()

##Boxplot Depth vs Bayesian D-max
#Boxplot all reads
subset_all_genus <- MetaCzechBayesian %>% filter(D_max > 0.0000, N_reads > 100, grepl("VM", sample), grepl("\\bgenus\\b", tax_rank))

p1<- subset_all_genus %>%
  mutate(Depth = fct_relevel(Depth,
                             "236","222","196","185","151","134","116","107","64", "26", "18")) %>%
  ggplot(aes(x=D_max, y=Depth)) + geom_boxplot(aes(x=D_max, y=Depth, col=Kingdom)) #+ geom_jitter() + 
theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#Boxplot parsed reads
subset_parsed_all <- MetaCzechBayesian %>% filter(N_reads > 100, Bayesian_D_max > 0.05, Bayesian_phi > 100, Bayesian_z > 2, Bayesian_D_max_std < 0.10)

p2 <- subset_parsed_all %>%
  mutate(Depth2 = fct_relevel(Depth2,
                              "196","151","134","116","107","64")) %>%
  ggplot(aes(x=Bayesian_D_max, y=Depth2)) + geom_boxplot(aes(x=Bayesian_D_max, y=Depth2, col=Kingdom)) #+ geom_jitter() +
theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_x_continuous(breaks =seq(0,0.4, by = 0.10), limits=c(0,0.4))

grid.arrange(p1, p2,
             ncol = 2, nrow = 1, top="Title")

##Scatter Plots Bayesian D-max and Mean Length most abundant taxa
#metaDMG_make wide table - abundance indexes
sample2 <- subset_dataset$sample
dmax <- subset_dataset$Bayesian_D_max
taxid2 <- subset_dataset$tax_name
meanl <- subset_dataset$mean_L

subset_wide1 <- data.frame(taxid2, sample2, dmax)
subset_wide2 <- data.frame(taxid2, sample2, meanl)
subset_wide_filtered1 <- subset_wide1 %>% filter(N_reads > 1000)
subset_wide_filtered2 <- subset_wide2 %>% filter(N_reads > 1000)

#Scatter plot + standard deviation D-max
pd <- position_dodge(0.5)
p1 <- subset_wide_filtered1 %>%
  mutate(sample = fct_relevel(sample,
                              "VM-19","VM-17","VM-15","VM-14"))%>%
  ggplot(aes(x=Bayesian_D_max, y=sample, colour=tax_name)) + 
  geom_errorbar(aes(xmin=Bayesian_D_max-Bayesian_D_max_std, xmax=Bayesian_D_max+Bayesian_D_max_std, colour=tax_name, group=tax_name),position=pd) +
  geom_point(size = 2.2, position=pd)+
  scale_colour_manual(values=c("#006400", "#00CED1", "#EEC900", "#8B4500", "#66CD00"))+
  #scale_x_continuous(breaks = seq(0,0.15, by=0.05), limits = c(0,0.15))+
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Bayesian D-max by depth") + xlab ("Bayesian D-max")+ ylab ("Sample") # labs(colour="Tax")
p1

#Scatter plot + standard deviation Mean length
pd <- position_dodge(0.5)
p2 <- subset_wide_filtered2 %>%
  mutate(sample = fct_relevel(sample,
                              "VM-19","VM-17","VM-15","VM-14"))%>%
  ggplot(aes(x=mean_L, y=sample, colour=tax_name)) + 
  geom_errorbar(aes(xmin=mean_L-std_L, xmax=mean_L+std_L, colour=tax_name, group=tax_name),position=pd) +
  geom_point(size = 2.2, position=pd)+
  scale_colour_manual(values=c("#006400", "#00CED1", "#EEC900", "#8B4500", "#66CD00"))+
  scale_x_continuous(breaks = seq(30,90, by=20), limits = c(30,90))+
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Mean Length") + xlab ("Mean Length")+ ylab ("Sample") # labs(colour="Tax")
p2

grid.arrange(p1, p2,
             ncol = 2, nrow = 1, top="Title")

###Supplementary Figures 5
#MawSpecies_BARPLOT_meanlength
metaDMG_aMaw_all <- read_csv("/path_to_metaDMG_output/VM_metaDMG-aMaw.csv")
subset100_species <- metaDMG_aMaw_all %>% filter(N_reads>100, grepl("\\bspecies\\b", tax_rank))

# Magic-palette
subset100_species$sample <- as.factor(subset100_species$sample)
length(levels(subset100_species$sample))  

# Create Color Palette
magic.palette <- c("#9A32CD", "#00CD00", "#050505", "#CD3333", "#008B00", "#1C86EE", "#00EEEE", "#EEB422", "#A3A3A3", "#00008B", "#EE6A50")    # defining 7 colours
names(magic.palette) <- levels(subset100_species$sample)                                                    

##Microbial Species Bayesian D-max vs Bayesian_z
# counting number of unique samples for number of different symbols
f <- length(unique(subset100_species$sample))
#plotting  D_max vs lambda_LR
p1 <- ggplot(subset100_species, aes(y=Bayesian_D_max, x=Bayesian_z)) + geom_point(aes(size=N_reads,col=sample)) +
  scale_colour_manual(values = magic.palette,name = "Sample") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust =1)) +
  scale_x_continuous(breaks =seq(-2,10, by=2), limits=c(-2,10))+ 
  ggtitle("Damage vs BayesianZ") + xlab("Bayesian Z") + ylab("D_max")
p1 + theme_minimal()
#+ scale_x_continuous(breaks =seq(0,150, by=10), limits=c(0,150)) + scale_y_continuous(breaks =seq(0,0.6, by=0.05))

##Microbial Species Barplot Dmax
Depth2 <- read.csv ("/path_to_metadata/VM_Sediment_context.csv", sep=";") 

subset100_species$new <- Depth2$depth_cm[match(subset100_species$sample, Depth2$sample_ID)]
names(subset100_species)[names(subset100_species) == 'new'] <- 'Depth'

write.csv(subset100_species, file = "subset100_species.csv")
subset100_species$Depth <- as.factor(subset100_species$Depth)

#Mean_Length
p2 <- subset100_species %>%
  mutate(Depth = fct_relevel(Depth,
                             "236","222","196","185","151","134","116","107","64", "26", "18")) %>%
  ggplot(aes(x=mean_L, y=Depth)) + geom_boxplot(aes(x=mean_L, y=Depth, colour=sample))+ scale_colour_manual(values = magic.palette,name = "Sample") + 
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#+ scale_x_continuous(breaks =seq(0,0.4, by = 0.10), limits=c(0,0.4))
p2 

#D-max
p3 <- subset100_species %>%
  mutate(Depth = fct_relevel(Depth,
                             "236","222","196","185","151","134","116","107","64", "26", "18")) %>%
  ggplot(aes(x=Bayesian_D_max, y=Depth)) + geom_boxplot(aes(x=Bayesian_D_max, y=Depth, colour=sample))+ scale_colour_manual(values = magic.palette,name = "Sample") + 
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p3

grid.arrange(p2, p3, 
             ncol = 2, nrow = 1)
