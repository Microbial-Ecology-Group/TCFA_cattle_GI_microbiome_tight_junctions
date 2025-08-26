# Start Here --------------------------------------------------------------
#Project: TCFA IHC
#Author: JDY Start: 17Oct24 Finished:
#Goal: Compare computer generated data to pathologist generated scores
#Note: The data sets used here have been already analyzed in different scripts, the focus here is the differences between methods.


# Library -----------------------------------------------------------------

setwd("C:/Users/danie/OneDrive - West Texas A and M University/Research/TCFA Challenge/IHC Scoring/Comparison of Methods/")

library(readxl)
library(ggplot2)
library(lme4)
library(car)
library(emmeans)
library(dplyr)
library(tidyr)
library(janitor)
library(cowplot)


# CLDN1 -------------------------------------------------------------------
#read in computer generated data
cldn1_cgd <- read_excel("Cleaned_measurements.xlsx", sheet = "CLDN1")

#read in meta data
cldn1_metadata <- read_excel("TCFA_IHC_Area_DAta.xlsx", sheet = "CLDN1")

#merge data
cldn1_merged <- merge(cldn1_metadata,cldn1_cgd)

#summarise all annotations from a single sample type to get an average score for each tissue
cldn1_data <- cldn1_merged %>% group_by(Sample_name,Tissue) %>% 
  summarise(H_negative= mean(`H-score: Negative %`),
            H_pixelwise = mean(`Pixelwise H-score`))


#read in the pathologist scores and clean the names
cldn1_dat <- read_excel("TCFA_IHC_PathologistHScores_Cleaned_for_R.xlsx", sheet = "CLDN1")
cldn1_dat <- clean_names(cldn1_dat)

#average the scores across the scorer
path_data <- cldn1_dat %>% group_by(sample_id, organ) %>% 
  summarise(epi_h= mean(epithelium_h_score, na.rm = T),
            lp_h= mean(lp_h_score, na.rm = T)) %>% 
  rename(Sample_name = sample_id)

#make names match
path_data$Sample_name <- gsub("Large Intestine-", "LI-", path_data$Sample_name)
path_data$Sample_name <- gsub("Small Intestine-", "SI-", path_data$Sample_name)
path_data$Sample_name <- gsub("Rumen-", "RU-", path_data$Sample_name)

# merge data sets
cldn1_final <- merge(path_data,cldn1_data)

##epi
#make calulations
cldn1_means <- (cldn1_final$epi_h + cldn1_final$H_pixelwise)/2
cldn1_diff <- cldn1_final$epi_h - cldn1_final$H_pixelwise
cldn1_mean_diff <- mean(cldn1_diff)
cldn1_sd <- sd(cldn1_diff)


cldn1_epi <- ggplot(cldn1_final, aes(x = cldn1_means, y = cldn1_diff))+
  geom_point()+
  geom_hline(yintercept = cldn1_mean_diff, color = "blue")+
  geom_hline(yintercept = cldn1_mean_diff + 1.96*cldn1_sd, color = "red", linetype = "dashed")+
  geom_hline(yintercept = cldn1_mean_diff - 1.96*cldn1_sd, color = "red", linetype = "dashed")+
  labs(title = "CLDN1: Epi H vs Pixelwise H", x = "Mean of Measurements", y = "Difference Between Measurements") +
  annotate("text", x = 195, y = 50, label = "Mean bias = 44.2")+
  theme_bw()+theme(legend.position = "none",
                   plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
                   panel.border = element_rect(colour = "black", size = 1.7),
                   axis.ticks = element_line(size = 1, colour = "black"),
                   plot.title = element_text(size = 30),
                   axis.title.y = element_text(size = 26, vjust = 2.5),
                   axis.text.x = element_text(size = 24, colour = "black"),
                   axis.text.y = element_text(size = 24, colour = "black"),
                   axis.title.x = element_text(size = 26, vjust = 2.5),
                   panel.grid.major.x = element_blank(),
                   legend.text = element_text(size = 40, hjust = 0.5),
                   legend.title = element_text(size = 45, hjust = 0.5),
                   legend.key.height = unit(2,"cm"))
cldn1_epi

#lp
#make calulations
cldn1_means2 <- (cldn1_final$lp_h + cldn1_final$H_pixelwise)/2
cldn1_diff2 <- cldn1_final$lp_h - cldn1_final$H_pixelwise
cldn1_mean_diff2 <- mean(cldn1_diff2)
cldn1_sd2 <- sd(cldn1_diff2)


cldn1_lp <- ggplot(cldn1_final, aes(x = cldn1_means2, y = cldn1_diff2))+
  geom_point()+
  geom_hline(yintercept = cldn1_mean_diff2, color = "blue")+
  geom_hline(yintercept = cldn1_mean_diff2 + 1.96*cldn1_sd2, color = "red", linetype = "dashed")+
  geom_hline(yintercept = cldn1_mean_diff2 - 1.96*cldn1_sd2, color = "red", linetype = "dashed")+
  labs(title = "CLDN1: LP H vs Pixelwise H", x = "Mean of Measurements", y = "Difference Between Measurements") +
  annotate("text", x = 132, y = -35, label = "Mean bias = -43.16")+
  theme_bw()+theme(legend.position = "none",
                   plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
                   panel.border = element_rect(colour = "black", size = 1.7),
                   axis.ticks = element_line(size = 1, colour = "black"),
                   plot.title = element_text(size = 30),
                   axis.title.y = element_text(size = 26, vjust = 2.5),
                   axis.text.x = element_text(size = 24, colour = "black"),
                   axis.text.y = element_text(size = 24, colour = "black"),
                   axis.title.x = element_text(size = 26, vjust = 2.5),
                   panel.grid.major.x = element_blank(),
                   legend.text = element_text(size = 40, hjust = 0.5),
                   legend.title = element_text(size = 45, hjust = 0.5),
                   legend.key.height = unit(2,"cm"))
cldn1_lp


#create an average of epi and lp score to more acurately match the pixelwise
cldn1_final$h_avg <- (cldn1_final$epi_h+cldn1_final$lp_h)/2

#make calulations
cldn1_means1 <- (cldn1_final$h_avg + cldn1_final$H_pixelwise)/2
cldn1_diff1 <- cldn1_final$h_avg - cldn1_final$H_pixelwise
cldn1_mean_diff1 <- mean(cldn1_diff1)
cldn1_sd1 <- sd(cldn1_diff1)

cldn1_avg <- ggplot(cldn1_final, aes(x = cldn1_means1, y = cldn1_diff1))+
  geom_point()+
  geom_hline(yintercept = cldn1_mean_diff1, color = "blue")+
  geom_hline(yintercept = cldn1_mean_diff1 + 1.96*cldn1_sd1, color = "red", linetype = "dashed")+
  geom_hline(yintercept = cldn1_mean_diff1 - 1.96*cldn1_sd1, color = "red", linetype = "dashed")+
  labs(title = "CLDN1: AVG H vs Pixelwise H", x = "Mean of Measurements", y = "Difference Between Measurements") +
  annotate("text", x = 156, y = 5, label = "Mean bias = 0.511")+
  theme_bw()+theme(legend.position = "none",
                   plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
                   panel.border = element_rect(colour = "black", size = 1.7),
                   axis.ticks = element_line(size = 1, colour = "black"),
                   plot.title = element_text(size = 30),
                   axis.title.y = element_text(size = 26, vjust = 2.5),
                   axis.text.x = element_text(size = 24, colour = "black"),
                   axis.text.y = element_text(size = 24, colour = "black"),
                   axis.title.x = element_text(size = 26, vjust = 2.5),
                   panel.grid.major.x = element_blank(),
                   legend.text = element_text(size = 40, hjust = 0.5),
                   legend.title = element_text(size = 45, hjust = 0.5),
                   legend.key.height = unit(2,"cm"))
cldn1_avg


# CLDN2 -------------------------------------------------------------------
#read in computer generated data
cldn2_cgd <- read_excel("Cleaned_measurements.xlsx", sheet = "CLDN2")

#read in meta data
cldn2_metadata <- read_excel("TCFA_IHC_Area_DAta.xlsx", sheet = "CLDN2")


#merge data
cldn2_merged <- merge(cldn2_metadata,cldn2_cgd)

#summarise all annotations from a single sample type to get an average score for each tissue
cldn2_data <- cldn2_merged %>% group_by(Sample_name,Tissue) %>% 
  summarise(H_negative= mean(`H-score: Negative %`),
            H_pixelwise = mean(`Pixelwise H-score`))

#read in pathologist data
cldn2_dat <- read_excel("TCFA_IHC_PathologistHScores_Cleaned_for_R.xlsx", sheet = "CLDN2")
cldn2_dat <- clean_names(cldn2_dat)

#average the scores across the scorer
cldn2_path_data <- cldn2_dat %>% group_by(sample_id, organ) %>% 
  summarise(epi_h= mean(epithelium_h_score, na.rm = T),
            lp_h= mean(lp_h_score, na.rm = T)) %>% 
  rename(Sample_name = sample_id)

#make names match
cldn2_path_data$Sample_name <- gsub("Large Intestine-", "LI-", cldn2_path_data$Sample_name)
cldn2_path_data$Sample_name <- gsub("Small Intestine-", "SI-", cldn2_path_data$Sample_name)
cldn2_path_data$Sample_name <- gsub("Rumen-", "RU-", cldn2_path_data$Sample_name)

# merge data sets
cldn2_final <- merge(cldn2_path_data,cldn2_data)

#make calulations epi
cldn2_means <- (cldn2_final$epi_h + cldn2_final$H_pixelwise)/2
cldn2_diff <- cldn2_final$epi_h - cldn2_final$H_pixelwise
cldn2_mean_diff <- mean(cldn2_diff)
cldn2_sd <- sd(cldn2_diff)

cldn2_epi <- ggplot(cldn2_final, aes(x = cldn2_means, y = cldn2_diff))+
  geom_point()+
  geom_hline(yintercept = cldn2_mean_diff, color = "blue")+
  geom_hline(yintercept = cldn2_mean_diff + 1.96*cldn2_sd, color = "red", linetype = "dashed")+
  geom_hline(yintercept = cldn2_mean_diff - 1.96*cldn2_sd, color = "red", linetype = "dashed")+
  labs(title = "CLDN2: Epi H vs Pixelwise H", x = "Mean of Measurements", y = "Difference Between Measurements") +
  annotate("text", x = 145, y = 15, label = "Mean bias = 7.85")+
  theme_bw()+theme(legend.position = "none",
                   plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
                   panel.border = element_rect(colour = "black", size = 1.7),
                   axis.ticks = element_line(size = 1, colour = "black"),
                   plot.title = element_text(size = 30),
                   axis.title.y = element_text(size = 26, vjust = 2.5),
                   axis.text.x = element_text(size = 24, colour = "black"),
                   axis.text.y = element_text(size = 24, colour = "black"),
                   axis.title.x = element_text(size = 26, vjust = 2.5),
                   panel.grid.major.x = element_blank(),
                   legend.text = element_text(size = 40, hjust = 0.5),
                   legend.title = element_text(size = 45, hjust = 0.5),
                   legend.key.height = unit(2,"cm"))
cldn2_epi

#make calulations lp
cldn2_means2 <- (cldn2_final$lp_h + cldn2_final$H_pixelwise)/2
cldn2_diff2 <- cldn2_final$lp_h - cldn2_final$H_pixelwise
cldn2_mean_diff2 <- mean(cldn2_diff2)
cldn2_sd2 <- sd(cldn2_diff2)

cldn2_lp <- ggplot(cldn2_final, aes(x = cldn2_means2, y = cldn2_diff2))+
  geom_point()+
  geom_hline(yintercept = cldn2_mean_diff2, color = "blue")+
  geom_hline(yintercept = cldn2_mean_diff2 + 1.96*cldn2_sd2, color = "red", linetype = "dashed")+
  geom_hline(yintercept = cldn2_mean_diff2 - 1.96*cldn2_sd2, color = "red", linetype = "dashed")+
  labs(title = "CLDN2: LP H vs Pixelwise H", x = "Mean of Measurements", y = "Difference Between Measurements") +
  annotate("text", x = 82, y = -73, label = "Mean bias = -75.41")+
  theme_bw()+theme(legend.position = "none",
                   plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
                   panel.border = element_rect(colour = "black", size = 1.7),
                   axis.ticks = element_line(size = 1, colour = "black"),
                   plot.title = element_text(size = 30),
                   axis.title.y = element_text(size = 26, vjust = 2.5),
                   axis.text.x = element_text(size = 24, colour = "black"),
                   axis.text.y = element_text(size = 24, colour = "black"),
                   axis.title.x = element_text(size = 26, vjust = 2.5),
                   panel.grid.major.x = element_blank(),
                   legend.text = element_text(size = 40, hjust = 0.5),
                   legend.title = element_text(size = 45, hjust = 0.5),
                   legend.key.height = unit(2,"cm"))
cldn2_lp



#create an average of epi and lp score to more acurately match the pixelwise
cldn2_final$h_avg <- (cldn2_final$epi_h+cldn2_final$lp_h)/2

#make calulations 
cldn2_means1 <- (cldn2_final$h_avg + cldn2_final$H_pixelwise)/2
cldn2_diff1 <- cldn2_final$h_avg - cldn2_final$H_pixelwise
cldn2_mean_diff1 <- mean(cldn2_diff1)
cldn2_sd1 <- sd(cldn2_diff1)

cldn2_avg <- ggplot(cldn2_final, aes(x = cldn2_means1, y = cldn2_diff1))+
  geom_point()+
  geom_hline(yintercept = cldn2_mean_diff1, color = "blue")+
  geom_hline(yintercept = cldn2_mean_diff1 + 1.96*cldn2_sd1, color = "red", linetype = "dashed")+
  geom_hline(yintercept = cldn2_mean_diff1 - 1.96*cldn2_sd1, color = "red", linetype = "dashed")+
  labs(title = "CLDN2: AVG H vs Pixelwise H", x = "Mean of Measurements", y = "Difference Between Measurements") +
  annotate("text", x = 110, y = -40, label = "Mean bias = -33.78")+
  theme_bw()+theme(legend.position = "none",
                   plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
                   panel.border = element_rect(colour = "black", size = 1.7),
                   axis.ticks = element_line(size = 1, colour = "black"),
                   plot.title = element_text(size = 30),
                   axis.title.y = element_text(size = 26, vjust = 2.5),
                   axis.text.x = element_text(size = 24, colour = "black"),
                   axis.text.y = element_text(size = 24, colour = "black"),
                   axis.title.x = element_text(size = 26, vjust = 2.5),
                   panel.grid.major.x = element_blank(),
                   legend.text = element_text(size = 40, hjust = 0.5),
                   legend.title = element_text(size = 45, hjust = 0.5),
                   legend.key.height = unit(2,"cm"))
cldn2_avg


# ECAD --------------------------------------------------------------------
#read in computer generated data
ecad_cgd <- read_excel("Cleaned_measurements.xlsx", sheet = "ECAD")

#read in meta data
ecad_metadata <- read_excel("TCFA_IHC_Area_DAta.xlsx", sheet = "ECAD")


#merge data
ecad_merged <- merge(ecad_metadata, ecad_cgd)

#summarise all annotations from a single sample type to get an average score for each tissue
ecad_data <- ecad_merged %>% group_by(Sample_Name,Tissue) %>% 
  summarise(H_negative= mean(`H-score: Negative %`),
            H_pixelwise = mean(`Pixelwise H-score`))

#read data and clean names
ecad_dat <- read_excel("TCFA_IHC_PathologistHScores_Cleaned_for_R.xlsx", sheet = "E-Cadherin")
ecad_dat <- clean_names(ecad_dat)

#average the scores across the scorer
ecad_path_data <- ecad_dat %>% group_by(sample_id, organ) %>% 
  summarise(epi_h= mean(epithelium_h_score, na.rm = T),
            lp_h= mean(lp_h_score, na.rm = T)) %>% 
  rename(Sample_Name = sample_id)

#make names match
ecad_path_data$Sample_Name <- gsub("Large Intestine-", "LI-", ecad_path_data$Sample_Name)
ecad_path_data$Sample_Name <- gsub("Small Intestine-", "SI-", ecad_path_data$Sample_Name)
ecad_path_data$Sample_Name <- gsub("Rumen-", "RU-", ecad_path_data$Sample_Name)

# merge data sets
ecad_final <- merge(ecad_path_data,ecad_data)

#make calculations epi
ecad_means <- (ecad_final$epi_h + ecad_final$H_pixelwise)/2
ecad_diff <- ecad_final$epi_h - ecad_final$H_pixelwise
ecad_mean_diff <- mean(ecad_diff)
ecad_sd <- sd(ecad_diff)

ecad_epi <- ggplot(ecad_final, aes(x = ecad_means, y = ecad_diff))+
  geom_point()+
  geom_hline(yintercept = ecad_mean_diff, color = "blue")+
  geom_hline(yintercept = ecad_mean_diff + 1.96*ecad_sd, color = "red", linetype = "dashed")+
  geom_hline(yintercept = ecad_mean_diff - 1.96*ecad_sd, color = "red", linetype = "dashed")+
  labs(title = "E-CAD: Epi H vs Pixelwise H", x = "Mean of Measurements", y = "Difference Between Measurements") +
  annotate("text", x = 260, y = 79, label = "Mean bias = 73.66")+
  theme_bw()+theme(legend.position = "none",
                   plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
                   panel.border = element_rect(colour = "black", size = 1.7),
                   axis.ticks = element_line(size = 1, colour = "black"),
                   plot.title = element_text(size = 30),
                   axis.title.y = element_text(size = 26, vjust = 2.5),
                   axis.text.x = element_text(size = 24, colour = "black"),
                   axis.text.y = element_text(size = 24, colour = "black"),
                   axis.title.x = element_text(size = 26, vjust = 2.5),
                   panel.grid.major.x = element_blank(),
                   legend.text = element_text(size = 40, hjust = 0.5),
                   legend.title = element_text(size = 45, hjust = 0.5),
                   legend.key.height = unit(2,"cm"))
ecad_epi

#make calculations lp
ecad_means2 <- (ecad_final$lp_h + ecad_final$H_pixelwise)/2
ecad_diff2 <- ecad_final$lp_h - ecad_final$H_pixelwise
ecad_mean_diff2 <- mean(ecad_diff2)
ecad_sd2 <- sd(ecad_diff2)

ecad_lp <- ggplot(ecad_final, aes(x = ecad_means2, y = ecad_diff2))+
  geom_point()+
  geom_hline(yintercept = ecad_mean_diff2, color = "blue")+
  geom_hline(yintercept = ecad_mean_diff2 + 1.96*ecad_sd2, color = "red", linetype = "dashed")+
  geom_hline(yintercept = ecad_mean_diff2 - 1.96*ecad_sd2, color = "red", linetype = "dashed")+
  labs(title = "E-CAD: LP H vs Pixelwise H", x = "Mean of Measurements", y = "Difference Between Measurements") +
  annotate("text", x = 140, y = -155, label = "Mean bias = -166.92")+
  theme_bw()+theme(legend.position = "none",
                   plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
                   panel.border = element_rect(colour = "black", size = 1.7),
                   axis.ticks = element_line(size = 1, colour = "black"),
                   plot.title = element_text(size = 30),
                   axis.title.y = element_text(size = 26, vjust = 2.5),
                   axis.text.x = element_text(size = 24, colour = "black"),
                   axis.text.y = element_text(size = 24, colour = "black"),
                   axis.title.x = element_text(size = 26, vjust = 2.5),
                   panel.grid.major.x = element_blank(),
                   legend.text = element_text(size = 40, hjust = 0.5),
                   legend.title = element_text(size = 45, hjust = 0.5),
                   legend.key.height = unit(2,"cm"))
ecad_lp


#create an average of epi and lp score to more acurately match the pixelwise
ecad_final$h_avg <- (ecad_final$epi_h+ecad_final$lp_h)/2

#make calulations
ecad_means1 <- (ecad_final$h_avg + ecad_final$H_pixelwise)/2
ecad_diff1 <- ecad_final$h_avg - ecad_final$H_pixelwise
ecad_mean_diff1 <- mean(ecad_diff1)
ecad_sd1 <- sd(ecad_diff1)

ecad_avg <- ggplot(ecad_final, aes(x = ecad_means1, y = ecad_diff1))+
  geom_point()+
  geom_hline(yintercept = ecad_mean_diff1, color = "blue")+
  geom_hline(yintercept = ecad_mean_diff1 + 1.96*ecad_sd1, color = "red", linetype = "dashed")+
  geom_hline(yintercept = ecad_mean_diff1 - 1.96*ecad_sd1, color = "red", linetype = "dashed")+
  labs(title = "E-CAD: AVG H vs Pixelwise H", x = "Mean of Measurements", y = "Difference Between Measurements") +
  annotate("text", x = 198, y = -40, label = "Mean bias = -46.63")+
  theme_bw()+theme(legend.position = "none",
                   plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
                   panel.border = element_rect(colour = "black", size = 1.7),
                   axis.ticks = element_line(size = 1, colour = "black"),
                   plot.title = element_text(size = 30),
                   axis.title.y = element_text(size = 26, vjust = 2.5),
                   axis.text.x = element_text(size = 24, colour = "black"),
                   axis.text.y = element_text(size = 24, colour = "black"),
                   axis.title.x = element_text(size = 26, vjust = 2.5),
                   panel.grid.major.x = element_blank(),
                   legend.text = element_text(size = 40, hjust = 0.5),
                   legend.title = element_text(size = 45, hjust = 0.5),
                   legend.key.height = unit(2,"cm"))
ecad_avg


# OCLDN -------------------------------------------------------------------
#read in computer generated data
ocldn_cgd <- read_excel("Cleaned_measurements.xlsx", sheet = "OCLDN")

#read in meta data
ocldn_metadata <- read_excel("TCFA_IHC_Area_DAta.xlsx", sheet = "OCLDN")

#merge data
ocldn_merged <- merge(ocldn_metadata,ocldn_cgd)

#summarise all annotations from a single sample type to get an average score for each tissue
ocldn_data <- ocldn_merged %>% group_by(Sample_name,Tissue) %>% 
  summarise(H_negative= mean(`H-score: Negative %`),
            H_pixelwise = mean(`Pixelwise H-score`))

#read data and clean names
ocldn_dat <- read_excel("TCFA_IHC_PathologistHScores_Cleaned_for_R.xlsx", sheet = "OCLDN")
ocldn_dat <- clean_names(ocldn_dat)

#average the scores across the scorer
ocldn_path_data <- ocldn_dat %>% group_by(sample_id, organ) %>% 
  summarise(epi_h= mean(epithelium_h_score, na.rm = T),
            lp_h= mean(lp_h_score, na.rm = T)) %>% 
  rename(Sample_name = sample_id)

#make names match
ocldn_path_data$Sample_name <- gsub("Large Intestine-", "LI-", ocldn_path_data$Sample_name)
ocldn_path_data$Sample_name <- gsub("Small Intestine-", "SI-", ocldn_path_data$Sample_name)
ocldn_path_data$Sample_name <- gsub("Rumen-", "RU-", ocldn_path_data$Sample_name)

# merge data sets
ocldn_final <- merge(ocldn_path_data,ocldn_data)

#make calculations epi
ocldn_means <- (ocldn_final$epi_h + ocldn_final$H_pixelwise)/2
ocldn_diff <- ocldn_final$epi_h - ocldn_final$H_pixelwise
ocldn_mean_diff <- mean(ocldn_diff)
ocldn_sd <- sd(ocldn_diff)

ocldn_epi <- ggplot(ocldn_final, aes(x = ocldn_means, y = ocldn_diff))+
  geom_point()+
  geom_hline(yintercept = ocldn_mean_diff, color = "blue")+
  geom_hline(yintercept = ocldn_mean_diff + 1.96*ocldn_sd, color = "red", linetype = "dashed")+
  geom_hline(yintercept = ocldn_mean_diff - 1.96*ocldn_sd, color = "red", linetype = "dashed")+
  labs(title = "OCLDN: Epi H vs Pixelwise H", x = "Mean of Measurements", y = "Difference Between Measurements") +
  annotate("text", x = 194, y = 115, label = "Mean bias = 106.47")+
  theme_bw()+theme(legend.position = "none",
                   plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
                   panel.border = element_rect(colour = "black", size = 1.7),
                   axis.ticks = element_line(size = 1, colour = "black"),
                   plot.title = element_text(size = 30),
                   axis.title.y = element_text(size = 26, vjust = 2.5),
                   axis.text.x = element_text(size = 24, colour = "black"),
                   axis.text.y = element_text(size = 24, colour = "black"),
                   axis.title.x = element_text(size = 26, vjust = 2.5),
                   panel.grid.major.x = element_blank(),
                   legend.text = element_text(size = 40, hjust = 0.5),
                   legend.title = element_text(size = 45, hjust = 0.5),
                   legend.key.height = unit(2,"cm"))
ocldn_epi

#make calculations lp
ocldn_means2 <- (ocldn_final$lp_h + ocldn_final$H_pixelwise)/2
ocldn_diff2 <- ocldn_final$lp_h - ocldn_final$H_pixelwise
ocldn_mean_diff2 <- mean(ocldn_diff2)
ocldn_sd2 <- sd(ocldn_diff2)

ocldn_lp <- ggplot(ocldn_final, aes(x = ocldn_means2, y = ocldn_diff2))+
  geom_point()+
  geom_hline(yintercept = ocldn_mean_diff2, color = "blue")+
  geom_hline(yintercept = ocldn_mean_diff2 + 1.96*ocldn_sd2, color = "red", linetype = "dashed")+
  geom_hline(yintercept = ocldn_mean_diff2 - 1.96*ocldn_sd2, color = "red", linetype = "dashed")+
  labs(title = "OCLDN: LP H vs Pixelwise H", x = "Mean of Measurements", y = "Difference Between Measurements") +
  annotate("text", x = 79, y = -98, label = "Mean bias = -101.93")+
  theme_bw()+theme(legend.position = "none",
                   plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
                   panel.border = element_rect(colour = "black", size = 1.7),
                   axis.ticks = element_line(size = 1, colour = "black"),
                   plot.title = element_text(size = 30),
                   axis.title.y = element_text(size = 26, vjust = 2.5),
                   axis.text.x = element_text(size = 24, colour = "black"),
                   axis.text.y = element_text(size = 24, colour = "black"),
                   axis.title.x = element_text(size = 26, vjust = 2.5),
                   panel.grid.major.x = element_blank(),
                   legend.text = element_text(size = 40, hjust = 0.5),
                   legend.title = element_text(size = 45, hjust = 0.5),
                   legend.key.height = unit(2,"cm"))
ocldn_lp



#create an average of epi and lp score to more acurately match the pixelwise
ocldn_final$h_avg <- (ocldn_final$epi_h+ocldn_final$lp_h)/2

#make calulations
ocldn_means1 <- (ocldn_final$h_avg + ocldn_final$H_pixelwise)/2
ocldn_diff1 <- ocldn_final$h_avg - ocldn_final$H_pixelwise
ocldn_mean_diff1 <- mean(ocldn_diff1)
ocldn_sd1 <- sd(ocldn_diff1)

ocldn_avg <- ggplot(ocldn_final, aes(x = ocldn_means1, y = ocldn_diff1))+
  geom_point()+
  geom_hline(yintercept = ocldn_mean_diff1, color = "blue")+
  geom_hline(yintercept = ocldn_mean_diff1 + 1.96*ocldn_sd1, color = "red", linetype = "dashed")+
  geom_hline(yintercept = ocldn_mean_diff1 - 1.96*ocldn_sd1, color = "red", linetype = "dashed")+
  labs(title = "OCLDN: AVG H vs Pixelwise H", x = "Mean of Measurements", y = "Difference Between Measurements") +
  annotate("text", x = 135, y = 7, label = "Mean bias = 2.27")+
  theme_bw()+theme(legend.position = "none",
                   plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
                   panel.border = element_rect(colour = "black", size = 1.7),
                   axis.ticks = element_line(size = 1, colour = "black"),
                   plot.title = element_text(size = 30),
                   axis.title.y = element_text(size = 26, vjust = 2.5),
                   axis.text.x = element_text(size = 24, colour = "black"),
                   axis.text.y = element_text(size = 24, colour = "black"),
                   axis.title.x = element_text(size = 26, vjust = 2.5),
                   panel.grid.major.x = element_blank(),
                   legend.text = element_text(size = 40, hjust = 0.5),
                   legend.title = element_text(size = 45, hjust = 0.5),
                   legend.key.height = unit(2,"cm"))
ocldn_avg


# MAC-387 -----------------------------------------------------------------
#read in computer generated data
mac_cgd <- read_excel("Cleaned_measurements.xlsx", sheet = "MAC387")

#read in meta data
mac_metadata <- read_excel("TCFA_IHC_Area_DAta.xlsx", sheet = "MAC387")


#merge data
mac_merged <- merge(mac_metadata,mac_cgd)

#summarise all annotations from a single sample type to get an average score for each tissue
mac_data <- mac_merged %>% group_by(Sample_name,Tissue) %>% 
  summarise(H_negative= mean(`H-score: Negative %`),
            H_pixelwise = mean(`Pixelwise H-score`))

#read data and clean names
mac_dat <- read_excel("TCFA_IHC_PathologistHScores_Cleaned_for_R.xlsx", sheet = "MAC387")
mac_dat <- clean_names(mac_dat)

#average the scores across the scorer
mac_path_data <- mac_dat %>% group_by(sample_id, organ) %>% 
  summarise(epi_h= mean(epithelium_h_score, na.rm = T),
            lp_h= mean(lp_h_score, na.rm = T)) %>% 
  rename(Sample_name = sample_id)

#make names match
mac_path_data$Sample_name <- gsub("Large Intestine-", "LI-", mac_path_data$Sample_name)
mac_path_data$Sample_name <- gsub("Small Intestine-", "SI-", mac_path_data$Sample_name)
mac_path_data$Sample_name <- gsub("Rumen-", "RU-", mac_path_data$Sample_name)

# merge data sets
mac_final <- merge(mac_path_data,mac_data)

#make calculations epi
mac_means <- (mac_final$epi_h + mac_final$H_pixelwise)/2
mac_diff <- mac_final$epi_h - mac_final$H_pixelwise
mac_mean_diff <- mean(mac_diff)
mac_sd <- sd(mac_diff)

mac_epi <- ggplot(mac_final, aes(x = mac_means, y = mac_diff))+
  geom_point()+
  geom_hline(yintercept = mac_mean_diff, color = "blue")+
  geom_hline(yintercept = mac_mean_diff + 1.96*mac_sd, color = "red", linetype = "dashed")+
  geom_hline(yintercept = mac_mean_diff - 1.96*mac_sd, color = "red", linetype = "dashed")+
  labs(title = "MAC-387: Epi H vs Pixelwise H", x = "Mean of Measurements", y = "Difference Between Measurements") +
  annotate("text", x = 205, y = -25, label = "Mean bias = -19.34")+
  theme_bw()+theme(legend.position = "none",
                   plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
                   panel.border = element_rect(colour = "black", size = 1.7),
                   axis.ticks = element_line(size = 1, colour = "black"),
                   plot.title = element_text(size = 30),
                   axis.title.y = element_text(size = 26, vjust = 2.5),
                   axis.text.x = element_text(size = 24, colour = "black"),
                   axis.text.y = element_text(size = 24, colour = "black"),
                   axis.title.x = element_text(size = 26, vjust = 2.5),
                   panel.grid.major.x = element_blank(),
                   legend.text = element_text(size = 40, hjust = 0.5),
                   legend.title = element_text(size = 45, hjust = 0.5),
                   legend.key.height = unit(2,"cm"))
mac_epi

#make calculations lp
mac_means2 <- (mac_final$lp_h + mac_final$H_pixelwise)/2
mac_diff2 <- mac_final$lp_h - mac_final$H_pixelwise
mac_mean_diff2 <- mean(mac_diff2)
mac_sd2 <- sd(mac_diff2)

mac_lp <- ggplot(mac_final, aes(x = mac_means2, y = mac_diff2))+
  geom_point()+
  geom_hline(yintercept = mac_mean_diff2, color = "blue")+
  geom_hline(yintercept = mac_mean_diff2 + 1.96*mac_sd2, color = "red", linetype = "dashed")+
  geom_hline(yintercept = mac_mean_diff2 - 1.96*mac_sd2, color = "red", linetype = "dashed")+
  labs(title = "MAC-387: LP H vs Pixelwise H", x = "Mean of Measurements", y = "Difference Between Measurements") +
  annotate("text", x = 123, y = -40, label = "Mean bias = -57.30")+
  theme_bw()+theme(legend.position = "none",
                   plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
                   panel.border = element_rect(colour = "black", size = 1.7),
                   axis.ticks = element_line(size = 1, colour = "black"),
                   plot.title = element_text(size = 30),
                   axis.title.y = element_text(size = 26, vjust = 2.5),
                   axis.text.x = element_text(size = 24, colour = "black"),
                   axis.text.y = element_text(size = 24, colour = "black"),
                   axis.title.x = element_text(size = 26, vjust = 2.5),
                   panel.grid.major.x = element_blank(),
                   legend.text = element_text(size = 40, hjust = 0.5),
                   legend.title = element_text(size = 45, hjust = 0.5),
                   legend.key.height = unit(2,"cm"))
mac_lp

#create an average of epi and lp score to more acurately match the pixelwise
mac_final$h_avg <- (mac_final$epi_h+mac_final$lp_h)/2

#make calulations
mac_means1 <- (mac_final$h_avg + mac_final$H_pixelwise)/2
mac_diff1 <- mac_final$h_avg - mac_final$H_pixelwise
mac_mean_diff1 <- mean(mac_diff1)
mac_sd1 <- sd(mac_diff1)

mac_avg <- ggplot(mac_final, aes(x = mac_means1, y = mac_diff1))+
  geom_point()+
  geom_hline(yintercept = mac_mean_diff1, color = "blue")+
  geom_hline(yintercept = mac_mean_diff1 + 1.96*mac_sd1, color = "red", linetype = "dashed")+
  geom_hline(yintercept = mac_mean_diff1 - 1.96*mac_sd1, color = "red", linetype = "dashed")+
  labs(title = "MAC-387: AVG H vs Pixelwise H", x = "Mean of Measurements", y = "Difference Between Measurements") +
  annotate("text", x = 165, y = -25, label = "Mean bias = -38.32")+
  theme_bw()+theme(legend.position = "none",
                   plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
                   panel.border = element_rect(colour = "black", size = 1.7),
                   axis.ticks = element_line(size = 1, colour = "black"),
                   plot.title = element_text(size = 30),
                   axis.title.y = element_text(size = 26, vjust = 2.5),
                   axis.text.x = element_text(size = 24, colour = "black"),
                   axis.text.y = element_text(size = 24, colour = "black"),
                   axis.title.x = element_text(size = 26, vjust = 2.5),
                   panel.grid.major.x = element_blank(),
                   legend.text = element_text(size = 40, hjust = 0.5),
                   legend.title = element_text(size = 45, hjust = 0.5),
                   legend.key.height = unit(2,"cm"))
mac_avg


# join plots --------------------------------------------------------------

joined_plot <- plot_grid(cldn1_epi,cldn1_lp, cldn1_avg,
                         cldn2_epi, cldn2_lp, cldn2_avg,
                         ocldn_epi, ocldn_lp, ocldn_avg,
                         ecad_epi,ecad_lp, ecad_avg,
                         mac_epi, mac_lp, mac_avg,
                         axis = "btrl", align = "hv", ncol = 3 )
joined_plot

