
# Start Here --------------------------------------------------------------
##TCFA Histology Data Analysis
##Author:JDY Started:30Aug2024 
#Notes: This data set contains H scores from three independent scorers. First step is to look for scorer bias then look at differences in tissue type

# Library -----------------------------------------------------------------
#set working directory 
setwd("C:/Users/danie/OneDrive - West Texas A and M University/Research/TCFA Challenge/IHC Scoring/R Things")

library(readxl)
library(ggplot2)
library(lme4)
library(lmerTest)
library(car)
library(emmeans)
library(dplyr)
library(tidyr)
library(janitor)
library(cowplot)

# Custom color palettes ---------------------------------------------------

sample_site_palette <- c("red4", "gold","navy","forestgreen", "orange", "violet")
scorer_palette <- c("forestgreen","orange","violet")

# CLDN1 -------------------------------------------------------------------
#read in data and clean the names
cldn1_dat <- read_excel("TCFA_IHC_PathologistHScores_Cleaned_for_R.xlsx", sheet = "CLDN1")
cldn1_dat <- clean_names(cldn1_dat)

####Estimate scorer influence
#look at data
#epithelium
hist(cldn1_dat$epithelium_h_score[cldn1_dat$scorer %in% c("JC")])
hist(cldn1_dat$epithelium_h_score[cldn1_dat$scorer %in% c("JE")])
hist(cldn1_dat$epithelium_h_score[cldn1_dat$scorer %in% c("KL")])

#lp
hist(cldn1_dat$lp_h_score[cldn1_dat$scorer %in% c("JC")])
hist(cldn1_dat$lp_h_score[cldn1_dat$scorer %in% c("JE")])
hist(cldn1_dat$lp_h_score[cldn1_dat$scorer %in% c("KL")])

#look for scorer bias
epi_scoreer_model <- lm(epithelium_h_score ~ scorer, data = cldn1_dat)
summary(epi_scoreer_model)
Anova(epi_scoreer_model)
kruskal.test(cldn1_dat$epithelium_h_score~cldn1_dat$scorer)

lp_scorer_model <- lm(lp_h_score ~ scorer, data = cldn1_dat)
summary(lp_scorer_model)
Anova(lp_scorer_model)
kruskal.test(cldn1_dat$lp_h_score~cldn1_dat$scorer)

#look at data
hist(cldn1_dat$epithelium_h_score)
hist(cldn1_dat$lp_h_score)


#order organ variable so that it matches what i want to plot
cldn1_dat$organ <- factor(cldn1_dat$organ, levels = c("Rumen","Small Intestine", "Large Intestine"))

#run a covariate model
cldn1_epi_model <- lm(epithelium_h_score ~ organ + scorer, data = cldn1_dat)
summary(cldn1_epi_model)
emmeans(cldn1_epi_model, pairwise~organ)
Anova(cldn1_epi_model)

#plot
cldn1_epi_box <- ggplot(cldn1_dat, aes(x= organ, y= epithelium_h_score, fill = organ)) +
  theme_bw() + labs(y= "H-Score", x="", title = "Claudin 1") +
  #scale_x_discrete(limits=c("Rumen", "Small Intestine", "Large Intestine"))+
  geom_boxplot(aes(color = organ),alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill = scorer), shape = 21, size = 3, position = position_jitter(seed = 1, w= .2))+
  scale_color_manual(values = sample_site_palette)+
  scale_fill_manual(values = sample_site_palette, labels = c("RU. TI", "S.I. TI", "L.I. TI", "1","2","3"))+
  coord_cartesian(ylim = c(0,300), xlim= c(1,3))+
  annotate("text", x = 3.3, y = -5, label = expression(paste(italic("P"), "= 0.02")), size = 5, color = "black")+
  annotate("rect", xmin = 3, xmax = 4, ymin = -50, ymax = -1, color = "black", alpha = 0)+
  annotate("text", x = 1, y= 300, label = "b", size = 5, color= "black")+
  annotate("text", x = 2, y= 300, label = "a", size = 5, color= "black")+
  annotate("text", x = 3, y= 300, label = "ab", size = 5, color= "black")+
  theme(legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 30),
    axis.title.y = element_text(size = 26, vjust = 2.5),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.text = element_text(size = 40, hjust = 0.5),
    legend.title = element_text(size = 45, hjust = 0.5),
    legend.key.height = unit(2,"cm"))
cldn1_epi_box

#lp
cldn1_lp_model <- lm(lp_h_score ~ organ + scorer, data = cldn1_dat)
summary(cldn1_lp_model)
emmeans(cldn1_lp_model, pairwise ~ organ)
anova(cldn1_lp_model)

cldn1_lp_box <- ggplot(cldn1_dat, aes(x= organ, y= lp_h_score, fill = organ)) +
  theme_bw() + labs(y= "H-Score", x="", title = "Claudin 1") +
  #scale_x_discrete(limits=c("Rumen", "Small Intestine", "Large Intestine"))+
  geom_boxplot(aes(color = organ),alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill = scorer), shape = 21, size = 3, position = position_jitter(seed = 1, w= .2))+
  scale_color_manual(values = sample_site_palette)+
  scale_fill_manual(values = sample_site_palette)+
  coord_cartesian(ylim = c(0,300), xlim= c(1,3))+
  annotate("text", x = 3.3, y = -5, label = expression(paste(italic("P"), " < 0.01")), size = 5, color = "black")+
  annotate("rect", xmin = 3, xmax = 4, ymin = -50, ymax = -1, color = "black", alpha = 0)+
  annotate("text", x = 1, y= 250, label = "b", size = 5, color= "black")+
  annotate("text", x = 2, y= 250, label = "a", size = 5, color= "black")+
  annotate("text", x = 3, y= 250, label = "a", size = 5, color= "black")+
  theme(legend.position = "none",
        plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 26, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 40, hjust = 0.5),
        legend.title = element_text(size = 45, hjust = 0.5),
        legend.key.height = unit(2,"cm"))

cldn1_lp_box

#read in meta data
metadata <- read_excel("metadata.xlsx")
metadata$animal <- factor(metadata$animal)

cldn1_data <- merge(metadata,cldn1_dat)

#check liver
cldn1_liver_epi <- lm(epithelium_h_score~ LA_Y_N + scorer, data = cldn1_data)
summary(cldn1_liver_epi)
emmeans(cldn1_liver_epi, pairwise~ LA_Y_N)
anova(cldn1_liver_epi)

cldn1_liver_lp <- lm(lp_h_score~ LA_Y_N + scorer, data = cldn1_data)
summary(cldn1_liver_lp)
emmeans(cldn1_liver_lp, pairwise~ LA_Y_N)

#create an average and median
cldn1_dat_avg <- cldn1_dat %>% group_by(sample_id, organ) %>% 
  summarise(mean_epi= mean(epithelium_h_score, na.rm = T),
            median_epi= median(epithelium_h_score, na.rm = T),
            mean_lp= mean(lp_h_score, na.rm = T),
            median_lp= median(lp_h_score, na.rm = T))


#calculate avg and 95% CI
cldn1_avg <-cldn1_dat_avg %>% group_by(organ) %>% 
  summarise(avg_epi = mean(mean_epi, na.rm = T),
            sd_epi = sd(mean_epi, na.rm = T),
            lower_ci_epi = avg_epi - qt(0.975, df = n() - 1) * sd_epi / sqrt(n()),
            upper_ci_epi = avg_epi + qt(0.975, df = n() - 1) * sd_epi / sqrt(n()),
            avg_lp = mean(mean_lp, na.rm = T),
            sd_lp = sd(mean_lp, na.rm = T),
            lower_ci_lp = avg_lp - qt(0.975, df = n() - 1) * sd_lp / sqrt(n()),
            upper_ci_lp = avg_lp + qt(0.975, df = n() - 1) * sd_lp / sqrt(n()))



### plot 
#Box plot
cldn1_epi_h_boxplot <- ggplot(cldn1_dat, aes(x= organ, y= epithelium_h_score, fill = organ)) +
  theme_bw() + labs(y= "", x="", title = "Claudin 1 Epi H Score") +
  scale_x_discrete(limits=c("Rumen", "Small Intestine", "Large Intestine"))+
  geom_boxplot() +
  geom_point()+
  scale_fill_manual(values = sample_site_palette)+
  theme(legend.position = "none",
        plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 30, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())
cldn1_epi_h_boxplot

#raincloud
#raincloud
cldn1_epi_h_rain_cloud<- ggplot(cldn1_dat, mapping=aes(x = organ, y = epithelium_h_score, group = organ, fill = organ))+
  ggdist::stat_halfeye(adjust = 0.8, width = 0.6, .width = 0, justification = -0.2,point_colour = NA)+
  geom_boxplot(width = 0.15,outlier.shape = NA)+
  ggdist::stat_dots(side = "left", justification = 1.1, binwidth = 3, dotsize = 1)+
  coord_flip(ylim = c(0,300), xlim = c(1.3,3.1))+theme_minimal()+scale_fill_manual(values = sample_site_palette)+
  ylab("H Score")+  xlab("Sample Type") +labs(title = "CLDN1")+ 
  scale_x_discrete(limits = c("Large Intestine", "Small Intestine","Rumen"))+
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"cm"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 30, vjust = 2.5),
        axis.title.x = element_text(size = 30, vjust = 1),
        axis.text.x =  element_text(size = 24, colour = "black"),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.ticks = element_line(size = 1, colour = "black"),
        panel.grid.major.x = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 2) )
cldn1_epi_h_rain_cloud


#boxplot

cldn1_lp_h_boxplot <- ggplot(cldn1_dat, aes(x= organ, y= lp_h_score, fill = organ)) +
  theme_bw() + labs(y= "", x="", title = "Claudin 1 LP H Score") +
  scale_x_discrete(limits=c("Rumen", "Small Intestine", "Large Intestine"))+
  geom_boxplot() +
  geom_point()+
  scale_fill_manual(values = sample_site_palette)+
  theme(legend.position = "none",
        plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 30, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())
cldn1_lp_h_boxplot

#rain cloud 
cldn1_lp_h_rain_cloud<- ggplot(cldn1_dat, mapping=aes(x = organ, y = lp_h_score, group = organ, fill = organ))+
  ggdist::stat_halfeye(adjust = 0.8, width = 0.6, .width = 0, justification = -0.2,point_colour = NA)+
  geom_boxplot(width = 0.15,outlier.shape = NA)+
  ggdist::stat_dots(side = "left", justification = 1.1, binwidth = 3, dotsize = 1)+
  coord_flip(ylim = c(0,300), xlim = c(1.3,3.1))+theme_minimal()+scale_fill_manual(values = sample_site_palette)+
  ylab("H Score")+  xlab("Sample Type") +labs(title = "CLDN1")+ 
  scale_x_discrete(limits = c("Large Intestine", "Small Intestine","Rumen"))+
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"cm"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 30, vjust = 2.5),
        axis.title.x = element_text(size = 30, vjust = 1),
        axis.text.x =  element_text(size = 24, colour = "black"),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.ticks = element_line(size = 1, colour = "black"),
        panel.grid.major.x = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 2) )
cldn1_lp_h_rain_cloud


# CLDN2 -------------------------------------------------------------------
cldn2_dat <- read_excel("TCFA_IHC_PathologistHScores_Cleaned_for_R.xlsx", sheet = "CLDN2")
cldn2_dat <- clean_names(cldn2_dat)

#look at data
#epithelium
hist(cldn2_dat$epithelium_h_score[cldn2_dat$scorer %in% c("JC")])
hist(cldn2_dat$epithelium_h_score[cldn2_dat$scorer %in% c("JE")])
hist(cldn2_dat$epithelium_h_score[cldn2_dat$scorer %in% c("KL")])

#lp
hist(cldn2_dat$lp_h_score[cldn2_dat$scorer %in% c("JC")])
hist(cldn2_dat$lp_h_score[cldn2_dat$scorer %in% c("JE")])
hist(cldn2_dat$lp_h_score[cldn2_dat$scorer %in% c("KL")])

#look for scorer bias
cldn2_epi_scorer_model <- lm(epithelium_h_score ~ scorer, data = cldn2_dat)
summary(cldn2_epi_scorer_model)
Anova(cldn2_epi_scorer_model)
kruskal.test(cldn2_dat$epithelium_h_score~cldn2_dat$scorer)

cldn2_lp_scorer_model <- lm(lp_h_score ~ scorer, data = cldn2_dat)
summary(cldn2_lp_scorer_model)
Anova(cldn2_lp_scorer_model)
kruskal.test(cldn2_dat$lp_h_score~cldn2_dat$scorer)

#look at full data
hist(cldn2_dat$epithelium_h_score)
hist(cldn2_dat$lp_h_score)

#set order of organ
cldn2_dat$organ <- factor(cldn2_dat$organ, levels = c("Rumen", "Small Intestine", "Large Intestine"))

# epithelium model
cldn2_epi_model <- lm(epithelium_h_score ~ organ + scorer, data = cldn2_dat)
summary(cldn2_epi_model)
emmeans(cldn2_epi_model, pairwise~ organ)
anova(cldn2_epi_model)

cldn2_epi_box <- ggplot(cldn2_dat, aes(x= organ, y= epithelium_h_score, fill = organ)) +
  theme_bw() + labs(y= "H-Score", x="", title = "Claudin 2") +
  #scale_x_discrete(limits=c("Rumen", "Small Intestine", "Large Intestine"))+
  geom_boxplot(aes(color = organ),alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill = scorer), shape = 21, size = 3, position = position_jitter(seed = 1, w= .2))+
  scale_color_manual(values = sample_site_palette)+
  scale_fill_manual(values = sample_site_palette)+
  coord_cartesian(ylim = c(0,300), xlim= c(1,3))+
  annotate("text", x = 3.3, y = -5, label = expression(paste(italic("P"), " < 0.01")), size = 5, color = "black")+
  annotate("rect", xmin = 3, xmax = 4, ymin = -50, ymax = -1, color = "black", alpha = 0)+
  annotate("text", x = 1, y= 300, label = "b", size = 5, color= "black")+
  annotate("text", x = 2, y= 300, label = "a", size = 5, color= "black")+
  annotate("text", x = 3, y= 300, label = "ab", size = 5, color= "black")+
  theme(legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 30),
    axis.title.y = element_text(size = 26, vjust = 2.5),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.text = element_text(size = 40, hjust = 0.5),
    legend.title = element_text(size = 45, hjust = 0.5),
    legend.key.height = unit(2,"cm"))
cldn2_epi_box

#LP
cldn2_lp_model <- lm(lp_h_score ~ organ + scorer, data = cldn2_dat)
summary(cldn2_lp_model)
emmeans(cldn2_lp_model, pairwise ~ organ)
anova(cldn2_lp_model)


cldn2_lp_plot <- ggplot(cldn2_dat, aes(x= organ, y= lp_h_score, fill = organ)) +
  theme_bw() + labs(y= "H-Score", x="", title = "Claudin 2") +
  #scale_x_discrete(limits=c("Rumen", "Small Intestine", "Large Intestine"))+
  geom_boxplot(aes(color = organ),alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill = scorer), shape = 21, size = 3, position = position_jitter(seed = 1, w= .2))+
  scale_color_manual(values = sample_site_palette)+
  scale_fill_manual(values = sample_site_palette)+
  coord_cartesian(ylim = c(0,300), xlim= c(1,3))+
  annotate("text", x = 3.3, y = -5, label = expression(paste(italic("P"), " < 0.01")), size = 5, color = "black")+
  annotate("rect", xmin = 3, xmax = 4, ymin = -50, ymax = -1, color = "black", alpha = 0)+
  annotate("text", x = 1, y= 100, label = "b", size = 5, color= "black")+
  annotate("text", x = 2, y= 100, label = "a", size = 5, color= "black")+
  annotate("text", x = 3, y= 100, label = "b", size = 5, color= "black")+
  theme(legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 30),
    axis.title.y = element_text(size = 26, vjust = 2.5),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.text = element_text(size = 40, hjust = 0.5),
    legend.title = element_text(size = 45, hjust = 0.5),
    legend.key.height = unit(2,"cm"))
cldn2_lp_plot

#merge
cldn2_data <- merge(metadata, cldn2_dat)

#check livers
cldn2_liver_epi <- lm(epithelium_h_score ~ LA_Y_N + scorer, cldn2_data)
summary(cldn2_liver_epi)
emmeans(cldn2_liver_epi, pairwise~ LA_Y_N)

cldn2_liver_lp <- lm(lp_h_score ~ LA_Y_N + scorer, data = cldn2_data)
summary(cldn2_liver_lp)
emmeans(cldn2_liver_lp, pairwise~LA_Y_N)

#create an average and median
cldn2_dat_avg<- cldn2_dat %>% group_by(sample_id, organ) %>% 
  summarise(mean_epi= mean(epithelium_h_score, na.rm = T),
            median_epi= median(epithelium_h_score, na.rm = T),
            mean_lp= mean(lp_h_score, na.rm = T),
            median_lp= median(lp_h_score, na.rm = T))

#calculate avg and 95% CI
cldn2_avg <-cldn2_dat_avg %>% group_by(organ) %>% 
  summarise(avg_epi = mean(mean_epi, na.rm = T),
            sd_epi = sd(mean_epi, na.rm = T),
            lower_ci_epi = avg_epi - qt(0.975, df = n() - 1) * sd_epi / sqrt(n()),
            upper_ci_epi = avg_epi + qt(0.975, df = n() - 1) * sd_epi / sqrt(n()),
            avg_lp = mean(mean_lp, na.rm = T),
            sd_lp = sd(mean_lp, na.rm = T),
            lower_ci_lp = avg_lp - qt(0.975, df = n() - 1) * sd_lp / sqrt(n()),
            upper_ci_lp = avg_lp + qt(0.975, df = n() - 1) * sd_lp / sqrt(n()))


# OCLDN -------------------------------------------------------------------
#read data and clean names
ocldn_dat <- read_excel("TCFA_IHC_PathologistHScores_Cleaned_for_R.xlsx", sheet = "OCLDN")
ocldn_dat <- clean_names(ocldn_dat)

####Estimate scorer influence
#look at data
#epithelium
hist(ocldn_dat$epithelium_h_score[ocldn_dat$scorer %in% c("JC")])
hist(ocldn_dat$epithelium_h_score[ocldn_dat$scorer %in% c("JE")])
hist(ocldn_dat$epithelium_h_score[ocldn_dat$scorer %in% c("KL")])

#lp
hist(ocldn_dat$lp_h_score[ocldn_dat$scorer %in% c("JC")])
hist(ocldn_dat$lp_h_score[ocldn_dat$scorer %in% c("JE")])
hist(ocldn_dat$lp_h_score[ocldn_dat$scorer %in% c("KL")])

#look for scorer bias
ocldn_epi_scorer_model <- lm(epithelium_h_score ~ scorer, data = ocldn_dat)
summary(ocldn_epi_scorer_model)
Anova(ocldn_epi_scorer_model)
kruskal.test(ocldn_dat$epithelium_h_score~ocldn_dat$scorer)

ocldn_lp_scorer_model <- lm(lp_h_score ~ scorer, data = ocldn_dat)
summary(ocldn_lp_scorer_model)
Anova(ocldn_lp_scorer_model)
kruskal.test(ocldn_dat$lp_h_score~ocldn_dat$scorer)

#visualize data
hist(ocldn_dat$epithelium_h_score)
hist(ocldn_dat$lp_h_score) 

#order variables
ocldn_dat$organ <- factor(ocldn_dat$organ, levels = c("Rumen", "Small Intestine", "Large Intestine"))

#epi
ocldn_epi_model <- lm(epithelium_h_score ~ organ + scorer, data = ocldn_dat)
summary(ocldn_epi_model)
lsmeans(ocldn_epi_model, pairwise ~ organ)
anova(ocldn_epi_model)

plot(residuals(ocldn_epi_model))
plot(ocldn_epi_model)

ocldn_epi_plot <- ggplot(ocldn_dat, aes(x= organ, y= epithelium_h_score, fill = organ)) +
  theme_bw() + labs(y= "H-Score", x="", title = "Occludin") +
  #scale_x_discrete(limits=c("Rumen", "Small Intestine", "Large Intestine"))+
  geom_boxplot(aes(color = organ),alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill = scorer), shape = 21, size = 3, position = position_jitter(seed = 1, w= .2))+
  scale_color_manual(values = sample_site_palette)+
  scale_fill_manual(values = sample_site_palette)+
  coord_cartesian(ylim = c(0,300), xlim= c(1,3))+
  annotate("text", x = 3.3, y = -5, label = expression(paste(italic("P"), " = 0.18")), size = 5, color = "black")+
  annotate("rect", xmin = 3, xmax = 4, ymin = -50, ymax = -1, color = "black", alpha = 0)+
  #annotate("text", x = 1, y= 300, label = "b", size = 5, color= "black")+
  #annotate("text", x = 2, y= 300, label = "a", size = 5, color= "black")+
  #annotate("text", x = 3, y= 300, label = "ab", size = 5, color= "black")+
  theme(legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 30),
    axis.title.y = element_text(size = 26, vjust = 2.5),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.text = element_text(size = 40, hjust = 0.5),
    legend.title = element_text(size = 45, hjust = 0.5),
    legend.key.height = unit(2,"cm"))
ocldn_epi_plot

#lp
ocldn_lp_model <- lm(lp_h_score ~ organ + scorer, data = ocldn_dat)
summary(ocldn_lp_model)
emmeans(ocldn_lp_model, pairwise~organ)
anova(ocldn_lp_model)

plot(residuals(ocldn_lp_model))
plot(ocldn_lp_model)

ocldn_lp_plot <- ggplot(ocldn_dat, aes(x= organ, y= lp_h_score, fill = organ)) +
  theme_bw() + labs(y= "H-Score", x="", title = "Occludin") +
  #scale_x_discrete(limits=c("Rumen", "Small Intestine", "Large Intestine"))+
  geom_boxplot(aes(color = organ),alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill = scorer), shape = 21, size = 3, position = position_jitter(seed = 1, w= .2))+
  scale_color_manual(values = sample_site_palette)+
  scale_fill_manual(values = sample_site_palette)+
  coord_cartesian(ylim = c(0,300), xlim= c(1,3))+
  annotate("text", x = 3.3, y = -5, label = expression(paste(italic("P"), " = 0.01")), size = 5, color = "black")+
  annotate("rect", xmin = 3, xmax = 4, ymin = -50, ymax = -1, color = "black", alpha = 0)+
  annotate("text", x = 1, y= 100, label = "b", size = 5, color= "black")+
  annotate("text", x = 2, y= 100, label = "a", size = 5, color= "black")+
  annotate("text", x = 3, y= 100, label = "a", size = 5, color= "black")+
  theme(legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 30),
    axis.title.y = element_text(size = 26, vjust = 2.5),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.text = element_text(size = 40, hjust = 0.5),
    legend.title = element_text(size = 45, hjust = 0.5),
    legend.key.height = unit(2,"cm"))
ocldn_lp_plot

#merge data
ocldn_data <- merge(metadata, ocldn_dat)

#check livers
ocldn_liver_epi <- lm(epithelium_h_score~ LA_Y_N + scorer, data = ocldn_data)
summary(ocldn_liver_epi)
emmeans(ocldn_liver_epi, pairwise~ LA_Y_N)


ocldn_liver_lp <- lm(lp_h_score ~ LA_Y_N + scorer, data = ocldn_data)
summary(ocldn_liver_lp)
emmeans(ocldn_liver_lp, pairwise~LA_Y_N)

#create an average and median
ocldn_dat_avg<- ocldn_dat %>% group_by(sample_id, organ) %>% 
  summarise(mean_epi= mean(epithelium_h_score, na.rm = T),
            median_epi= median(epithelium_h_score, na.rm = T),
            mean_lp= mean(lp_h_score, na.rm = T),
            median_lp= median(lp_h_score, na.rm = T))

#calculate avg and 95% CI
ocldn_avg <-ocldn_dat_avg %>% group_by(organ) %>% 
  summarise(avg_epi = mean(mean_epi, na.rm = T),
            sd_epi = sd(mean_epi, na.rm = T),
            lower_ci_epi = avg_epi - qt(0.975, df = n() - 1) * sd_epi / sqrt(n()),
            upper_ci_epi = avg_epi + qt(0.975, df = n() - 1) * sd_epi / sqrt(n()),
            avg_lp = mean(mean_lp, na.rm = T),
            sd_lp = sd(mean_lp, na.rm = T),
            lower_ci_lp = avg_lp - qt(0.975, df = n() - 1) * sd_lp / sqrt(n()),
            upper_ci_lp = avg_lp + qt(0.975, df = n() - 1) * sd_lp / sqrt(n()))


# MAC387 ------------------------------------------------------------------
#read data and clean names
mac_dat <- read_excel("TCFA_IHC_PathologistHScores_Cleaned_for_R.xlsx", sheet = "MAC387")
mac_dat <- clean_names(mac_dat)

####Estimate scorer influence
#look at data
#epithelium
hist(mac_dat$epithelium_h_score[mac_dat$scorer %in% c("JC")])
hist(mac_dat$epithelium_h_score[mac_dat$scorer %in% c("JE")])
hist(mac_dat$epithelium_h_score[mac_dat$scorer %in% c("KL")])

#lp
hist(mac_dat$lp_h_score[mac_dat$scorer %in% c("JC")])
hist(mac_dat$lp_h_score[mac_dat$scorer %in% c("JE")])
hist(mac_dat$lp_h_score[mac_dat$scorer %in% c("KL")])

#look for scorer bias
mac_epi_scorer_model <- lm(epithelium_h_score ~ scorer, data = mac_dat)
summary(mac_epi_scorer_model)
Anova(mac_epi_scorer_model)
kruskal.test(mac_dat$epithelium_h_score~mac_dat$scorer)

mac_lp_scorer_model <- lm(lp_h_score ~ scorer, data = mac_dat)
summary(mac_lp_scorer_model)
Anova(mac_lp_scorer_model)
kruskal.test(mac_dat$lp_h_score~mac_dat$scorer)

#visualize the data
hist(mac_dat$epithelium_h_score) #lots of zeros
hist(mac_dat$lp_h_score)# lots of zeros

#order variable
mac_dat$organ <- factor(mac_dat$organ, levels = c("Rumen","Small Intestine", "Large Intestine"))

#epi model
mac_epi_model <- lm(epithelium_h_score ~ organ + scorer, data = mac_dat)
summary(mac_epi_model)
emmeans(mac_epi_model, pairwise ~ organ)
anova(mac_epi_model)

plot(residuals(mac_epi_model))
plot(mac_epi_model)

mac_epi_plot <- ggplot(mac_dat, aes(x= organ, y= epithelium_h_score, fill = organ)) +
  theme_bw() + labs(y= "H-Score", x="", title = "MAC-387") +
  #scale_x_discrete(limits=c("Rumen", "Small Intestine", "Large Intestine"))+
  geom_boxplot(aes(color = organ),alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill = scorer), shape = 21, size = 3, position = position_jitter(seed = 1, w= .2))+
  scale_color_manual(values = sample_site_palette)+
  scale_fill_manual(values = sample_site_palette)+
  coord_cartesian(ylim = c(0,300), xlim= c(1,3))+
  annotate("text", x = 3.3, y = -5, label = expression(paste(italic("P"), " < 0.01")), size = 5, color = "black")+
  annotate("rect", xmin = 3, xmax = 4, ymin = -50, ymax = -1, color = "black", alpha = 0)+
  annotate("text", x = 1, y= 300, label = "a", size = 5, color= "black")+
  annotate("text", x = 2, y= 300, label = "b", size = 5, color= "black")+
  annotate("text", x = 3, y= 300, label = "b", size = 5, color= "black")+
  theme(legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 30),
    axis.title.y = element_text(size = 26, vjust = 2.5),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.text = element_text(size = 40, hjust = 0.5),
    legend.title = element_text(size = 45, hjust = 0.5),
    legend.key.height = unit(2,"cm"))
mac_epi_plot

#lp
mac_lp_model <- lm(lp_h_score ~ organ + scorer, data = mac_dat)
summary(mac_lp_model)
emmeans(mac_lp_model, pairwise ~ organ)
anova(mac_lp_model)

plot(residuals(mac_lp_model))
plot(mac_lp_model)

mac_lp_plot <- ggplot(mac_dat, aes(x= organ, y= lp_h_score, fill = organ)) +
  theme_bw() + labs(y= "H-Score", x="", title = "MAC-387") +
  #scale_x_discrete(limits=c("Rumen", "Small Intestine", "Large Intestine"))+
  geom_boxplot(aes(color = organ), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill = scorer), shape = 21, size = 3, position = position_jitter(seed = 1, w= .2))+
  scale_color_manual(values = sample_site_palette)+
  scale_fill_manual(values = sample_site_palette)+
  coord_cartesian(ylim = c(0,300), xlim= c(1,3))+
  annotate("text", x = 3.3, y = -5, label = expression(paste(italic("P"), " < 0.01")), size = 5, color = "black")+
  annotate("rect", xmin = 3, xmax = 4, ymin = -50, ymax = -1, color = "black", alpha = 0)+
  annotate("text", x = 1, y= 150, label = "b", size = 5, color= "black")+
  annotate("text", x = 2, y= 150, label = "a", size = 5, color= "black")+
  annotate("text", x = 3, y= 150, label = "ab", size = 5, color= "black")+
  theme(legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 30),
    axis.title.y = element_text(size = 26, vjust = 2.5),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.text = element_text(size = 40, hjust = 0.5),
    legend.title = element_text(size = 45, hjust = 0.5),
    legend.key.height = unit(2,"cm"))
mac_lp_plot

#merge metadata
mac_data <- merge(metadata, mac_dat)

#check livers
mac_liver_epi <- lm(epithelium_h_score~ LA_Y_N + scorer, data = mac_data)
summary(mac_liver_epi)
emmeans(mac_liver_epi, pairwise~LA_Y_N)

mac_liver_lp <- lm(lp_h_score~ LA_Y_N + scorer, data = mac_data)
summary(mac_liver_lp)
emmeans(mac_liver_lp, pairwise~ LA_Y_N)

#create an average and median
mac_dat_avg<- mac_dat %>% group_by(sample_id, organ) %>% 
  summarise(mean_epi= mean(epithelium_h_score, na.rm = T),
            median_epi= median(epithelium_h_score, na.rm = T),
            mean_lp= mean(lp_h_score, na.rm = T),
            median_lp= median(lp_h_score, na.rm = T))

#calculate avg and 95% CI
mac_avg <-mac_dat_avg %>% group_by(organ) %>% 
  summarise(avg_epi = mean(mean_epi, na.rm = T),
            sd_epi = sd(mean_epi, na.rm = T),
            lower_ci_epi = avg_epi - qt(0.975, df = n() - 1) * sd_epi / sqrt(n()),
            upper_ci_epi = avg_epi + qt(0.975, df = n() - 1) * sd_epi / sqrt(n()),
            avg_lp = mean(mean_lp, na.rm = T),
            sd_lp = sd(mean_lp, na.rm = T),
            lower_ci_lp = avg_lp - qt(0.975, df = n() - 1) * sd_lp / sqrt(n()),
            upper_ci_lp = avg_lp + qt(0.975, df = n() - 1) * sd_lp / sqrt(n()))

# E-cad -------------------------------------------------------------------
#read data and clean names
ecad_dat <- read_excel("TCFA_IHC_PathologistHScores_Cleaned_for_R.xlsx", sheet = "E-Cadherin")
ecad_dat <- clean_names(ecad_dat)

####Estimate scorer influence
#look at data
#epithelium
hist(ecad_dat$epithelium_h_score[ecad_dat$scorer %in% c("JC")])
hist(ecad_dat$epithelium_h_score[ecad_dat$scorer %in% c("JE")])
hist(ecad_dat$epithelium_h_score[ecad_dat$scorer %in% c("KL")])

#lp
hist(ecad_dat$lp_h_score[ecad_dat$scorer %in% c("JC")])
hist(ecad_dat$lp_h_score[ecad_dat$scorer %in% c("JE")])
hist(ecad_dat$lp_h_score[ecad_dat$scorer %in% c("KL")])

#look for scorer bias
ecad_epi_scorer_model <- lm(epithelium_h_score ~ scorer, data = ecad_dat)
summary(ecad_epi_scorer_model)
Anova(ecad_epi_scorer_model)
kruskal.test(ecad_dat$epithelium_h_score~ecad_dat$scorer)

ecad_lp_scorer_model <- lm(lp_h_score ~ scorer, data = ecad_dat)
summary(ecad_lp_scorer_model)
Anova(ecad_lp_scorer_model)
kruskal.test(ecad_dat$lp_h_score~ecad_dat$scorer)

#visualize the data
hist(ecad_dat$epithelium_h_score)
hist(ecad_dat$lp_h_score)

#order variables
ecad_dat$organ <- factor(ecad_dat$organ, levels = c("Rumen", "Small Intestine", "Large Intestine"))

#epi model
ecad_epi_model <- lm(epithelium_h_score ~ organ + scorer, data = ecad_dat)
summary(ecad_epi_model)
emmeans(ecad_epi_model, pairwise~organ)
anova(ecad_epi_model)

plot(residuals(ecad_epi_model))
plot(ecad_epi_model)

ecad_epi_plot <- ggplot(ecad_dat, aes(x= organ, y= epithelium_h_score, fill = organ)) +
  theme_bw() + labs(y= "H-Score", x="", title = "E-Cadherin") +
  #scale_x_discrete(limits=c("Rumen", "Small Intestine", "Large Intestine"))+
  geom_boxplot(aes(color = organ),alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill = scorer), shape = 21, size = 3, position = position_jitter(seed = 1, w= .2))+
  scale_color_manual(values = sample_site_palette)+
  scale_fill_manual(values = sample_site_palette)+
  coord_cartesian(ylim = c(0,300), xlim= c(1,3))+
  annotate("text", x = 3.3, y = -5, label = expression(paste(italic("P"), " < 0.01")), size = 5, color = "black")+
  annotate("rect", xmin = 3, xmax = 4, ymin = -50, ymax = -1, color = "black", alpha = 0)+
  annotate("text", x = 1, y= 310, label = "b", size = 5, color= "black")+
  annotate("text", x = 2, y= 310, label = "a", size = 5, color= "black")+
  annotate("text", x = 3, y= 310, label = "ab", size = 5, color= "black")+
  theme(legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 30),
    axis.title.y = element_text(size = 26, vjust = 2.5),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.text = element_text(size = 40, hjust = 0.5),
    legend.title = element_text(size = 45, hjust = 0.5),
    legend.key.height = unit(2,"cm"))
ecad_epi_plot

#lp
ecad_lp_model <- lm(lp_h_score ~ organ + scorer, data = ecad_dat)
summary(ecad_lp_model)
emmeans(ecad_lp_model, pairwise ~ organ)
anova(ecad_lp_model)

plot(residuals(ecad_lp_model))
plot(ecad_lp_model)

ecad_lp_plot <- ggplot(ecad_dat, aes(x= organ, y= lp_h_score, fill = organ)) +
  theme_bw() + labs(y= "H-Score", x="", title = "E-Cadherin") +
  #scale_x_discrete(limits=c("Rumen", "Small Intestine", "Large Intestine"))+
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill = scorer), shape = 21, size = 3, position = position_jitter(seed = 1, w= .2))+
  scale_color_manual(values = sample_site_palette)+
  scale_fill_manual(values = sample_site_palette)+
  coord_cartesian(ylim = c(0,300), xlim= c(1,3))+
  annotate("text", x = 3.3, y = -5, label = expression(paste(italic("P"), " < 0.01")), size = 5, color = "black")+
  annotate("rect", xmin = 3, xmax = 4, ymin = -50, ymax = -1, color = "black", alpha = 0)+
  annotate("text", x = 1, y= 150, label = "b", size = 5, color= "black")+
  annotate("text", x = 2, y= 150, label = "a", size = 5, color= "black")+
  annotate("text", x = 3, y= 150, label = "a", size = 5, color= "black")+
  theme(legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 30),
    axis.title.y = element_text(size = 26, vjust = 2.5),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.text = element_text(size = 40, hjust = 0.5),
    legend.title = element_text(size = 45, hjust = 0.5),
    legend.key.height = unit(2,"cm"))
ecad_lp_plot

#merge meta data
ecad_data <- merge(metadata, ecad_dat)

#check livers
ecad_liver_epi <- lm(epithelium_h_score ~ LA_Y_N + scorer, data = ecad_data )
summary(ecad_liver_epi)
emmeans(ecad_liver_epi, pairwise~LA_Y_N)

ecad_liver_lp <- lm(lp_h_score ~ LA + scorer, data = ecad_data)
summary(ecad_liver_lp)
emmeans(ecad_liver_lp, pairwise~ LA)
#only flukes dont worry about it

#create an average and median
ecad_dat_avg<- ecad_dat %>% group_by(sample_id, organ) %>% 
  summarise(mean_epi= mean(epithelium_h_score, na.rm = T),
            median_epi= median(epithelium_h_score, na.rm = T),
            mean_lp= mean(lp_h_score, na.rm = T),
            median_lp= median(lp_h_score, na.rm = T))

#calculate avg and 95% CI
ecad_avg <-ecad_dat_avg %>% group_by(organ) %>% 
  summarise(avg_epi = mean(mean_epi, na.rm = T),
            sd_epi = sd(mean_epi, na.rm = T),
            lower_ci_epi = avg_epi - qt(0.975, df = n() - 1) * sd_epi / sqrt(n()),
            upper_ci_epi = avg_epi + qt(0.975, df = n() - 1) * sd_epi / sqrt(n()),
            avg_lp = mean(mean_lp, na.rm = T),
            sd_lp = sd(mean_lp, na.rm = T),
            lower_ci_lp = avg_lp - qt(0.975, df = n() - 1) * sd_lp / sqrt(n()),
            upper_ci_lp = avg_lp + qt(0.975, df = n() - 1) * sd_lp / sqrt(n()))


# make legends ------------------------------------------------------------

sample_legend <- ggplot(ecad_dat, aes(x= organ, y= lp_h_score, fill = organ)) +
  theme_bw() + labs(y= "H-Score", x="", title = "E-Cadherin") +
  #scale_x_discrete(limits=c("Rumen", "Small Intestine", "Large Intestine"))+
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  #geom_point(aes(fill = scorer), shape = 21, size = 3, position = position_jitter(seed = 1, w= .2))+
  scale_color_manual(values = sample_site_palette, labels = c("RU TI", "SI TI","LI TI"))+
  guides(fill = guide_legend(title = "Sample Type", byrow = TRUE), color = guide_legend(title = "Sample Type", byrow = TRUE))+
  scale_fill_manual(values = sample_site_palette,labels = c("RU TI", "SI TI","LI TI"))+
  coord_cartesian(ylim = c(0,300), xlim= c(1,3))+
  annotate("text", x = 3.3, y = -5, label = expression(paste(italic("P"), " < 0.01")), size = 5, color = "black")+
  annotate("rect", xmin = 3, xmax = 4, ymin = -50, ymax = -1, color = "black", alpha = 0)+
  annotate("text", x = 1, y= 150, label = "b", size = 5, color= "black")+
  annotate("text", x = 2, y= 150, label = "a", size = 5, color= "black")+
  annotate("text", x = 3, y= 150, label = "a", size = 5, color= "black")+
  theme(#legend.position = "none",
        plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 26, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 40, hjust = 0.5),
        legend.title = element_text(size = 45, hjust = 0.5),
        legend.key.height = unit(2,"cm"))
sample_legend

sample_legend <- get_legend(sample_legend)


scorer_legend_plot <- ggplot(ecad_dat, aes(x= organ, y= lp_h_score, fill = organ)) +
  theme_bw() + labs(y= "H-Score", x="", title = "E-Cadherin") +
  #scale_x_discrete(limits=c("Rumen", "Small Intestine", "Large Intestine"))+
  #geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill = scorer), shape = 21, size = 10, position = position_jitter(seed = 1, w= .2))+
  scale_color_manual(values = scorer_palette, labels = c("1","2","3"))+
  guides(fill = guide_legend(title = "Pathologist", byrow = TRUE), color = guide_legend(title = "Pathologist", byrow = TRUE))+
  scale_fill_manual(values = scorer_palette, labels = c("1","2","3"))+
  coord_cartesian(ylim = c(0,300), xlim= c(1,3))+
  annotate("text", x = 3.3, y = -5, label = expression(paste(italic("P"), " < 0.01")), size = 5, color = "black")+
  annotate("rect", xmin = 3, xmax = 4, ymin = -50, ymax = -1, color = "black", alpha = 0)+
  annotate("text", x = 1, y= 150, label = "b", size = 5, color= "black")+
  annotate("text", x = 2, y= 150, label = "a", size = 5, color= "black")+
  annotate("text", x = 3, y= 150, label = "a", size = 5, color= "black")+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 30),
    axis.title.y = element_text(size = 26, vjust = 2.5),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.text = element_text(size = 40, hjust = 0.5),
    legend.title = element_text(size = 45, hjust = 0.5),
    legend.key.height = unit(2,"cm"))
scorer_legend_plot

scorer_legend <- get_legend(scorer_legend_plot)

# Joining plots -----------------------------------------------------------
#join legend
legend <- plot_grid(sample_legend, scorer_legend, align = "hv", axis = "btrl", ncol = 1)
legend

#join epi plots
plot_grid(cldn1_epi_box, cldn2_epi_box, ocldn_epi_plot, ecad_epi_plot, mac_epi_plot, legend, align = "hv", axis = "btrl", ncol = 3)

#join lp plots
plot_grid(cldn1_lp_box, cldn2_lp_plot, ocldn_lp_plot, ecad_lp_plot, mac_lp_plot, legend, align = "hv", axis = "btrl", ncol = 3 )
