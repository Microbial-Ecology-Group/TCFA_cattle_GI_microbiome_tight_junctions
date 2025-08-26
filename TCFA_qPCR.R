# Start Here --------------------------------------------------------------
#qPCR Data, TCFA Project
#Author: JDY Started: 29Aug24 Finished: 17OCt24
#Description: this data set contains qPCR results from 5 different tight junction proteins.
#Note: the get legend function only works if you plot the plot one time with the legend, extract the legend and then re-run the plot with out the legend and then merge the plots later on


# Library -----------------------------------------------------------------
library(readxl)
library(ggplot2)
library(lme4)
library(car)
library(emmeans)
library(dplyr)
library(tidyr)
library(janitor)
library(cowplot)

#set working directory
setwd("C:/Users/danie/OneDrive - West Texas A and M University/Research/TCFA Challenge/qPCR Data")

#build a palette
sample_site_palette <- c( "red4", "gold", "navy")

#set levels to data for consistency
dat$sample_site <- factor(dat$sample_site, levels = c("RUTI", "SITI", "LITI"))
levels(dat$sample_site)
# Import Data -------------------------------------------------------------

dat <- read_excel("Extracted_qPCR_Data.xlsx")
#View(dat)

#look at data
hist(dat$CLDN1_mean)
hist(dat$CLDN2_mean)
hist(dat$OCLDN_mean)
hist(dat$TJP1_mean)
hist(dat$CDH1_mean)

########CLDN1
#test for normality
shapiro.test(dat$CLDN1_mean)  #not normal

#log transform
dat$cldn1_log <- log(dat$CLDN1_mean)
#re-test normality
shapiro.test(dat$cldn1_log) #still not normal

#kruskal wallis 
kruskal.test(dat$CLDN1_mean, dat$sample_site) #p = 0.017 Significant

pairwise.wilcox.test(dat$CLDN1_mean, dat$sample_site, p.adjust.method = "BH") #small intestine is different from others

#interquartile range
summary(dat$CLDN1_mean)
IQR(dat$CLDN1_mean, na.rm = T)

#calculate threshold
CLD1min <- 13.82-(1.5*1452.191)
CLDN1max <- 1466.01+(1.5*1452.191)
#find outlier
dat$CLDN1_mean[which(dat$CLDN1_mean < CLD1min | dat$CLDN1_mean > CLDN1max)]

#Box plot
cldn1_box <- ggplot(dat, aes(x= sample_site, y= CLDN1_mean, fill = sample_site)) +
  theme_bw() + labs(y= "Copies/20 ng RNA\n (log10)", x="", title = "Claudin 1") +
  scale_x_discrete(limits=c("RUTI", "SITI", "LITI"))+
  geom_boxplot(aes(color= sample_site), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(color = sample_site), size = 2, position = position_jitter(seed = 3, w= .2))+
  scale_y_log10(breaks= c(1,10,100,1000,10000,100000),
                labels=c("1","10","100","1,000","10,000","100,000"))+
  scale_fill_manual(values = sample_site_palette)+
  scale_color_manual(values = sample_site_palette)+
  annotate("text", x = 3.7, y = 0.1, label = expression(paste(italic("P"), " = 0.01")), size = 5, color = "black")+
  annotate("rect", xmin = 3.4, xmax = 4, ymin = 0, ymax = 0.16, color = "black", alpha = 0)+
  annotate("text", x = 1, y= 0.18, label = "n = 21", size = 5, color= "black")+
  annotate("text", x = 2, y= 0.18, label = "n = 20", size = 5, color= "black")+
  annotate("text", x = 3, y= 0.18, label = "n = 19", size = 5, color= "black")+
  annotate("text", x = 1, y= 9650, label = "b", size = 5, color= "black")+
  annotate("text", x = 2, y= 9650, label = "a", size = 5, color= "black")+
  annotate("text", x = 3, y= 9650, label = "ab", size = 5, color= "black")+
  theme(legend.position = "none",
        plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 26, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())
cldn1_box

#raincloud
cldn1_rain_cloud<- ggplot(dat, mapping=aes(x = sample_site, y = CLDN1_mean, group = sample_site, fill = sample_site))+
  ggdist::stat_halfeye(adjust = 0.8, width = 0.6, .width = 0, justification = -0.2,point_colour = NA)+
  geom_boxplot(width = 0.15,outlier.shape = NA)+
  ggdist::stat_dots(side = "left", justification = 1.1, binwidth = 3, dotsize = 0.017)+
  scale_y_log10()+
  coord_flip( xlim = c(1.3,3.1))+theme_minimal()+scale_fill_manual(values = sample_site_palette)+
  ylab("Copies/ 20 ng of RNA")+  xlab("Sample Type") +labs(title = "CLDN1")+ 
  scale_x_discrete(limits = c("LITI", "SITI","RUTI"))+
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
cldn1_rain_cloud

#read excel for subscript
cldn1_cld <- read_excel("CLD.xlsx")
cldn1_cld <- as.data.frame(cldn1_cld)


#do it without the log 
cldn1.2_rain_cloud<- ggplot(dat, mapping=aes(x = sample_site, y = CLDN1_mean, group = sample_site, fill = sample_site))+
  ggdist::stat_halfeye(adjust = 0.8, width = 0.6, .width = 0, justification = -0.2,point_colour = NA)+
  geom_boxplot(width = 0.15,outlier.shape = NA)+
  gghalves::geom_half_point(side = "l", range_scale = .4, alpha = .2, aes(color = sample_site) ) +
  scale_y_log10()+
  coord_flip(xlim = c(1,3.1))+theme_minimal()+scale_fill_manual(values = sample_site_palette)+
  scale_color_manual(values = sample_site_palette)+
  ylab("Copies/ 20 ng of RNA")+  xlab("Sample Type") +labs(title = "CLDN1")+ 
  scale_x_discrete(limits = c("LITI", "SITI","RUTI"))+
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
cldn1.2_rain_cloud



#test LA Y or NO
kruskal.test(dat$CLDN1_mean, dat$LA_Y_N)#no difference p = 0.2212

#test animal ID
kruskal.test(dat$CLDN1_mean, dat$animal)#no difference p = 0.1581

#test plate
kruskal.test(dat$CLDN1_mean, dat$CLDN1_plate_number)#no difference p = 0.2102

########CLDN2
#test for normality
shapiro.test(dat$CLDN2_mean)# not normal

#log transform
dat$cldn2_log <- log(dat$CLDN2_mean)

#re-test normality
shapiro.test(dat$cldn2_log)#p = 0.09 normal but will procede with nonparametric test for time being.

#kruskal wallis
kruskal.test(dat$CLDN2_mean, dat$sample_site)#p = 0.000308 Significant 

pairwise.wilcox.test(dat$CLDN2_mean, dat$sample_site, p.adjust.method = "BH") #rumen different than all other samples 

cldn2_box <- ggplot(dat, aes(x= sample_site, y= CLDN2_mean, fill = sample_site)) +
  theme_bw() + labs(y= "Copies/20 ng RNA\n (log10)", x="", title = "Claudin 2") +
  scale_x_discrete(limits=c("RUTI", "SITI", "LITI"))+
  geom_boxplot(aes(color= sample_site), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(color = sample_site), size = 2, position = position_jitter(seed = 2, w= .2))+
  scale_y_log10(breaks= c(1,10,100,1000,10000,100000),
                labels=c("1","10","100","1,000","10,000","100,000"))+
  scale_fill_manual(values = sample_site_palette)+
  scale_color_manual(values = sample_site_palette)+
  annotate("text", x = 3.7, y = 0.1, label = expression(paste(italic("P"), " < 0.01")), size = 5, color = "black")+
  annotate("rect", xmin = 3.4, xmax = 4, ymin = 0, ymax = 0.16, color = "black", alpha = 0)+
  annotate("text", x = 1, y= 0.18, label = "n = 11", size = 5, color= "black")+
  annotate("text", x = 2, y= 0.18, label = "n = 20", size = 5, color= "black")+
  annotate("text", x = 3, y= 0.18, label = "n = 17", size = 5, color= "black")+
  annotate("text", x = 1, y= 9650, label = "b", size = 5, color= "black")+
  annotate("text", x = 2, y= 9650, label = "a", size = 5, color= "black")+
  annotate("text", x = 3, y= 9650, label = "a", size = 5, color= "black")+
  theme(legend.position = "none",
        plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 26, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())
cldn2_box

#raincloud
cldn2_rain_cloud <- ggplot(dat, mapping=aes(x = sample_site, y = CLDN2_mean, group = sample_site, fill = sample_site))+
  ggdist::stat_halfeye(adjust = 0.8, width = 0.6, .width = 0, justification = -0.2,point_colour = NA)+
  geom_boxplot(width = 0.15,outlier.shape = NA)+
  gghalves::geom_half_point(side = "l", range_scale = .4, alpha = .2, aes(color = sample_site) ) +
  scale_y_log10()+
  coord_flip(xlim = c(1,3.1))+theme_minimal()+scale_fill_manual(values = sample_site_palette)+
  scale_color_manual(values = sample_site_palette)+
  ylab("Copies/ 20 ng of RNA")+  xlab("Sample Type") +labs(title = "CLDN2")+ 
  scale_x_discrete(limits = c("LITI", "SITI","RUTI"))+
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
cldn2_rain_cloud

#test LA y or n
kruskal.test(dat$CLDN2_mean, dat$LA_Y_N)# no difference p = 0.2328

#test animal
kruskal.test(dat$CLDN2_mean, dat$animal)# no difference p = 0.3482

#test plate
kruskal.test(dat$CLDN2_mean, dat$CLDN2_plate_number)# no difference p = 0.7332

########OCLDN
#test for normality
shapiro.test(dat$OCLDN_mean)#not normal

#log transform
dat$ocldn_log <- log(dat$OCLDN_mean)

#re-test normality
shapiro.test(dat$ocldn_log)#not nomral

#kruskal wallis
kruskal.test(dat$OCLDN_mean, dat$sample_site)# not significant p = 0.07633

ocldn_box <- ggplot(dat, aes(x= sample_site, y= OCLDN_mean, fill = sample_site)) +
  theme_bw() + labs(y= "Copies/20 ng RNA\n (log10)", x="", title = "Occludin") +
  scale_x_discrete(limits=c("RUTI", "SITI", "LITI"))+
  geom_boxplot(aes(color= sample_site), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(color = sample_site), size = 2, position = position_jitter(seed = 1, w= .2))+
  scale_y_log10(breaks= c(1,10,100,1000,10000,100000),
                labels=c("1","10","100","1,000","10,000","100,000"))+
  scale_fill_manual(values = sample_site_palette)+
  scale_color_manual(values = sample_site_palette)+
  annotate("text", x = 3.7, y = 0.1, label = expression(paste(italic("P"), " = 0.08")), size = 5, color = "black")+
  annotate("rect", xmin = 3.4, xmax = 4, ymin = 0, ymax = 0.16, color = "black", alpha = 0)+
  annotate("text", x = 1, y= 0.18, label = "n = 21", size = 5, color= "black")+
  annotate("text", x = 2, y= 0.18, label = "n = 21", size = 5, color= "black")+
  annotate("text", x = 3, y= 0.18, label = "n = 19", size = 5, color= "black")+
  theme(legend.position = "none",
        plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 26, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())
ocldn_box

#raincloud
ocldn_rain_cloud <- ggplot(dat, mapping=aes(x = sample_site, y = OCLDN_mean, group = sample_site, fill = sample_site))+
  ggdist::stat_halfeye(adjust = 0.8, width = 0.6, .width = 0, justification = -0.2,point_colour = NA)+
  geom_boxplot(width = 0.15,outlier.shape = NA)+
  gghalves::geom_half_point(side = "l", range_scale = .4, alpha = .2, aes(color = sample_site) ) +
  scale_y_log10()+
  coord_flip(xlim = c(1,3.1))+theme_minimal()+scale_fill_manual(values = sample_site_palette)+
  scale_color_manual(values = sample_site_palette)+
  ylab("Copies/ 20 ng of RNA")+  xlab("Sample Type") +labs(title = "OCLDN")+ 
  scale_x_discrete(limits = c("LITI", "SITI","RUTI"))+
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
ocldn_rain_cloud

#test LA y or n
kruskal.test(dat$OCLDN_mean, dat$LA_Y_N)# no difference p = 0.5598

#test animal
kruskal.test(dat$OCLDN_mean, dat$animal)# no difference p = 0.2892

#test plate
kruskal.test(dat$OCLDN_mean, dat$OCLDN_plate_number)# no difference p = 0.07363

########TJP1
#test for normality
shapiro.test(dat$TJP1_mean)#not normal

#log transform
dat$TJP1_log <- log(dat$TJP1_mean)

#re-test
shapiro.test(dat$TJP1_log)#0.04468 could round. will procede with non-parametric tests.

#kruskal wallace
kruskal.test(dat$TJP1_mean,dat$sample_site)#p = 0.02787

#pairwise
pairwise.wilcox.test(dat$TJP1_mean, dat$sample_site, p.adjust.method = "BH") # SI different from LI p = 0.01

tjp1_box <- ggplot(dat, aes(x= sample_site, y= TJP1_mean, fill = sample_site)) +
  theme_bw() + labs(y= "Copies/20 ng RNA\n (log10)", x="", title = "Tight Junction Protein 1") +
  scale_x_discrete(limits=c("RUTI", "SITI", "LITI"))+
  geom_boxplot(aes(color= sample_site), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(color = sample_site), size = 2, position = position_jitter(seed = 1, w= .2))+
  scale_y_log10(breaks= c(1,10,100,1000,10000,100000),
                labels=c("1","10","100","1,000","10,000","100,000"))+
  scale_fill_manual(values = sample_site_palette)+
  scale_color_manual(values = sample_site_palette)+
  annotate("text", x = 3.7, y = 0.1, label = expression(paste(italic("P"), " = 0.03")), size = 5, color = "black")+
  annotate("rect", xmin = 3.4, xmax = 4, ymin = 0, ymax = 0.16, color = "black", alpha = 0)+
  annotate("text", x = 1, y= 0.18, label = "n = 21", size = 5, color= "black")+
  annotate("text", x = 2, y= 0.18, label = "n = 20", size = 5, color= "black")+
  annotate("text", x = 3, y= 0.18, label = "n = 20", size = 5, color= "black")+
  annotate("text", x = 1, y= 9650, label = "ab", size = 5, color= "black")+
  annotate("text", x = 2, y= 9650, label = "b", size = 5, color= "black")+
  annotate("text", x = 3, y= 9650, label = "a", size = 5, color= "black")+
  theme(legend.position = "none",
        plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 26, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())
tjp1_box

tjp1_rain_cloud <- ggplot(dat, mapping=aes(x = sample_site, y = TJP1_mean, group = sample_site, fill = sample_site))+
  ggdist::stat_halfeye(adjust = 0.8, width = 0.6, .width = 0, justification = -0.2,point_colour = NA)+
  geom_boxplot(width = 0.15,outlier.shape = NA)+
  gghalves::geom_half_point(side = "l", range_scale = .4, alpha = .2, aes(color = sample_site) ) +
  scale_y_log10(breaks= c(1,10,100,1000,10000,100000),
                labels=c("1","10","100","1,000","10,000","100,000"))+
  coord_flip(xlim = c(1,3.1))+theme_minimal()+scale_fill_manual(values = sample_site_palette)+
  scale_color_manual(values = sample_site_palette)+
  ylab("Copies/ 20 ng of RNA")+  xlab("Sample Type") +labs(title = "TJP1")+ 
  scale_x_discrete(limits = c("LITI", "SITI","RUTI"))+
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
tjp1_rain_cloud

#test LA y or n
kruskal.test(dat$TJP1_mean, dat$LA_Y_N) # no difference p = 0.8504

#test animal
kruskal.test(dat$TJP1_mean, dat$animal)# no difference p = 0.6147

#test plate
kruskal.test(dat$TJP1_mean, dat$TJP1_plate_number)# no difference p = 007643

########CDH1
#test for normality
shapiro.test(dat$CDH1_mean)# not normal

#log transform
dat$CDH1_log <- log(dat$CDH1_mean)

#re-test
shapiro.test(dat$CDH1_log)# not normal

#krsuskal wallace
kruskal.test(dat$CDH1_mean,dat$sample_site)#p = 0.01178

#pairwise
pairwise.wilcox.test(dat$CDH1_mean, dat$sample_site, p.adjust.method = "BH")# SI different from RU and LI tends to be different than RU

cdh1_box <- ggplot(dat, aes(x= sample_site, y= CDH1_mean, fill = sample_site)) +
  theme_bw() + labs(y= "Copies/20 ng RNA\n (log10)", x="", title = "CDH1") +
  scale_x_discrete(limits=c("RUTI", "SITI", "LITI"))+
  geom_boxplot(aes(color= sample_site), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(color = sample_site), size = 2, position = position_jitter(seed = 1, w= .2))+
  scale_y_log10(breaks= c(1,10,100,1000,10000,100000),
                labels=c("1","10","100","1,000","10,000","100,000"))+
  scale_fill_manual(values = sample_site_palette, labels = c("RU TI", "SI TI", "LI TI"))+
  scale_color_manual(values = sample_site_palette, labels = c("RU TI", "SI TI", "LI TI"))+
  annotate("text", x = 3.7, y = 0.1, label = expression(paste(italic("P"), " = 0.01")), size = 5, color = "black")+
  annotate("rect", xmin = 3.4, xmax = 4, ymin = 0, ymax = 0.16, color = "black", alpha = 0)+
  annotate("text", x = 1, y= 0.18, label = "n = 19", size = 5, color= "black")+
  annotate("text", x = 2, y= 0.18, label = "n = 20", size = 5, color= "black")+
  annotate("text", x = 3, y= 0.18, label = "n = 20", size = 5, color= "black")+
  annotate("text", x = 1, y= 23000, label = "b", size = 5, color= "black")+
  annotate("text", x = 2, y= 23000, label = "a", size = 5, color= "black")+
  annotate("text", x = 3, y= 23000, label = "ab", size = 5, color= "black")+
  guides(fill = guide_legend(title = "Sample Type", byrow = TRUE), color = guide_legend(title = "Sample Type", byrow = TRUE))+
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
cdh1_box

# Function to extract legend
get_legend <- function(my_plot) {
  tmp <- ggplot_gtable(ggplot_build(my_plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Extract the legend
legend <- get_legend(cdh1_box)


CDH1_rain_cloud <- ggplot(dat, mapping=aes(x = sample_site, y = CDH1_mean, group = sample_site, fill = sample_site))+
  ggdist::stat_halfeye(adjust = 0.8, width = 0.6, .width = 0, justification = -0.2,point_colour = NA)+
  geom_boxplot(width = 0.15,outlier.shape = NA)+
  gghalves::geom_half_point(side = "l", range_scale = .4, alpha = .2, aes(color = sample_site) ) +
  scale_y_log10()+
  coord_flip(xlim = c(1,3.1))+theme_minimal()+scale_fill_manual(values = sample_site_palette)+
  scale_color_manual(values = sample_site_palette)+
  ylab("Copies/ 20 ng of RNA")+  xlab("Sample Type") +labs(title = "CDH1")+ 
  scale_x_discrete(limits = c("LITI", "SITI","RUTI"))+
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
CDH1_rain_cloud

#test LA y or n
kruskal.test(dat$CDH1_mean, dat$LA_Y_N)# no difference p = 0.4645

#test animal
kruskal.test(dat$CDH1_mean, dat$animal)# no difference p = 0.3173

# test plate
kruskal.test(dat$CDH1_mean, dat$CDH1_Plate_number)# no difference p = 0.307


# Join plots --------------------------------------------------------------
#rain cloud
plot_grid(cldn1.2_rain_cloud, cldn2_rain_cloud, ocldn_rain_cloud, CDH1_rain_cloud, tjp1_rain_cloud, align = "hv", axis = "btrl" )

#boxplots
final <- plot_grid(cldn1_box, cldn2_box, ocldn_box, cdh1_box, tjp1_box, legend, rel_widths = c(1,1,1,1,1,1), rel_heights = c(1,1,1,1,1,1),align = "hv", axis = "btrl")
final

final_long <- plot_grid(cldn1_box, cldn2_box, ocldn_box, cdh1_box, tjp1_box, legend, ncol = 2, rel_widths = c(1,1,1,1,1,1), rel_heights = c(1,1,1,1,1,1),align = "hv", axis = "btrl")
final_long


# RIN Data ----------------------------------------------------------------
#import new data sheet
RIN_Data <- read_excel("RIN_Data.xlsx")
names(RIN_Data)

#count total samples by region
RIN_Data %>% group_by(Sample_type) %>% 
  count(Sample_type)

#count NA
sum(is.na(RIN_Data$RIN))

#count NA by region
RIN_Data %>% group_by(Sample_type) %>% 
  summarise(sum_na = sum(is.na(RIN)))

#plot
hist(RIN_Data$RIN)

#average
mean(RIN_Data$RIN, na.rm = T)

#calculate average and 95% CI
RIN_avg <- RIN_Data %>% group_by(Sample_type) %>% 
  summarize(
  avg = mean(RIN, na.rm = TRUE),
  sd = sd(RIN, na.rm = TRUE),
  lower_ci = avg - qt(0.975, df = n() - 1) * sd / sqrt(n()),
  upper_ci = avg + qt(0.975, df = n() - 1) * sd / sqrt(n())
)

#test for normality
shapiro.test(RIN_Data$RIN) #not normal

#test for a difference
kruskal.test(RIN_Data$RIN, RIN_Data$Sample_type) #p = 0.2722



# Calprotectin qPCR REsults -----------------------------------------------

#read in new data set
dat2 <- read_excel("Extracted_qpcr_Data_Calpro.xlsx")
View(dat2)

#count total samples by region
dat2 %>% group_by(sample_site) %>% 
  count(sample_site)

#count NA
sum(is.na(dat2$A8_mean))

#####A8
histogram(dat2$A8_mean)

#count NA by region
dat2 %>% group_by(sample_site) %>% 
  summarise(sum_na = sum(is.na(A8_mean)))
#test for normality
shapiro.test(dat2$A8_mean)#not normal

#kruskal wallace
kruskal.test(dat2$A8_mean, dat2$sample_site)

#pairwise
pairwise.wilcox.test(dat2$A8_mean, dat2$sample_site, p.adjust.method = "BH") # RU different than LI

ggplot(dat2, aes(x= sample_site, y= A8_mean, fill = sample_site)) +
  theme_bw() + labs(y= "Copies/20 ng RNA", x="", title = "A8 Mean") +
  scale_x_discrete(limits=c("RUTI", "SITI", "LITI"))+
  geom_boxplot() +
  geom_point()+
  scale_y_log10()+
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

######A9

histogram(dat2$A9_mean)

#count NA by region
dat2 %>% group_by(sample_site) %>% 
  summarise(sum_na = sum(is.na(A9_mean)))
#test for normality
shapiro.test(dat2$A9_mean)#not normal

#kruskal wallace
kruskal.test(dat2$A9_mean, dat2$sample_site)# p = 0.00452

#pairwise
pairwise.wilcox.test(dat2$A9_mean, dat2$sample_site, p.adjust.method = "BH") # RU different than Intestines

ggplot(dat2, aes(x= sample_site, y= A9_mean, fill = sample_site)) +
  theme_bw() + labs(y= "Copies/20 ng RNA", x="", title = "A9 Mean") +
  scale_x_discrete(limits=c("RUTI", "SITI", "LITI"))+
  geom_boxplot() +
  geom_point()+
  scale_y_log10()+
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

####HG A8
histogram(dat2$HG_A8_mean_CT)

#count NA by region
dat2 %>% group_by(sample_site) %>% 
  summarise(sum_na = sum(is.na(HG_A8_mean_CT)))
#test for normality
shapiro.test(dat2$HG_A8_mean_CT)#is normal but procede non parametrically to match everything else

#kruskal wallace
kruskal.test(dat2$HG_A8_mean_CT, dat2$sample_site) #p = 0.066


ggplot(dat2, aes(x= sample_site, y= HG_A8_mean_CT, fill = sample_site)) +
  theme_bw() + labs(y= "Mean CT", x="", title = "HG A8 Mean") +
  scale_x_discrete(limits=c("RUTI", "SITI", "LITI"))+
  geom_boxplot() +
  geom_point()+
  scale_y_log10()+
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
