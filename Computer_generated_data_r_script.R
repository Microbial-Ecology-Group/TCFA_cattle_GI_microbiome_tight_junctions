
# Start Here --------------------------------------------------------------
#Project Name: TCFA 
#Project Author: JDY
#Start Date: 10Oct24
#Script Data: Computer generated H scores from IHC of Intestinal samples
#General Info: This data set contains all results from computer generated H scores of bovine IHC samples. 5 different targets were used. The first step for each target will be to join the computer generated data with human generated data to allow for statistical analysis
#Note: to use the get_legend function, run the plot one time with the legend showing, run get legend, then re-run the plot with legend.positon = "none"

# Library -----------------------------------------------------------------
setwd("C:/Users/danie/OneDrive - West Texas A and M University/Research/TCFA Challenge/IHC Scoring/Computer Assisted Scoring/R Things/")

library(readxl)
library(dplyr)
library(tidyr)
library(janitor)
library(ggplot2)
library(cowplot)

# Homemade functions ------------------------------------------------------

# Function to extract legend
get_legend <- function(my_plot) {
  tmp <- ggplot_gtable(ggplot_build(my_plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Color palettes ----------------------------------------------------------
sample_site_palette <- c( "red4", "gold","navy")

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

#look at data
hist(cldn1_data$H_negative)
hist(cldn1_data$H_pixelwise)

#test normality
shapiro.test(cldn1_data$H_pixelwise) #p = 0.07392 is normal but will use non parametric to match the others

#test for difference
kruskal.test(cldn1_data$H_pixelwise, cldn1_data$Tissue) #p = 0.00006033
pairwise.wilcox.test(cldn1_data$H_pixelwise, cldn1_data$Tissue, p.adjust.method = "BH") # all different

#set order for plot
cldn1_data$Tissue <- factor(cldn1_data$Tissue, levels = c("RU", "SI", "LI"))

#plot
cldn1_boxplot <- ggplot(cldn1_data, aes(x= Tissue, y= H_pixelwise, fill = Tissue)) +
  theme_bw() + labs(y= "Pixelwise H Score", x="", title = "Claudin 1") +
  scale_x_discrete(limits=c("RU", "SI", "LI"))+
  geom_boxplot(aes(color= Tissue), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(color = Tissue), size = 2, position = position_jitter(seed = 1, w= .2))+
  scale_fill_manual(values = sample_site_palette, label= c("RU TI", "SI TI", "LI TI"))+
  scale_color_manual(values = sample_site_palette,label= c("RU TI", "SI TI", "LI TI"))+
  coord_cartesian(y = c(0,300), x = c(1,3.2))+
  guides(fill = guide_legend(title = "Sample Type", byrow = TRUE), color = guide_legend(title = "Sample Type", byrow = TRUE))+
  annotate("text", x = 3.4, y = -5, label = expression(paste(italic("P"), " < 0.01")), size = 5, color = "black")+
  annotate("rect", xmin = 3, xmax = 4, ymin = -50, ymax = 5, color = "black", alpha = 0)+
  annotate("text", x = 1, y= 200, label = "a", size = 5, color= "black")+
  annotate("text", x = 2, y= 200, label = "b", size = 5, color= "black")+
  annotate("text", x = 3, y= 200, label = "b", size = 5, color= "black")+
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
cldn1_boxplot


#get legend to merge at the end 
legend <- get_legend(cldn1_boxplot)

#calculate avg and 95% CI

cldn1_avg <- cldn1_data %>% group_by(Tissue) %>% 
  summarise(avg = mean(H_pixelwise, na.rm = T),
            sd = sd(H_pixelwise, na.rm = T),
            lower_ci = avg - qt(0.975, df = n() - 1) * sd / sqrt(n()),
            upper_ci = avg + qt(0.975, df = n() - 1) * sd / sqrt(n()))


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

#look at data
hist(cldn2_data$H_negative)
hist(cldn2_data$H_pixelwise)

#test normality
shapiro.test(cldn2_data$H_pixelwise) #p = 0.201 is normal but will use non parametric to match the others

#test for difference
kruskal.test(cldn2_data$H_pixelwise, cldn2_data$Tissue) #p = 0.9739

#set order for plot
cldn2_data$Tissue <- factor(cldn2_data$Tissue, levels = c("RU", "SI", "LI"))

#plot
cldn2_boxplot <- ggplot(cldn2_data, aes(x= Tissue, y= H_pixelwise, fill = Tissue)) +
  theme_bw() + labs(y= "Pixelwise H Score", x="", title = "Claudin 2") +
  scale_x_discrete(limits=c("RU", "SI", "LI"))+
  geom_boxplot(aes(color= Tissue), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(color = Tissue), size = 2, position = position_jitter(seed = 1, w= .2))+
  scale_fill_manual(values = sample_site_palette, label= c("RUTI", "SITI", "LITI"))+
  scale_color_manual(values = sample_site_palette,label= c("RUTI", "SITI", "LITI"))+
  coord_cartesian(y = c(0,300), x = c(1,3.2))+
  guides(fill = guide_legend(title = "Sample Type", byrow = TRUE), color = guide_legend(title = "Sample Type", byrow = TRUE))+
  annotate("text", x = 3.4, y = -5, label = expression(paste(italic("P"), " = 0.97")), size = 5, color = "black")+
  annotate("rect", xmin = 3, xmax = 4, ymin = -50, ymax = 5, color = "black", alpha = 0)+
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
cldn2_boxplot

#calculate avg and 95% CI

cldn2_avg <- cldn2_data %>% group_by(Tissue) %>% 
  summarise(avg = mean(H_pixelwise, na.rm = T),
            sd = sd(H_pixelwise, na.rm = T),
            lower_ci = avg - qt(0.975, df = n() - 1) * sd / sqrt(n()),
            upper_ci = avg + qt(0.975, df = n() - 1) * sd / sqrt(n()))

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

#look at data
hist(ecad_data$H_negative)
hist(ecad_data$H_pixelwise)

#test normality
shapiro.test(ecad_data$H_pixelwise) #p = 0.5715 normal

#test for difference
kruskal.test(ecad_data$H_pixelwise, ecad_data$Tissue) #p = 0.001062
pairwise.wilcox.test(ecad_data$H_pixelwise, ecad_data$Tissue, p.adjust.method = "BH") #rumen different than all si li not different

#set order for plot
ecad_data$Tissue <- factor(ecad_data$Tissue, levels = c("RU", "SI", "LI"))

#plot
ecad_boxplot <- ggplot(ecad_data, aes(x= Tissue, y= H_pixelwise, fill = Tissue)) +
  theme_bw() + labs(y= "Pixelwise H Score", x="", title = "E-Cadherin") +
  scale_x_discrete(limits=c("RU", "SI", "LI"))+
  geom_boxplot(aes(color= Tissue), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(color = Tissue), size = 2, position = position_jitter(seed = 1, w= .2))+
  scale_fill_manual(values = sample_site_palette, label= c("RUTI", "SITI", "LITI"))+
  scale_color_manual(values = sample_site_palette,label= c("RUTI", "SITI", "LITI"))+
  coord_cartesian(y = c(0,300), x = c(1,3.2))+
  guides(fill = guide_legend(title = "Sample Type", byrow = TRUE), color = guide_legend(title = "Sample Type", byrow = TRUE))+
  annotate("text", x = 3.4, y = -5, label = expression(paste(italic("P"), " < 0.01")), size = 5, color = "black")+
  annotate("rect", xmin = 3, xmax = 4, ymin = -50, ymax = 5, color = "black", alpha = 0)+
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
ecad_boxplot

#calculate avg and 95% CI

ecad_avg <- ecad_data %>% group_by(Tissue) %>% 
  summarise(avg = mean(H_pixelwise, na.rm = T),
            sd = sd(H_pixelwise, na.rm = T),
            lower_ci = avg - qt(0.975, df = n() - 1) * sd / sqrt(n()),
            upper_ci = avg + qt(0.975, df = n() - 1) * sd / sqrt(n()))

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

#look at data
hist(ocldn_data$H_negative)
hist(ocldn_data$H_pixelwise)

#test normality
shapiro.test(ocldn_data$H_pixelwise) #p = 0.7478 is normal but will use non parametric to match the others

#test for difference
kruskal.test(ocldn_data$H_pixelwise, ocldn_data$Tissue) #p = 0.0000000413
pairwise.wilcox.test(ocldn_data$H_pixelwise, ocldn_data$Tissue, p.adjust.method = "BH")# all different

#set order for plot
ocldn_data$Tissue <- factor(ocldn_data$Tissue, levels = c("RU", "SI", "LI"))

#plot
ocldn_boxplot <- ggplot(ocldn_data, aes(x= Tissue, y= H_pixelwise, fill = Tissue)) +
  theme_bw() + labs(y= "Pixelwise H Score", x="", title = "Occludin") +
  scale_x_discrete(limits=c("RU", "SI", "LI"))+
  geom_boxplot(aes(color= Tissue), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(color = Tissue), size = 2, position = position_jitter(seed = 1, w= .2))+
  scale_fill_manual(values = sample_site_palette, label= c("RUTI", "SITI", "LITI"))+
  scale_color_manual(values = sample_site_palette,label= c("RUTI", "SITI", "LITI"))+
  coord_cartesian(y = c(0,300), x = c(1,3.2))+
  guides(fill = guide_legend(title = "Sample Type", byrow = TRUE), color = guide_legend(title = "Sample Type", byrow = TRUE))+
  annotate("text", x = 3.4, y = -5, label = expression(paste(italic("P"), " < 0.01")), size = 5, color = "black")+
  annotate("rect", xmin = 3, xmax = 4, ymin = -50, ymax = 5, color = "black", alpha = 0)+
  annotate("text", x = 1, y= 200, label = "a", size = 5, color= "black")+
  annotate("text", x = 2, y= 200, label = "c", size = 5, color= "black")+
  annotate("text", x = 3, y= 200, label = "b", size = 5, color= "black")+
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
ocldn_boxplot

#calculate avg and 95%CI
ocldn_avg <- ocldn_data %>% group_by(Tissue) %>% 
  summarise(avg = mean(H_pixelwise, na.rm = T),
            sd = sd(H_pixelwise, na.rm = T),
            lower_ci = avg - qt(0.975, df = n() - 1) * sd / sqrt(n()),
            upper_ci = avg + qt(0.975, df = n() - 1) * sd / sqrt(n()))


# MAC 387 -----------------------------------------------------------------
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

#look at data
hist(mac_data$H_negative)
hist(mac_data$H_pixelwise)

#test normality
shapiro.test(mac_data$H_pixelwise) #p = 0.0.00000004388 not normal

#test for difference
kruskal.test(mac_data$H_pixelwise, mac_data$Tissue) #p = 0.00005872
pairwise.wilcox.test(mac_data$H_pixelwise, mac_data$Tissue, p.adjust.method = "BH") #rumen diff from all si and li not diff

#set order for plot
mac_data$Tissue <- factor(mac_data$Tissue, levels = c("RU", "SI", "LI"))

#plot
mac_boxplot <- ggplot(mac_data, aes(x= Tissue, y= H_pixelwise, fill = Tissue)) +
  theme_bw() + labs(y= "Pixelwise H Score", x="", title = "MAC-387") +
  scale_x_discrete(limits=c("RU", "SI", "LI"))+
  geom_boxplot(aes(color= Tissue), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(color = Tissue), size = 2, position = position_jitter(seed = 1, w= .2))+
  scale_fill_manual(values = sample_site_palette, label= c("RUTI", "SITI", "LITI"))+
  scale_color_manual(values = sample_site_palette,label= c("RUTI", "SITI", "LITI"))+
  coord_cartesian(y = c(0,300), x = c(1,3.2))+
  guides(fill = guide_legend(title = "Sample Type", byrow = TRUE), color = guide_legend(title = "Sample Type", byrow = TRUE))+
  annotate("text", x = 3.4, y = -5, label = expression(paste(italic("P"), " < 0.01")), size = 5, color = "black")+
  annotate("rect", xmin = 3, xmax = 4, ymin = -50, ymax = 5, color = "black", alpha = 0)+
  annotate("text", x = 1, y= 250, label = "a", size = 5, color= "black")+
  annotate("text", x = 2, y= 250, label = "b", size = 5, color= "black")+
  annotate("text", x = 3, y= 250, label = "b", size = 5, color= "black")+
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
mac_boxplot

#calculate avg and 95%CI
mac_avg <- mac_data %>% group_by(Tissue) %>% 
  summarise(avg = mean(H_pixelwise, na.rm = T),
            sd = sd(H_pixelwise, na.rm = T),
            lower_ci = avg - qt(0.975, df = n() - 1) * sd / sqrt(n()),
            upper_ci = avg + qt(0.975, df = n() - 1) * sd / sqrt(n()))



# Merge Plots -------------------------------------------------------------

final_boxplot <- plot_grid(cldn1_boxplot,cldn2_boxplot, ocldn_boxplot, ecad_boxplot, mac_boxplot,legend,
                           ncol = 3, axis = "btrl", align = "hv")
final_boxplot
