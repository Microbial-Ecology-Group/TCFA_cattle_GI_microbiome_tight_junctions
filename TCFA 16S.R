# Library -----------------------------------------------------------------
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(vegan)
library(data.table)
library(metagenomeSeq)
library(stringr)
library(ggdendro)
library(cowplot)
library(pairwiseAdonis)
library(btools)
library(scales)
library(metagMisc)
library(GUniFrac)
library(viridis)
library(randomcoloR)
library(multcompView)
library(rcompanion)
library(microbiome)
library(devtools)


source("C:/Users/danie/OneDrive - West Texas A and M University/Bioninformatics/Source Scripts/changeSILVATaxaNames.R")
source("C:/Users/danie/OneDrive - West Texas A and M University/Bioninformatics/Source Scripts/g_unifrac.R")
source("C:/Users/danie/OneDrive - West Texas A and M University/Bioninformatics/Source Scripts/MergeLowAbund.R")
source("C:/Users/danie/OneDrive - West Texas A and M University/Bioninformatics/Source Scripts/uw_unifrac.R")
source("C:/Users/danie/OneDrive - West Texas A and M University/Bioninformatics/Source Scripts/w_unifrac.R")

# Example of a function that we can use (and re-use) to remove unwanted taxa
pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(allTaxa, physeq))
}

changeSILVAtaxa <- function(x) {
  # remove the D__ etc...
  tax.clean <- data.frame(row.names = row.names(x),
                          Kingdom = str_replace(x[,1], "k__",""),
                          Phylum = str_replace(x[,2], "p__",""),
                          Class = str_replace(x[,3], "c__",""),
                          Order = str_replace(x[,4], "o__",""),
                          Family = str_replace(x[,5], "f__",""),
                          Genus = str_replace(x[,6], "g__",""),
                          Species = str_replace(x[,7],"s_",""),
                          stringsAsFactors = FALSE)}


setwd("C:/Users/danie/OneDrive - West Texas A and M University/Research/TCFA Challenge/16S/R")


# Palettes ----------------------------------------------------------------

sample_type_palette <- c("green", "dodgerblue", "dodgerblue4", "darkviolet", "firebrick1", "firebrick4", "gold", "goldenrod2")
sample_site <- c("navy", "green", "gold", "darkviolet", "red4")
no_la_palette <- c("green", "dodgerblue", "dodgerblue4", "firebrick1", "firebrick4", "gold", "goldenrod2")
rumen_core_palette <- c("firebrick1","firebrick4")
si_core_palette <- c("gold", "goldenrod2")
li_core_palette <- c("dodgerblue", "dodgerblue4")
fecal_core_palette <- c("green")

combo_palette <-c("green", "dodgerblue", "dodgerblue4", "darkviolet", "firebrick1", "firebrick4", "gold", "goldenrod2", animal_palette) 

# Import Data -------------------------------------------------------------

#read in microbiome data
tcfadata <- import_biom('table-with-taxonomy.biom', 'tree.nwk', 'dna-sequences.fasta')

#Write Sample names to create metadata
#write.csv(sample_names(tcfadata), "samplenames1.csv")

#read in metadata
metadata <- import_qiime_sample_data('samplenames1.txt')

#merge metadata
data <- merge_phyloseq(tcfadata, metadata)
data #165 samples 42247 taxa


# Explore Data ------------------------------------------------------------

#names pre processing
#rank then paste better names
rank_names(data)
colnames(tax_table(data)) <- c("Kingdom","Phylum","Class","Order","Family","Genus", "Species")
rank_names(data)

#change the silva style names
tax.data <- data.frame(tax_table(data)) 
#changed this from changesilvanames, because it was reverting back to domain instead of kingdom
tax.data.names <- changeSILVAtaxa(tax.data)

#change NA's to better names
#convert to charcters
for(i in 1:7){ tax.data.names[,i] <- as.character(tax.data.names[,i])}
#replace with empty string
tax.data.names[is.na(tax.data.names)] <- ""

head(tax.data)

# now filling in the empty slots with the highest assigned taxonomy
for (i in 1:nrow(tax.data.names)){
  if (tax.data.names[i,2] == ""){
    domain <- paste("unclassified ", tax.data.names[i,1], sep = "")
    tax.data.names[i, 2:7] <- domain
  } else if (tax.data.names[i,3] == ""){
    kingdom <- paste("unclassified ", tax.data.names[i,2], sep = "")
    tax.data.names[i, 3:7] <- kingdom
  } else if (tax.data.names[i,4] == ""){
    phylum <- paste("unclassified ", tax.data.names[i,3], sep = "")
    tax.data.names[i, 4:7] <- phylum
  } else if (tax.data.names[i,5] == ""){
    class <- paste("unclassified ", tax.data.names[i,4], sep = "")
    tax.data.names[i, 5:7] <- class
  } else if (tax.data.names[i,6] == ""){
    order <- paste("unclassified ", tax.data.names[i,5], sep = "")
    tax.data.names[i, 6:7] <- order
  } else if (tax.data.names[i,7] == ""){
    tax.data.names$Genus[i] <- paste("unclassified ",tax.data.names$Family[i], sep = "")
  }
}

#check headers re-insert matrix and check tail
head(tax.data.names)
tax_table(data) <- as.matrix(tax.data.names)
tail(tax_table(data))

#prune for taxa within Eukaryota (We dont have domain with this classifier)
data <- subset_taxa(data, Domain!="Eukaryota")
data


#split samples and controls
samples <- subset_samples(data, sample_type1 == "Sample")
samples <- prune_taxa(taxa_sums(samples) > 0, samples)
samples #40969 taxa 152 samples

controls <- subset_samples(data, sample_type1 =="Control")
controls <- prune_taxa(taxa_sums(controls) > 0, controls)
controls #4532 taxa 12 sample
sort(sample_sums(controls)) #EB2 LIFlEB6   NTC15    NTC9   NTC12   NTC11    NTC3    NTC7    NTC2    NTC4    NTC6    NTC5 
#                             1       1       1       1       3      30      58      86   12305   81981  116164  507039 
controls <- prune_samples(sample_sums(controls)> 100, controls)

#investigate the ntc with reads
controls_genus <-  tax_glom(controls, taxrank = "Genus", NArm = F) %>% 
  psmelt()
length(unique(controls_genus$Genus)) #328 genera
ggplot(controls_genus, aes(x=sample_name, y=Abundance, fill= Genus))+
  geom_bar(stat = "identity")
#Need Lee to look at this.NTC5 has a lot of reads probably the one that got switched in library prep.

#look at number of reads per sample and distributions
sample_sum_df <- data.frame(sum= sample_sums(samples)) 
ggplot(sample_sum_df, aes(x= sum))+
  geom_histogram(color= "black", fill = "indianred", binwidth = 5000)+
  ggtitle("Distribution of sample sequencing depth")+ 
  xlab("Read Counts")+
  theme(axis.title.y = element_blank())

mean(sample_sums(samples)) #561,001.9
median(sample_sums(samples)) #574,024.5
sort(sample_sums(samples)) 

#cut off below 250,000 Need to talk with morley and lee about this 
samples <- prune_samples(sample_sums(samples) > 250000, samples)
samples #cuts us down to 145 samples
samples <- prune_taxa(taxa_sums(samples) > 0, samples)
samples #40661 taxa and 145 samples remain
min(sample_sums(samples)) #270,390
max(sample_sums(samples))#886,923
mean(sample_sums(samples))#581,594.4
sort(sample_sums(samples))
sum(sample_sums(samples))

# Comparing Between Samples -----------------------------------------------
sample_sum_df_new <- data.frame(ASV_count = sample_sums(samples))
metadata.df <- as(sample_data(samples), "data.frame") 
seqDepth_metadata <- cbind(metadata.df, sample_sum_df_new)

#comparing livers
kruskal.test(seqDepth_metadata$ASV_count, seqDepth_metadata$liver_score) # 0.04

#comparing body site
kruskal.test(seqDepth_metadata$ASV_count, seqDepth_metadata$body_site) #<0.001
pairwise.wilcox.test(seqDepth_metadata$ASV_count, seqDepth_metadata$body_site)

ggplot(seqDepth_metadata, aes(x= body_site, y = ASV_count))+
  geom_boxplot()

#comparing for animal
kruskal.test(seqDepth_metadata$ASV_count, seqDepth_metadata$Animal)#NS

#comparing sample type (FL vs Epi and by region)
kruskal.test(seqDepth_metadata$ASV_count, seqDepth_metadata$sample_type)# < 0.001

# Checking classification efficiency --------------------------------------

# First we have to extract the tax_table and convert it to a "data.frame"
taxa.df <- as.data.frame(tax_table(samples))
taxa.df
#KINGDOM
# Now let's search for which taxa start with "unclassified"
unclassified_kingdom.df <- taxa.df %>% filter(grepl('Unassigned', Kingdom))
# Pull out the unique 
unclassified_kingdom <- row.names(unclassified_kingdom.df)

# We can run the "pop_taxa" function here using the genus phyloseq object and the list of features we want to remove
trimmed_microbiome_kingdom.ps = pop_taxa(samples, unclassified_kingdom)

# Sum the counts for each sample
sample_sums(trimmed_microbiome_kingdom.ps)

# To sum across all samples, we can place that command inside the "sum()" function
sum(sample_sums(trimmed_microbiome_kingdom.ps))

# Now, we calculate the proportion of reads mapped to the species level, out of all microbiome mapped reads
sum(sample_sums(trimmed_microbiome_kingdom.ps)) / sum(sample_sums(samples)) * 100

###Kingdom 99.90254%

#Phylum
# Now let's search for which taxa start with "unclassified"
unclassified_phy.df <- taxa.df %>% filter(grepl('unclassified', Phylum))
# Pull out the unique 
unclassified_phy <- row.names(unclassified_phy.df)

# We can run the "pop_taxa" function here using the genus phyloseq object and the list of features we want to remove
trimmed_microbiome_phy.ps = pop_taxa(samples, unclassified_phy)

# Sum the counts for each sample
sample_sums(trimmed_microbiome_phy.ps)

# To sum across all samples, we can place that command inside the "sum()" function
sum(sample_sums(trimmed_microbiome_phy.ps))

# Now, we calculate the proportion of reads mapped to the species level, out of all microbiome mapped reads
sum(sample_sums(trimmed_microbiome_phy.ps)) / sum(sample_sums(samples)) * 100

###Phylum 99.86914%

#CLASS
# Now let's search for which taxa start with "unclassified"
unclassified_class.df <- taxa.df %>% filter(grepl('unclassified', Class))
# Pull out the unique 
unclassified_class <- row.names(unclassified_class.df)

# We can run the "pop_taxa" function here using the genus phyloseq object and the list of features we want to remove
trimmed_microbiome_class.ps = pop_taxa(samples, unclassified_class)

# Sum the counts for each sample
sample_sums(trimmed_microbiome_class.ps)

# To sum across all samples, we can place that command inside the "sum()" function
sum(sample_sums(trimmed_microbiome_class.ps))

# Now, we calculate the proportion of reads mapped to the species level, out of all microbiome mapped reads
sum(sample_sums(trimmed_microbiome_class.ps)) / sum(sample_sums(samples)) * 100

###Class 99.86158%

#ORDER
# Now let's search for which taxa start with "unclassified"
unclassified_order.df <- taxa.df %>% filter(grepl('unclassified', Order))
# Pull out the unique 
unclassified_order <- row.names(unclassified_order.df)

# We can run the "pop_taxa" function here using the genus phyloseq object and the list of features we want to remove
trimmed_microbiome_order.ps = pop_taxa(samples, unclassified_order)

# Sum the counts for each sample
sample_sums(trimmed_microbiome_order.ps)

# To sum across all samples, we can place that command inside the "sum()" function
sum(sample_sums(trimmed_microbiome_order.ps))

# Now, we calculate the proportion of reads mapped to the species level, out of all microbiome mapped reads
sum(sample_sums(trimmed_microbiome_order.ps)) / sum(sample_sums(samples)) * 100

###Order 99.82529

#FAMILY
# Now let's search for which taxa start with "unclassified"
unclassified_family.df <- taxa.df %>% filter(grepl('unclassified', Family))
# Pull out the unique 
unclassified_family <- row.names(unclassified_family.df)

# We can run the "pop_taxa" function here using the genus phyloseq object and the list of features we want to remove
trimmed_microbiome_family.ps = pop_taxa(samples, unclassified_family)

# Sum the counts for each sample
sample_sums(trimmed_microbiome_family.ps)

# To sum across all samples, we can place that command inside the "sum()" function
sum(sample_sums(trimmed_microbiome_family.ps))

# Now, we calculate the proportion of reads mapped to the species level, out of all microbiome mapped reads
sum(sample_sums(trimmed_microbiome_family.ps)) / sum(sample_sums(samples)) * 100

###Family 99.6918%

#GENUS
# Now let's search for which taxa start with "unclassified"
unclassified_genus.df <- taxa.df %>% filter(grepl('unclassified', Genus))
# Pull out the unique 
unclassified_genus <- row.names(unclassified_genus.df)

# We can run the "pop_taxa" function here using the genus phyloseq object and the list of features we want to remove
trimmed_microbiome_genus.ps = pop_taxa(samples, unclassified_genus)

# Sum the counts for each sample
sample_sums(trimmed_microbiome_genus.ps)

# To sum across all samples, we can place that command inside the "sum()" function
sum(sample_sums(trimmed_microbiome_genus.ps))

# Now, we calculate the proportion of reads mapped to the species level, out of all microbiome mapped reads
sum(sample_sums(trimmed_microbiome_genus.ps)) / sum(sample_sums(samples)) * 100

###Genus 63.65056% ??? Need to ask about this

# Alpha Diversity ---------------------------------------------------------

alphadiv1 <- estimate_richness(samples, measures = c("Observed", "Shannon", "Simpson", "InvSimpson"))
alphadiv2 <- estimate_pd(samples)

alpha_div <- cbind(alphadiv1, alphadiv2)
alpha_div

alpha_div <- alpha_div[,c(1:5)]
alpha_div
alpha_div_meta <- cbind(metadata.df,alpha_div)
alpha_div_meta

#pairwise for observed
obs_anova <- pairwise.wilcox.test(alpha_div_meta$Observed, alpha_div_meta$sample_type, p.adjust.method = "BH")
#superscript p-values
obs_pvalues <- obs_anova$p.value %>% 
  fullPTable()
multcompLetters(obs_pvalues)

#pairwise for shannon
shan_anova <- pairwise.wilcox.test(alpha_div_meta$Shannon, alpha_div_meta$sample_type, p.adjust.method = "BH")

#superscript p values
shan_pvalues <- shan_anova$p.value %>% 
  fullPTable()
multcompLetters(shan_pvalues)

#plots
microbiome_richness_boxplot <- ggplot(alpha_div_meta, aes(x= sample_type, y= Observed, fill = sample_type)) +
  theme_bw() + labs(y= "", x= "", title = "") +
  scale_x_discrete(limits=c("RUFL", "RUTI", "SIFL", "SITI", "LIFL", "LITI", "Fecal", "Liver_Abscess"))+
  geom_boxplot() +
  geom_point()+
  scale_fill_manual(values = sample_type_palette)+
  theme(#legend.position = "none",
        plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 30, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())

microbiome_richness_boxplot

microbiome_shannon_boxplot <- ggplot(alpha_div_meta, aes(x= sample_type, y= Shannon, fill = sample_type)) +
  theme_bw() + labs(y= "", x="", title = "") +
  scale_x_discrete(limits=c("RUFL", "RUTI", "SIFL", "SITI", "LIFL", "LITI", "Fecal", "Liver_Abscess"))+
  scale_y_continuous(limits = c(2,6.5))+
  scale_fill_manual(values= sample_typle_palette)+
  geom_boxplot() +
  geom_point()+
  theme(#legend.position = "none",
        plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 30, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())

microbiome_shannon_boxplot

plot_grid(microbiome_richness_boxplot,microbiome_shannon_boxplot, align = "h", axis = "rl", ncol = 2)

microbiome_PD_boxplot <- ggplot(alpha_div_meta, aes(x= sample_type, y= PD, fill = sample_type)) +
  theme_bw() + labs(y= "PD", x="", title = "Faiths Distance") +
  scale_x_discrete(limits=c("RUFL", "RUTI", "SIFL", "SITI", "LIFL", "LITI", "Fecal", "Liver_Abscess"))+
  geom_boxplot() +
  geom_point()+ 
  scale_fill_manual(values = sample_type_palette)+
  theme(#legend.position = "none",
        panel.border = element_rect(colour = "black", size = 1),
        title = element_text(size = 28),
        axis.ticks = element_line(colour = "black", size = 0.75),
        axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(size = 24),
        axis.text.x = element_blank())
microbiome_PD_boxplot

# Beta Diversity ----------------------------------------------------------

####CSS Transform###
any(taxa_sums(samples)==0)
samples.css <- phyloseq_transform_css(samples, log = F)
samples.css.df <- as(sample_data(samples.css), "data.frame")

###split by body part
#RUFL
RUFL.css <- subset_samples(samples.css, sample_type == "RUFL",)
RUFL.css <- prune_taxa(taxa_sums(RUFL.css)>0,RUFL.css)
RUFL.css #21 samples 12930 taxa
RUFL.css.df <- as(sample_data(RUFL.css), "data.frame")

#RUTI
RUTI.css <- subset_samples(samples.css, sample_type == "RUTI",)
RUTI.css <- prune_taxa(taxa_sums(RUTI.css)>0,RUTI.css)
RUTI.css #19 samples 12956 taxa
RUTI.css.df <- as(sample_data(RUTI.css), "data.frame")

#SIFL
SIFL.css <- subset_samples(samples.css, sample_type == "SIFL",)
SIFL.css <- prune_taxa(taxa_sums(SIFL.css)>0,SIFL.css)
SIFL.css #16 samples 6060 taxa
SIFL.css.df <- as(sample_data(SIFL.css), "data.frame")

#SITI
SITI.css <- subset_samples(samples.css, sample_type == "SITI",)
SITI.css <- prune_taxa(taxa_sums(SITI.css)>0,SITI.css)
SITI.css #19 samples 14528 taxa
SITI.css.df <- as(sample_data(SITI.css), "data.frame")


#LIFL
LIFL.css <- subset_samples(samples.css, sample_type == "LIFL",)
LIFL.css <- prune_taxa(taxa_sums(LIFL.css)>0,LIFL.css)
LIFL.css #20 samples 15406 taxa
LIFL.css.df <- as(sample_data(LIFL.css), "data.frame")

#LITI
LITI.css <- subset_samples(samples.css, sample_type == "LITI",)
LITI.css <- prune_taxa(taxa_sums(LITI.css)>0,LITI.css)
LITI.css #21 samples 14428 taxa
LITI.css.df <- as(sample_data(LITI.css), "data.frame")

#Fecal
fecal.css <- subset_samples(samples.css, sample_type == "Fecal")
fecal.css <- prune_taxa(taxa_sums(fecal.css)>0,fecal.css)
fecal.css #19 samples 14436 taxa

#LA
LA.css <- subset_samples(samples.css, sample_type == "Liver_Abscess",)
LA.css <- prune_taxa(taxa_sums(LA.css)>0,LA.css)
LA.css #10 samples 136 taxa
LA.css.df <- as(sample_data(LA.css), "data.frame")

# Unifrac and Ordination by body site -------------------------------------

samples.dist <- gunifrac(samples.css)
samples.ord <- ordinate(samples.css, method = "NMDS", distance = samples.dist)

#Ordination plots
#sample type
plot_ordination(samples.css, samples.ord, type = "samples", color = "sample_type")+
  theme_bw()+
  labs(title = "", x= "NMDS1", y= "NMDS2")+
  geom_point(size= 5, shape= 18)+
  #geom_text(label = samples.css@sam_data$Sample_names)+
  stat_ellipse(geom = "polygon", aes(fill= sample_type), alpha = 0.1, lty= 1, size= 1)+
  scale_colour_manual(values = sample_type_palette)+
  scale_fill_manual(values = sample_type_palette)+
  #(limits=c("Oral", "Rumen", "Abomasum", "Duodenum", "Jejunum", "Ileum", "Cecum", "Spiral_Colon", "Distal_Colon","Fecal"))+
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", size = 1),
        title = element_text(size = 28),
        axis.ticks = element_line(colour = "black", size = 0.75),
        axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(size = 24))

#ordinate by sample type without LA
no_la <- subset_samples(samples.css, sample_type %in% c("RUTI", "RUFL", "SITI","SIFL","LITI", "LIFL", "Fecal"))
no_la_samples.dist <- gunifrac(no_la)
no_la.ord <- ordinate(no_la, method = "NMDS", distance = no_la_samples.dist)
no_la
samples.css

plot_ordination(no_na_physeq, no_la.ord, type = "samples", color = "sample_type")+
  theme_bw()+
  labs(title = "", x= "NMDS1", y= "NMDS2")+
  geom_point(size= 5, shape= 18)+
  #geom_text(label = samples.css@sam_data$Sample_names)+
  stat_ellipse(geom = "polygon", aes(fill= sample_type), alpha = 0.1, lty= 1, size= 1)+
  scale_colour_manual(values = no_la_palette)+
  scale_fill_manual(values = no_la_palette)+
  theme(legend.position = "none",
    panel.border = element_rect(colour = "black", size = 1),
    title = element_text(size = 28),
    axis.ticks = element_line(colour = "black", size = 0.75),
    axis.text = element_text(colour = "black", size = 12),
    axis.title = element_text(size = 24))


#by animal
animal_palette <- randomColor(21)

plot_ordination(samples.css, samples.ord, type = "samples", color = "Animal")+
  theme_bw()+
  labs(title = "", x= "NMDS1", y= "NMDS2")+
  geom_point(size= 5, shape= 18)+
  #geom_text(label = samples.css@sam_data$Sample_names)+
  stat_ellipse(geom = "polygon", aes(fill= Animal), alpha = .1, lty= 1, size= 1)+
  scale_colour_manual(values = animal_palette)+
  scale_fill_manual(values = animal_palette)+
  theme(legend.position = "none",
    panel.border = element_rect(colour = "black", size = 1),
    title = element_text(size = 28),
    axis.ticks = element_line(colour = "black", size = 0.75),
    axis.text = element_text(colour = "black", size = 12),
    axis.title = element_text(size = 24))

#body site
plot_ordination(samples.css, samples.ord, type = "samples", color = "body_site")+
  theme_bw()+
  labs(title = "", x= "NMDS1", y= "NMDS2")+
  geom_point(size= 5, shape= 18)+
  #geom_text(label = samples.css@sam_data$Sample_names)+
  stat_ellipse(geom = "polygon", aes(fill= body_site), alpha = 0.1, lty= 1, size= 1)+
  scale_colour_manual(values = sample_site)+
  scale_fill_manual(values = sample_site)+
  #(limits=c("Oral", "Rumen", "Abomasum", "Duodenum", "Jejunum", "Ileum", "Cecum", "Spiral_Colon", "Distal_Colon","Fecal"))+
  theme(legend.position = "none",
    panel.border = element_rect(colour = "black", size = 1),
    title = element_text(size = 28),
    axis.ticks = element_line(colour = "black", size = 0.75),
    axis.text = element_text(colour = "black", size = 12),
    axis.title = element_text(size = 24))


# Relative Abundance ------------------------------------------------------
#calculations
rel_abund <- transform_sample_counts(samples.css, function(x) {x/sum(x)} * 100)
ra_family <- tax_glom(rel_abund, taxrank = "Family", NArm = F)
ra_family_filt <- merge_low_abundance(ra_family, threshold = 0.001) 
ra_family_filt #145 samples 165 taxa



#test for a difference in sample type
microbiome_df = as.data.frame(as(sample_data(samples.css),"matrix"))
adonis_microbiome <- adonis2(samples.dist ~sample_type, microbiome_df)
adonis_microbiome

permutest(betadisper(samples.dist, microbiome_df$sample_type), pairwise = T)

#pairwise
adonis_microbiome_body_site <- pairwise.adonis2(samples.dist ~ sample_type, microbiome_df)
# view the results
adonis_microbiome_body_site

#by animal 
animal <- adonis2(samples.dist ~Animal, microbiome_df)
animal
pairwise.adonis2(samples.dist~Animal, microbiome_df)
permutest(betadisper(samples.dist, microbiome_df$Animal), pairwise = T)


#check dispersion
permutest(betadisper(samples.dist, microbiome_df$sample_type), pairwise = TRUE)

# Check individual sample RA for odd things by body site ------------------

#RUFL
ra_RUFL <- subset_samples(rel_abund, sample_type =="RUFL")
ra_RUFL <- prune_taxa(taxa_sums(ra_RUFL) > 0, ra_RUFL)
ra_RUFL #12930 taxa, 21 samples
RUFL_phy <- tax_glom(ra_RUFL, taxrank = "Phylum", NArm = F)
RUFL_class <- tax_glom(ra_RUFL, taxrank = "Class", NArm = F)
RUFL_family <- tax_glom(ra_RUFL, taxrank = "Family", NArm = F)
RUFL_genus <- tax_glom(ra_RUFL, taxrank = "Genus", NArm = F)

RUFL_family_filt <- merge_low_abundance(RUFL_family, threshold = 0.1)
RUFL_family_filt  #39 taxa
RUFL_family_melt <- psmelt(RUFL_family_filt)

RUFL_genus_melt <- merge_low_abundance(RUFL_genus, threshold = 0.1) %>% 
  psmelt()


RUFLPlot <- ggplot(RUFL_family_melt, aes(x= Sample, y= Abundance, fill = Family)) +
  theme_minimal() + #facet_wrap(~Body_Site, scales = "free_x")+
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_fill_manual(values = randomColor(89) ) +
  theme(#legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.y = element_line(size = 0.75, colour = "black"),
    axis.title.y = element_text(size = 24),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_line(colour = "black", size = 0.7))
RUFLPlot

#RUTI
ra_RUTI <- subset_samples(rel_abund, sample_type =="RUTI")
ra_RUTI <- prune_taxa(taxa_sums(ra_RUTI) > 0, ra_RUTI)
ra_RUTI #12956 taxa, 19 samples
RUTI_phy <- tax_glom(ra_RUTI, taxrank = "Phylum", NArm = F)
RUTI_class <- tax_glom(ra_RUTI, taxrank = "Class", NArm = F)
RUTI_family <- tax_glom(ra_RUTI, taxrank = "Family", NArm = F)
RUTI_genus <- tax_glom(ra_RUTI, taxrank = "Genus", NArm = F)

RUTI_family_filt <- merge_low_abundance(RUTI_family, threshold = 0.1)
RUTI_family_filt  #49 taxa
RUTI_family_melt <- psmelt(RUTI_family_filt)

RUTI_genus_melt <- merge_low_abundance(RUTI_genus, threshold = 0.1) %>% 
  psmelt()

RUTIPlot <- ggplot(RUTI_family_melt, aes(x= Sample, y= Abundance, fill = Family)) +
  theme_minimal() + #facet_wrap(~Body_Site, scales = "free_x")+
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_fill_manual(values = randomColor(99) ) +
  theme(#legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.y = element_line(size = 0.75, colour = "black"),
    axis.title.y = element_text(size = 24),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_line(colour = "black", size = 0.7))
RUTIPlot


#SIFL
ra_SIFL <- subset_samples(rel_abund, sample_type =="SIFL")
ra_SIFL <- prune_taxa(taxa_sums(ra_SIFL) > 0, ra_SIFL)
ra_SIFL #6060 taxa, 16 samples
SIFL_phy <- tax_glom(ra_SIFL, taxrank = "Phylum", NArm = F)
SIFL_class <- tax_glom(ra_SIFL, taxrank = "Class", NArm = F)
SIFL_family <- tax_glom(ra_SIFL, taxrank = "Family", NArm = F)
SIFL_genus <- tax_glom(ra_SIFL, taxrank = "Genus", NArm = F)

SIFL_family_filt <- merge_low_abundance(SIFL_family, threshold = 0.1)
SIFL_family_filt  #32 taxa
SIFL_family_melt <- psmelt(SIFL_family_filt)

SIFL_genus_melt <- merge_low_abundance(SIFL_genus, threshold = 0.1) %>% 
  psmelt()


SIFLPlot <- ggplot(SIFL_genus_melt, aes(x= Sample, y= Abundance, fill = Genus)) +
  theme_minimal() + #facet_wrap(~Body_Site, scales = "free_x")+
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_fill_manual(values = randomColor(61) ) +
  theme(#legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.y = element_line(size = 0.75, colour = "black"),
    axis.title.y = element_text(size = 24),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_line(colour = "black", size = 0.7))
SIFLPlot

#SITI
ra_SITI <- subset_samples(rel_abund, sample_type =="SITI")
ra_SITI <- prune_taxa(taxa_sums(ra_SITI) > 0, ra_SITI)
ra_SITI #14528 taxa, 19 samples
SITI_phy <- tax_glom(ra_SITI, taxrank = "Phylum", NArm = F)
SITI_class <- tax_glom(ra_SITI, taxrank = "Class", NArm = F)
SITI_family <- tax_glom(ra_SITI, taxrank = "Family", NArm = F)
SITI_genus <- tax_glom(ra_SITI, taxrank = "Genus", NArm = F)

SITI_family_filt <- merge_low_abundance(SITI_family, threshold = 0.1)
SITI_family_filt  #51 taxa
SITI_family_melt <- psmelt(SITI_family_filt)

SITI_genus_melt <- merge_low_abundance(SITI_genus, threshold = 0.1) %>% 
  psmelt()

SITIPlot <- ggplot(SITI_genus_melt, aes(x= Sample, y= Abundance, fill = Genus)) +
  theme_minimal() + #facet_wrap(~Body_Site, scales = "free_x")+
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_fill_manual(values = randomColor(105) ) +
  theme(legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.y = element_line(size = 0.75, colour = "black"),
    axis.title.y = element_text(size = 24),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_line(colour = "black", size = 0.7))
SITIPlot

#LIFL
ra_LIFL <- subset_samples(rel_abund, sample_type =="LIFL")
ra_LIFL <- prune_taxa(taxa_sums(ra_LIFL) > 0, ra_LIFL)
ra_LIFL #15406 taxa, 20 samples
LIFL_phy <- tax_glom(ra_LIFL, taxrank = "Phylum", NArm = F)
LIFL_class <- tax_glom(ra_LIFL, taxrank = "Class", NArm = F)
LIFL_family <- tax_glom(ra_LIFL, taxrank = "Family", NArm = F)
LIFL_genus <- tax_glom(ra_LIFL, taxrank = "Genus", NArm = F)

LIFL_family_filt <- merge_low_abundance(LIFL_family, threshold = 0.1)
LIFL_family_filt  #41 taxa
LIFL_family_melt <- psmelt(LIFL_family_filt)

LIFL_genus_melt <- merge_low_abundance(LIFL_genus, threshold = 0.1) %>% 
  psmelt()


LIFLPlot <- ggplot(SIFL_family_melt, aes(x= Sample, y= Abundance, fill = Family)) +
  theme_minimal() + #facet_wrap(~Body_Site, scales = "free_x")+
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_fill_manual(values = randomColor(61) ) +
  theme(#legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.y = element_line(size = 0.75, colour = "black"),
    axis.title.y = element_text(size = 24),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_line(colour = "black", size = 0.7))
LIFLPlot

#LITI
ra_LITI <- subset_samples(rel_abund, sample_type =="LITI")
ra_LITI <- prune_taxa(taxa_sums(ra_LITI) > 0, ra_LITI)
ra_LITI #14428 taxa, 19 samples
LITI_phy <- tax_glom(ra_LITI, taxrank = "Phylum", NArm = F)
LITI_class <- tax_glom(ra_LITI, taxrank = "Class", NArm = F)
LITI_family <- tax_glom(ra_LITI, taxrank = "Family", NArm = F)
LITI_genus <- tax_glom(ra_LITI, taxrank = "Genus", NArm = F)

LITI_family_filt <- merge_low_abundance(LITI_family, threshold = 0.1)
LITI_family_filt  #51 taxa
LITI_family_melt <- psmelt(LITI_family_filt)

LITI_genus_melt <- merge_low_abundance(SITI_genus, threshold = 0.1) %>% 
  psmelt()

LITIPlot <- ggplot(LITI_family_melt, aes(x= Sample, y= Abundance, fill = Family)) +
  theme_minimal() + #facet_wrap(~Body_Site, scales = "free_x")+
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_fill_manual(values = randomColor(105) ) +
  theme(#legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.y = element_line(size = 0.75, colour = "black"),
    axis.title.y = element_text(size = 24),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_line(colour = "black", size = 0.7))
LITIPlot

#Fecal
ra_fecal <-  subset_samples(rel_abund, sample_type =="Fecal")
ra_fecal <- prune_taxa(taxa_sums(ra_fecal)> 0, ra_fecal)
ra_fecal#14436 taxa 19 samples
fecal_phy <- tax_glom(ra_fecal, taxrank = "Phylum", NArm = F)
fecal_class <- tax_glom(ra_fecal, taxrank = "Class", NArm = F)
fecal_family <- tax_glom(ra_fecal, taxrank = "Family", NArm = F)
fecal_genus <- tax_glom(ra_fecal, taxrank = "Genus", NArm = F)

fecal_family_melt <- merge_low_abundance(fecal_family, threshold = 0.1) %>% 
  psmelt()

fecal_genus_melt <- merge_low_abundance(fecal_genus, threshold = 0.1) %>% 
  psmelt()

fecalPlot <- ggplot(fecal_genus_melt, aes(x= Sample, y= Abundance, fill = Genus)) +
  theme_minimal() + #facet_wrap(~Body_Site, scales = "free_x")+
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_fill_manual(values = randomColor(105) ) +
  theme(#legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.y = element_line(size = 0.75, colour = "black"),
    axis.title.y = element_text(size = 24),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_line(colour = "black", size = 0.7))
fecalPlot

#plot all fecal/colon for comparison to re-runs
ra_hindgut <-  subset_samples(rel_abund, sample_type %in% c("Fecal", "LITI", "LIFL"))
ra_ <- prune_taxa(taxa_sums(ra_hindgut)> 0, ra_hindgut)
ra_hindgut #40661 taxa and 16 60 samples
hindgut_phy <- tax_glom(ra_hindgut, taxrank = "Phylum", NArm = F)
hindgut_class <- tax_glom(ra_hindgut, taxrank = "Class", NArm = F)
hindgut_family <- tax_glom(ra_hindgut, taxrank = "Family", NArm = F)
hindgut_genus <- tax_glom(ra_hindgut, taxrank = "Genus", NArm = F)

hindgut_family_melt <- merge_low_abundance(hindgut_family, threshold = 0.1) %>% 
  psmelt()

hindgut_phy_melt <- merge_low_abundance(hindgut_phy, threshold = 0.1) %>% 
  psmelt()

hindgut_phy_Plot <- ggplot(hindgut_phy_melt, aes(x= Sample, y= Abundance, fill = Phylum)) +
  theme_minimal() + #facet_wrap(~Body_Site, scales = "free_x")+
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_fill_manual(values = randomColor(105) ) +
  theme(#legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.y = element_line(size = 0.75, colour = "black"),
    axis.title.y = element_text(size = 24),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_line(colour = "black", size = 0.7))
hindgut_phy_Plot

hindgut_fam_Plot <- ggplot(hindgut_family_melt, aes(x= Sample, y= Abundance, fill = Family)) +
  theme_minimal() + #facet_wrap(~Body_Site, scales = "free_x")+
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_fill_manual(values = randomColor(105) ) +
  theme(#legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.y = element_line(size = 0.75, colour = "black"),
    axis.title.y = element_text(size = 24),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_line(colour = "black", size = 0.7))
hindgut_fam_Plot

#LA
ra_LA <-  subset_samples(rel_abund, sample_type =="Liver_Abscess")
ra_LA <- prune_taxa(taxa_sums(ra_LA)> 0, ra_LA)
ra_LA#136 taxa 10 samples
LA_phy <- tax_glom(ra_LA, taxrank = "Phylum", NArm = F)
LA_class <- tax_glom(ra_LA, taxrank = "Class", NArm = F)
LA_family <- tax_glom(ra_LA, taxrank = "Family", NArm = F)
LA_genus <- tax_glom(ra_LA, taxrank = "Genus", NArm = F)

LA_family_melt <- merge_low_abundance(LA_family, threshold = 0.1) %>% 
  psmelt()

LA_genus_melt <- merge_low_abundance(LA_genus, threshold = 0.1) %>% 
  psmelt()

LAPlot <- ggplot(LA_family_melt, aes(x= Sample, y= Abundance, fill = Family)) +
  theme_minimal() + #facet_wrap(~Body_Site, scales = "free_x")+
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_fill_manual(values = randomColor(105) ) +
  geom_text(aes(label=sample_name), data = LA_family_melt)+
  theme(#legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.y = element_line(size = 0.75, colour = "black"),
    axis.title.y = element_text(size = 24),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_line(colour = "black", size = 0.7))
LAPlot


# RA phylum and Dendrogram ------------------------------------------------

ra_samples <- prune_taxa(taxa_sums(samples.css) > 0, samples.css)
ra_samples#40661 taxa, 145 samples
samples_phy <- tax_glom(ra_samples, taxrank = "Phylum", NArm = F)
samples_class <- tax_glom(ra_samples, taxrank = "Class", NArm = F)
Samples_order <- tax_glom(ra_samples, taxrank = "Order", NArm = F)
samples_family <- tax_glom(ra_samples, taxrank = "Family", NArm = F)
samples_genus <- tax_glom(ra_samples, taxrank = "Genus", NArm = F)


#cluster
samples.hclust <- hclust(samples.dist, method = "ward.D2")
samples.dendro.plot <- plot(samples.hclust)
samples.dendro <- as.dendrogram(samples.hclust)
samples.dendro.data <- dendro_data(samples.dendro, type = "rectangle")
metadata_for_dendro.samples <- as_tibble(samples.css@sam_data)
samples.dendro.data$labels <- samples.dendro.data$labels %>% 
  left_join(metadata_for_dendro.samples, by = c("label" = "sample_name"))
samples.dendro.data$labels
samples.dendro.order <- samples.dendro.data$labels$label
samples.dendro.order

#cluster w/ no LA
no_samples.hclust <- hclust(no_la_samples.dist, method = "ward.D2")
no_samples.dendro.plot <- plot(no_samples.hclust)
no_samples.dendro <- as.dendrogram(no_samples.hclust)
no_samples.dendro.data <- dendro_data(no_samples.dendro, type = "rectangle")
no_metadata_for_dendro.samples <- as_tibble(no_la@sam_data)
no_samples.dendro.data$labels <- no_samples.dendro.data$labels %>% 
  left_join(no_metadata_for_dendro.samples, by = c("label" = "sample_name"))
no_samples.dendro.data$labels
no_samples.dendro.order <- no_samples.dendro.data$labels$label
no_samples.dendro.order

#look at phylum
samples_phy #33 taxa 145 samples
samples_phy_filt <- merge_low_abundance(samples_phy, threshold = 0.1)
samples_phy_filt #15 taxa
samples_phy_melt <- psmelt(samples_phy_filt)

#phy no LA
no_ra_samples <- prune_taxa(taxa_sums(no_la) > 0, no_la)
no_ra_samples#40578 taxa, 135 samples
no_samples_phy <- tax_glom(no_ra_samples, taxrank = "Phylum", NArm = F)
no_samples_phy
no_samples_phy_filt <- merge_low_abundance(no_samples_phy, threshold = 0.1)
no_samples_phy_filt #14 taxa
no_samples_phy_melt <- psmelt(no_samples_phy_filt)

#look at family
samples_family #399 taxa 145 samples
samples_family_filt <- merge_low_abundance(samples_family, threshold = 0.1)
samples_family_filt #54 taxa
samples_family_melt <- psmelt(samples_family_filt)


#calcualte abundance
datphy <- samples_phy_melt %>% group_by(Phylum) %>% 
  summarize(mean(Abundance, na.rm=TRUE),
            sd(Abundance, na.rm = T))
write.csv(datphy, "phylum abundance.csv")


#Plots
microbiome_dendro_plot <- ggplot(samples.dendro.data$segments) +
  theme_minimal() +
  labs(y= "") +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_point(data = samples.dendro.data$labels, aes(x=x,y=y, colour = body_site), size = 9, shape = 15, position = position_nudge()) +
  geom_point(data = samples.dendro.data$labels, aes(x=x, y=y, colour = Animal), size = 9, shape = 15, position = position_nudge(y = -0.5))+
  #geom_text(aes(x=x, y=y, label=label), data = samples.dendro.data$labels)+
  scale_colour_manual(values = combo_palette)+
  theme(plot.title = element_text(size = 20),
        #legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7))
microbiome_dendro_plot

#Dendro no LA
no_la_dendro_plot <-  ggplot(no_samples.dendro.data$segments) +
  theme_minimal() +
  labs(y= "") +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_point(data = no_samples.dendro.data$labels, aes(x=x,y=y, colour = sample_type), size = 9, shape = 15, position = position_nudge(y=-.02)) +
  #geom_text(aes(x=x, y=y, label=label), data = no_samples.dendro.data$labels)+
  scale_colour_manual(values = no_la_palette)+
  theme(plot.title = element_text(size = 20),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7))
no_la_dendro_plot

sample_type_dendro<- ggplot(samples.dendro.data$segments) +
  theme_minimal() +
  labs(y= "") +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_point(data = samples.dendro.data$labels, aes(x=x,y=y, colour = sample_type), size = 9, shape = 15, position = position_nudge(y=-.02)) +
  #geom_text(aes(x=x, y=y, label=label), data = samples.dendro.data$labels)+
  scale_x_discrete(limits = samples.dendro.order, expand = c(0.035,0,0,0))+
  scale_colour_manual(values = sample_type_palette)+
  theme(plot.title = element_text(size = 20),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7))
sample_type_dendro

ra_phy <- ggplot(samples_phy_melt, aes(x= Sample, y= Abundance, fill = Phylum)) +
  theme_minimal() + #facet_wrap(~Body_Site, scales = "free_x")+
  labs(y= "") +
  geom_bar(stat = "identity", colour = "black") +
  scale_x_discrete(limits = samples.dendro.order, expand = c(0.035,0,0,0)) +
  scale_fill_manual(values = randomColor(33)) +
  theme(#legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7))
ra_phy

#phylum no LA

no_ra_phy <- ggplot(no_samples_phy_melt, aes(x= Sample, y= Abundance, fill = Phylum)) +
  theme_minimal() + #facet_wrap(~Body_Site, scales = "free_x")+
  labs(y= "") +
  geom_bar(stat = "identity", colour = "black") +
  scale_x_discrete(limits = no_samples.dendro.order, expand = c(0.035,0,0,0)) +
  scale_fill_manual(values = randomColor(33)) +
  theme(#legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.y = element_line(size = 0.75, colour = "black"),
    axis.title.y = element_text(size = 24),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_line(colour = "black", size = 0.7))
no_ra_phy

family_pallete <- randomColor(66)

allsamples_ra <- ggplot(samples_family_melt, aes(x= Sample, y= Abundance, fill = Family)) +
  theme_minimal() + #facet_wrap(~Body_Site, scales = "free_x")+
  labs(y= "") +
  geom_bar(stat = "identity", colour = "black") +
  scale_x_discrete(limits = samples.dendro.order, expand = c(0.035,0,0,0)) +
  scale_fill_manual(values = family_pallete ) +
  theme(#legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7))
allsamples_ra

family_stats <- samples_family_melt %>% group_by(Family) %>% 
  summarize(mean(Abundance, na.rm=TRUE),
            sd(Abundance, na.rm = T))


# Compare Firmicutes vs Bacteroidetes -------------------------------------
#Firmicutes
firmicutes_sub <- subset_taxa(rel_abund, Phylum=="Firmicutes")

firmicutes_phyla <- tax_glom(firmicutes_sub, taxrank = "Phylum", NArm = F) %>%
  psmelt()
firmicutes_family <- tax_glom(firmicutes_sub, taxrank = "Family", NArm = F) %>%
  psmelt() #118

#make a custom palette
firmicutes_family1 <- tax_glom(firmicutes_sub, taxrank = "Family", NArm = F)
firm_pallete <- randomColor(118)
write.csv(firm_pallete, "RA_firm_palette.csv")
write.csv(tax_table(firmicutes_family1), "RA_firm_filt.csv")

ggplot(firmicutes_family, aes(x= sample_type, y= Abundance)) +
  theme_bw() + #coord_flip() +
  labs(y= "", x = "") +
  geom_bar(aes(fill= Family), stat = "summary", colour = "black") +
  geom_errorbar(firmicutes_phyla, mapping= aes(x= sample_type, y= Abundance), stat = "summary", width = 0.55, size = 0.6) +
  scale_fill_manual(values = randomColor(115)) +
  scale_x_discrete(limits = c("RUFL", "RUTI", "SIFL", "SITI", "LIFL", "LITI", "Fecal", "Liver_Abscess"),
                   labels = c("","","","","","","","","","","")) +
  scale_y_continuous(expand = c(0.002,0,0.1,0)) +
  ylim(0,100)+
  theme(#legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        panel.border = element_rect(size = 1.5, colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 28, vjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 28, colour = "black"),
        axis.ticks = element_line(size = 0.9, colour = "black"))

firm_anova <- pairwise.wilcox.test(firmicutes_phyla$Abundance, firmicutes_phyla$sample_type, p.adjust.method = "BH")

#superscript
firm_pvalues <- firm_anova$p.value %>% 
  fullPTable()
multcompLetters(firm_pvalues)

datfirm <- firmicutes_family %>% group_by(Family) %>% 
  summarize(mean(Abundance, na.rm=TRUE),
            sd(Abundance, na.rm = T))
write.csv(datfirm,"RA Firm.csv")
#Bacteroidetes
bacteroidetes_sub <- subset_taxa(rel_abund, Phylum=="Bacteroidota")

bacteroidetes_phyla <- tax_glom(bacteroidetes_sub, taxrank = "Phylum", NArm = F) %>%
  psmelt()
bacteroidetes_family <- tax_glom(bacteroidetes_sub, taxrank = "Family", NArm = F) %>%
  psmelt()

#build palette
bacteroidetes_family1 <- tax_glom(bacteroidetes_sub, taxrank = "Family", NArm = F)
bac_pallete <- randomColor(57)
write.csv(bac_pallete, "RA_bac_palette.csv")
write.csv(tax_table(bacteroidetes_family1), "RA_bac_filt.csv")

ggplot(bacteroidetes_family, aes(x= sample_type, y= Abundance)) +
  theme_bw() + #coord_flip() +
  labs(y="", x= "") +
  geom_bar(aes(fill= Family), stat = "summary", colour = "black") +
  geom_errorbar(bacteroidetes_phyla, mapping= aes(x= sample_type, y= Abundance), stat = "summary", width = 0.55, size = 0.6) +
  scale_fill_manual(values = randomColor(115)) +
  scale_x_discrete(limits = c("RUFL", "RUTI", "SIFL", "SITI", "LIFL", "LITI", "Fecal", "Liver_Abscess"),
                   labels = c("","","","","","","","","","","")) +
  scale_y_continuous(expand = c(0.002,0,0.1,0)) +
  ylim(0,100)+
  theme(#legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        panel.border = element_rect(size = 1.5, colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 28, vjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 28, colour = "black"),
        axis.ticks = element_line(size = 0.9, colour = "black"))

bac_anova <- pairwise.wilcox.test(bacteroidetes_phyla$Abundance, bacteroidetes_phyla$sample_type, p.adjust.method = "BH")
#superscript
bac_pvalues <- bac_anova$p.value %>% 
  fullPTable()
multcompLetters(bac_pvalues)

datbacteroidetes <- bacteroidetes_family %>% group_by(Family) %>% 
  summarize(mean(Abundance, na.rm=TRUE),
            sd(Abundance, na.rm = T))
write.csv(datbacteroidetes, "RA Bacter.csv")

#Actinobacteriota
actino_sub <- subset_taxa(rel_abund, Phylum=="Actinobacteriota")

actino_phy <- tax_glom(actino_sub, taxrank = "Phylum", NArm = F) %>% 
  psmelt()
actino_phy
actino_family <- tax_glom(actino_sub, taxrank = "Family", NArm = F) %>% 
  psmelt()

#make palette
actino_pallete <- randomColor(73)
actino_family1 <- tax_glom(actino_sub, taxrank = "Family", NArm = F)
write.csv(actino_pallete, "RA_actino_palette.csv")
write.csv(tax_table(actino_family1), "RA_actino_filt.csv")


ggplot(actino_family, aes(x= sample_type, y= Abundance)) +
  theme_bw() + #coord_flip() +
  labs(y="", x= "") +
  geom_bar(aes(fill= Family), stat = "summary", colour = "black") +
  geom_errorbar(actino_phy, mapping= aes(x= sample_type, y= Abundance), stat = "summary", width = 0.55, size = 0.6) +
  scale_fill_manual(values = randomColor(115)) +
  scale_x_discrete(limits = c("RUFL", "RUTI", "SIFL", "SITI", "LIFL", "LITI", "Fecal", "Liver_Abscess"),
                   labels = c("", "", "", "", "", "", "", "", "","", "")) +
  scale_y_continuous(expand = c(0.002,0,0.1,0)) +
  ylim(0,100)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
    panel.border = element_rect(size = 1.5, colour = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text = element_text(size = 14, colour = "black"),
    axis.title.y = element_text(size = 28, vjust = 0.5),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 28, colour = "black"),
    axis.ticks = element_line(size = 0.9, colour = "black"))

act_anova <- pairwise.wilcox.test(actino_phy$Abundance, actino_phy$sample_type, p.adjust.method = "BH")
act_anova
#superscript
act_pvalues <- act_anova$p.value %>% 
  fullPTable()
multcompLetters(act_pvalues)

datactino <- actino_family %>% group_by(Family) %>% 
  summarize(mean(Abundance, na.rm=TRUE),
            sd(Abundance, na.rm = T))
write.csv(datactino, "RA Actino.csv")



# search for important families -------------------------------------------

#fuso
fuso <- subset_taxa(ra_family, Family=="Fusobacteriaceae") %>% 
  psmelt()

ggplot(fuso, aes(x= sample_type, y= Abundance, fill = sample_type)) +
  theme_bw() + labs(y= "Relative Abundance (%)", x="", title = "") +
  scale_x_discrete(limits=c("RUFL", "RUTI", "SIFL", "SITI", "LIFL", "LITI", "Fecal", "Liver_Abscess"))+
  geom_boxplot() +
  geom_point()+
  scale_fill_manual(values = sample_type_palette )+
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

kruskal.test(fuso$Abundance, fuso$sample_type)
fuso_anova <- pairwise.wilcox.test(fuso$Abundance,fuso$sample_type, p.adjust.method = "BH")
#superscript
fuso_pvalues <- fuso_anova$p.value %>% 
  fullPTable()
multcompLetters(fuso_pvalues)

datfuso <- fuso %>% group_by(sample_type) %>% 
  summarize(mean(Abundance, na.rm=T),
            sd(Abundance, na.rm = T))


#actino
actino <- subset_taxa(ra_family, Family=="Actinomycetaceae") %>% 
  psmelt()

ggplot(actino, aes(x= sample_type, y= Abundance, fill = sample_type)) +
  theme_bw() + labs(y= "Relative Abundance (%)", x="", title = "") +
  scale_x_discrete(limits=c("RUFL", "RUTI", "SIFL", "SITI", "LIFL", "LITI", "Fecal", "Liver_Abscess"))+
  geom_boxplot() +
  geom_point()+
  scale_fill_manual(values = sample_type_palette )+
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

kruskal.test(actino$Abundance, actino$sample_type)
actino_anova <- pairwise.wilcox.test(actino$Abundance,actino$sample_type, p.adjust.method = "BH")
#superscript
actino_pvalues <- actino_anova$p.value %>% 
  fullPTable()
multcompLetters(actino_pvalues)

datactino <- actino %>% group_by(sample_type) %>% 
  summarize(mean(Abundance, na.rm=T),
            sd(Abundance, na.rm = T))

#Bacteroiodaceae

bac <- subset_taxa(ra_family, Family=="Bacteroidaceae") %>% 
  psmelt()

ggplot(bac, aes(x= sample_type, y= Abundance, fill = sample_type)) +
  theme_bw() + labs(y= "Relative Abundance (%)", x="", title = "") +
  scale_x_discrete(limits=c("RUFL", "RUTI", "SIFL", "SITI", "LIFL", "LITI", "Fecal", "Liver_Abscess"))+
  geom_boxplot() +
  geom_point()+
  scale_fill_manual(values = sample_type_palette )+
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

kruskal.test(bac$Abundance, bac$sample_type)
bac_anova <- pairwise.wilcox.test(bac$Abundance,bac$sample_type, p.adjust.method = "BH")
#superscript
bac_pvalues <- bac_anova$p.value %>% 
  fullPTable()
multcompLetters(bac_pvalues)

datbac <- bac %>% group_by(sample_type) %>% 
  summarize(mean(Abundance, na.rm=T),
            sd(Abundance, na.rm = T))


# Core Analysis -----------------------------------------------------------

#transform data to check for core

core_ra <- transform(samples_family,"compositional")
core_ra_family <- aggregate_taxa(core_ra,"Family")
core_ra_family_filt <- core(core_ra_family, detection= 0.001, prevalence =.9)
core_ra_family_filt

#plot
# core with compositionals
prevalences <- seq(0.1,1,0.1)
detections <- 10^seq(log10(1e-3), log10(0.1), length = 6)

detections
detections <- trunc(detections*10^3)/10^3

p1 <- plot_core(core_ra_family,
                plot.type = "heatmap",
                colours = rev(brewer.pal(10, "Spectral")),
                prevalences = prevalences,
                detections = detections, min.prevalence = .90)
p1 + theme_bw() +
  labs(y= "", x= "Detection Threshold \n (Relative Abundance (%))", title = "") +
  theme(legend.key.size = unit(4, "lines"),
        plot.title = element_text(size = 24),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.border = element_rect(colour = "black", size = 1),
        axis.ticks = element_line(size = 0.75, colour = "black"),
        axis.text.y = element_text(colour = "black", size = 26),
        axis.text.x = element_text(colour = "black", size = 24),
        axis.title.x = element_text(size = 40, vjust = -0.75))


# Core new way ------------------------------------------------------------

core_no <- subset_samples(core_ra_family, body_site %in% c("Rumen", "Jejunum", "Colon", "Fecal"))
core_no_filt <- core(core_no, detection = 0.001, prevalence = 0.9)
core_no_filt@tax_table

core_taxa_no_la <- core_no_filt %>%
  psmelt() %>%
  group_by(OTU, sample_type) %>%
  summarize(
    avg = mean(Abundance, na.rm = TRUE),
    sd = sd(Abundance, na.rm = TRUE),
    lower_ci = avg - qt(0.975, df = n() - 1) * sd / sqrt(n()),
    upper_ci = avg + qt(0.975, df = n() - 1) * sd / sqrt(n())
  )
View(core_taxa_no_la)

all_gi <- ggplot(core_taxa_no_la, aes(x= OTU, y = avg, fill = sample_type))+
  geom_point(aes(),position = position_dodge(width = .7), size = 7, shape = 21, alpha = 1.5)+
  scale_fill_manual(values = no_la_palette )+
  geom_linerange(mapping =aes(ymin=lower_ci,ymax=upper_ci),position = position_dodge(width = .7), linetype = 1)+
  labs(x = "", y = "" )+
  scale_y_continuous(breaks = seq(0, 0.35, by = 0.05))+
  theme_bw()+
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", size = 1),
        title = element_text(size = 28),
        axis.ticks = element_line(colour = "black", size = 0.75),
        axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(size = 24),
        axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
        axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1))
all_gi


# Core taxa by region -----------------------------------------------------

#rumen
rumen_core <- subset_samples(core_ra_family, body_site=="Rumen")
rumen_ra <- transform(rumen_core, "compositional")
core_rumen_filt <- core(rumen_ra, detection= 0.001, prevalence =.9)
core_rumen_filt


p3 <- plot_core(rumen_ra,
                plot.type = "heatmap",
                colours = rev(brewer.pal(10, "Spectral")),
                prevalences = prevalences,
                detections = detections, min.prevalence = .90)
p3 + theme_bw() +
  labs(y= "", x= "Detection Threshold \n (Relative Abundance (%))", title = "") +
  theme(legend.key.size = unit(4, "lines"),
        legend.position = "none",
        plot.title = element_text(size = 24),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.border = element_rect(colour = "black", size = 1),
        axis.ticks = element_line(size = 0.75, colour = "black"),
        axis.text.y = element_text(colour = "black", size = 26),
        axis.text.x = element_text(colour = "black", size = 24),
        axis.title.x = element_text(size = 40, vjust = -0.75))
#no labels
p3 + theme_bw() +
  labs(y = "", x = "", title = "") +
  theme(
    legend.key.size = unit(4, "lines"),
    legend.position = "none",
    plot.title = element_text(size = 24),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    panel.border = element_rect(colour = "black", size = 1),
    axis.ticks = element_line(size = 0.75, colour = "black"),
    axis.text.y = element_text(colour = "black", size = 26),
    axis.text.x = element_text(colour = "black", size = 24),
    axis.title.x = element_text(size = 40, vjust = -0.75)
  )

#bubble plot
rumen_core_new <- core_rumen_filt %>%
  psmelt() %>%
  group_by(OTU, sample_type)

rumen_core_plot<- ggplot(rumen_core_new, aes(x= OTU, y = Abundance, fill = sample_type))+
  geom_point(aes(),position = position_jitterdodge(), size = 5, shape = 21, alpha = .9)+
  scale_fill_manual(values = rumen_core_palette )+
  #geom_linerange(mapping =aes(ymin=lower_ci,ymax=upper_ci),position = position_dodge(width = .5), linetype = 1)+
  labs(x = "", y = "" )+
  scale_y_continuous(breaks = seq(0, 1, by = 0.05))+
  theme_bw()+
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", size = 1),
        title = element_text(size = 28),
        axis.ticks = element_line(colour = "black", size = 0.75),
        axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(size = 24),
        axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1))

rumen_core_plot



#SI
si_core <- subset_samples(core_ra_family, body_site=="Jejunum")
si_ra <- transform(si_core, "compositional")
core_si <- core(si_ra, detection = 0.001, prevalence = 0.90)
core_si

p6 <- plot_core(si_ra,
                plot.type = "heatmap",
                colours = rev(brewer.pal(10, "Spectral")),
                prevalences = prevalences,
                detections = detections, min.prevalence = .90)
p6 + theme_bw() +
  labs(y= "", x= "Detection Threshold \n (Relative Abundance (%))", title = "") +
  theme(legend.key.size = unit(4, "lines"),
        plot.title = element_text(size = 24),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.border = element_rect(colour = "black", size = 1),
        axis.ticks = element_line(size = 0.75, colour = "black"),
        axis.text.y = element_text(colour = "black", size = 26),
        axis.text.x = element_text(colour = "black", size = 24),
        axis.title.x = element_text(size = 40, vjust = -0.75))
#no lables
p6 + theme_bw() +
  labs(y = "", x = "", title = "") +
  theme(
    legend.key.size = unit(4, "lines"),
    legend.position = "none",
    plot.title = element_text(size = 24),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    panel.border = element_rect(colour = "black", size = 1),
    axis.ticks = element_line(size = 0.75, colour = "black"),
    axis.text.y = element_text(colour = "black", size = 26),
    axis.text.x = element_text(colour = "black", size = 24),
    axis.title.x = element_text(size = 40, vjust = -0.75)
  )


#si bubble
#bubble plot
si_core_new <- core_si %>%
  psmelt() %>%
  group_by(OTU, sample_type)

si_core_plot<- ggplot(si_core_new, aes(x= OTU, y = Abundance, fill = sample_type))+
  geom_point(aes(),position = position_jitterdodge(), size = 5, shape = 21, alpha = .9)+
  scale_fill_manual(values = si_core_palette )+
  #geom_linerange(mapping =aes(ymin=lower_ci,ymax=upper_ci),position = position_dodge(width = .5), linetype = 1)+
  labs(x = "", y = "" )+
  scale_y_continuous(breaks = seq(0, 1, by = 0.05))+
  theme_bw()+
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", size = 1),
        title = element_text(size = 28),
        axis.ticks = element_line(colour = "black", size = 0.75),
        axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(size = 24),
        axis.text.x = element_blank())

si_core_plot


hindgut_core
#hindgut
hindgut_core <- subset_samples(core_ra_family, body_site %in% c("Colon", "Fecal"))
hindgut_ra <- transform(hindgut_core, "compositional")
core_hindgut <- core(hindgut_ra, detection = 0.001, prevalence = 0.90)
core_hindgut

p5 <- plot_core(hindgut_ra,
                plot.type = "heatmap",
                colours = rev(brewer.pal(10, "Spectral")),
                prevalences = prevalences,
                detections = detections, min.prevalence = .90)
p5 + theme_bw() +
  labs(y= "", x= "Detection Threshold \n (Relative Abundance (%))", title = "") +
  theme(legend.key.size = unit(4, "lines"),
        plot.title = element_text(size = 24),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.border = element_rect(colour = "black", size = 1),
        axis.ticks = element_line(size = 0.75, colour = "black"),
        axis.text.y = element_text(colour = "black", size = 26),
        axis.text.x = element_text(colour = "black", size = 24),
        axis.title.x = element_text(size = 40, vjust = -0.75))

#colon
colon_core <- subset_samples(core_ra_family, body_site=="Colon")
colon_ra <- transform(colon_core, "compositional")
core_colon <- core(colon_ra, detection = 0.001, prevalence = 0.90)
core_colon

p7 <- plot_core(colon_ra,
                plot.type = "heatmap",
                colours = rev(brewer.pal(10, "Spectral")),
                prevalences = prevalences,
                detections = detections, min.prevalence = .90)
p7 + theme_bw() +
  labs(y= "", x= "Detection Threshold \n (Relative Abundance (%))", title = "") +
  theme(legend.key.size = unit(4, "lines"),
        plot.title = element_text(size = 24),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.border = element_rect(colour = "black", size = 1),
        axis.ticks = element_line(size = 0.75, colour = "black"),
        axis.text.y = element_text(colour = "black", size = 26),
        axis.text.x = element_text(colour = "black", size = 24),
        axis.title.x = element_text(size = 40, vjust = -0.75))
#no lables
p7 + theme_bw() +
  labs(y = "", x = "", title = "") +
  theme(
    legend.key.size = unit(4, "lines"),
    legend.position = "none",
    plot.title = element_text(size = 24),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    panel.border = element_rect(colour = "black", size = 1),
    axis.ticks = element_line(size = 0.75, colour = "black"),
    axis.text.y = element_text(colour = "black", size = 26),
    axis.text.x = element_text(colour = "black", size = 24),
    axis.title.x = element_text(size = 40, vjust = -0.75)
  )

#new colon core
#bubble plot
colon_core_new <- core_colon %>%
  psmelt() %>%
  group_by(OTU, sample_type)

colon_core_plot<- ggplot(colon_core_new, aes(x= OTU, y = Abundance, fill = sample_type))+
  geom_point(aes(),position = position_jitterdodge(), size = 5, shape = 21, alpha = .9)+
  scale_fill_manual(values = li_core_palette )+
  #geom_linerange(mapping =aes(ymin=lower_ci,ymax=upper_ci),position = position_dodge(width = .5), linetype = 1)+
  labs(x = "", y = "" )+
  scale_y_continuous(breaks = seq(0, 1, by = 0.05))+
  theme_bw()+
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", size = 1),
        title = element_text(size = 28),
        axis.ticks = element_line(colour = "black", size = 0.75),
        axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(size = 24),
        axis.text.x = element_blank())

colon_core_plot



# Split lumen and mucosal -------------------------------------------------
##rumen tissue
rumen_tis_core <- subset_samples(core_ra_family, sample_type=="RUTI")
core_ra_rumen_ti <- transform(rumen_tis_core, "compositional")
core_rumen_tis_filt <- core(core_ra_rumen_ti, detection= 0.001, prevalence =.9)
core_rumen_tis_filt

p2 <- plot_core(core_ra_rumen_ti,
                plot.type = "heatmap",
                colours = rev(brewer.pal(10, "Spectral")),
                prevalences = prevalences,
                detections = detections, min.prevalence = .90)
p2 + theme_bw() +
  labs(y= "", x= "Detection Threshold \n (Relative Abundance (%))", title = "") +
  theme(legend.key.size = unit(4, "lines"),
        plot.title = element_text(size = 24),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.border = element_rect(colour = "black", size = 1),
        axis.ticks = element_line(size = 0.75, colour = "black"),
        axis.text.y = element_text(colour = "black", size = 26),
        axis.text.x = element_text(colour = "black", size = 24),
        axis.title.x = element_text(size = 40, vjust = -0.75))

##rumen fluid
rufl_core <- subset_samples(core_ra_family, sample_type=="RUFL")
core_ra_rumen_fl <- transform(rufl_core, "compositional")
core_rumen_fl_filt <- core(core_ra_rumen_fl, detection = 0.001, prevalence = .9)
core_rumen_fl_filt

p9 <- plot_core(core_ra_rumen_fl,
                plot.type = "heatmap",
                colours = rev(brewer.pal(10, "Spectral")),
                prevalences = prevalences,
                detections = detections, min.prevalence = .90)
p9 + theme_bw() +
  labs(y= "", x= "Detection Threshold \n (Relative Abundance (%))", title = "") +
  theme(legend.key.size = unit(4, "lines"),
        plot.title = element_text(size = 24),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.border = element_rect(colour = "black", size = 1),
        axis.ticks = element_line(size = 0.75, colour = "black"),
        axis.text.y = element_text(colour = "black", size = 26),
        axis.text.x = element_text(colour = "black", size = 24),
        axis.title.x = element_text(size = 40, vjust = -0.75))

##SITI
siti_core <- subset_samples(core_ra_family, sample_type=="SITI")
core_ra_siti <- transform(siti_core, "compositional")
core_si_filt <- core(core_ra_siti, detection = 0.001, prevalence = 0.90)
core_si_filt

p10 <- plot_core(core_ra_siti,
                plot.type = "heatmap",
                colours = rev(brewer.pal(10, "Spectral")),
                prevalences = prevalences,
                detections = detections, min.prevalence = .90)
p10 + theme_bw() +
  labs(y= "", x= "Detection Threshold \n (Relative Abundance (%))", title = "") +
  theme(legend.key.size = unit(4, "lines"),
        plot.title = element_text(size = 24),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.border = element_rect(colour = "black", size = 1),
        axis.ticks = element_line(size = 0.75, colour = "black"),
        axis.text.y = element_text(colour = "black", size = 26),
        axis.text.x = element_text(colour = "black", size = 24),
        axis.title.x = element_text(size = 40, vjust = -0.75))
##SIFL
sitfl_core <- subset_samples(core_ra_family, sample_type=="SIFL")
core_ra_sifl <- transform(sitfl_core, "compositional")
core_sifl_filt <- core(core_ra_sifl, detection = 0.001, prevalence = 0.90)
core_sifl_filt

p11 <- plot_core(core_ra_sifl,
                plot.type = "heatmap",
                colours = rev(brewer.pal(10, "Spectral")),
                prevalences = prevalences,
                detections = detections, min.prevalence = .90)
p11 + theme_bw() +
  labs(y= "", x= "Detection Threshold \n (Relative Abundance (%))", title = "") +
  theme(legend.key.size = unit(4, "lines"),
        plot.title = element_text(size = 24),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.border = element_rect(colour = "black", size = 1),
        axis.ticks = element_line(size = 0.75, colour = "black"),
        axis.text.y = element_text(colour = "black", size = 26),
        axis.text.x = element_text(colour = "black", size = 24),
        axis.title.x = element_text(size = 40, vjust = -0.75))
#colon ti
liti_core <- subset_samples(core_ra_family, sample_type=="LITI")
core_ra_liti <- transform(liti_core, "compositional")
core_liti_fil <- core(core_ra_liti, detection = 0.001, prevalence = 0.90)
core_liti_fil

p12 <- plot_core(core_ra_liti,
                plot.type = "heatmap",
                colours = rev(brewer.pal(10, "Spectral")),
                prevalences = prevalences,
                detections = detections, min.prevalence = .90)
p12 + theme_bw() +
  labs(y= "", x= "Detection Threshold \n (Relative Abundance (%))", title = "") +
  theme(legend.key.size = unit(4, "lines"),
        plot.title = element_text(size = 24),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.border = element_rect(colour = "black", size = 1),
        axis.ticks = element_line(size = 0.75, colour = "black"),
        axis.text.y = element_text(colour = "black", size = 26),
        axis.text.x = element_text(colour = "black", size = 24),
        axis.title.x = element_text(size = 40, vjust = -0.75))
#colon fl
lifl_core <- subset_samples(core_ra_family, sample_type== "LIFL")
core_ra_lifl <- transform(lifl_core, "compositional")
core_lifl_fil <- core(core_ra_lifl, detection = 0.001, prevalence = 0.90)
core_lifl_fil

p13 <- plot_core(core_ra_lifl,
                 plot.type = "heatmap",
                 colours = rev(brewer.pal(10, "Spectral")),
                 prevalences = prevalences,
                 detections = detections, min.prevalence = .90)
p13 + theme_bw() +
  labs(y= "", x= "Detection Threshold \n (Relative Abundance (%))", title = "") +
  theme(legend.key.size = unit(4, "lines"),
        plot.title = element_text(size = 24),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.border = element_rect(colour = "black", size = 1),
        axis.ticks = element_line(size = 0.75, colour = "black"),
        axis.text.y = element_text(colour = "black", size = 26),
        axis.text.x = element_text(colour = "black", size = 24),
        axis.title.x = element_text(size = 40, vjust = -0.75))
#fecal
fecal_core<- subset_samples(core_ra_family, body_site== "Fecal")
fecal_ra <- transform(fecal_core, "compositional")
core_fecal <- core(fecal_ra, detection = 0.001, prevalence = 0.90)
core_fecal

p8 <- plot_core(fecal_ra,
                plot.type = "heatmap",
                colours = rev(brewer.pal(10, "Spectral")),
                prevalences = prevalences,
                detections = detections, min.prevalence = .90)
p8 + theme_bw() +
  labs(y= "", x= "Detection Threshold \n (Relative Abundance (%))", title = "") +
  theme(legend.key.size = unit(4, "lines"),
        plot.title = element_text(size = 24),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.border = element_rect(colour = "black", size = 1),
        axis.ticks = element_line(size = 0.75, colour = "black"),
        axis.text.y = element_text(colour = "black", size = 26),
        axis.text.x = element_text(colour = "black", size = 24),
        axis.title.x = element_text(size = 40, vjust = -0.75))

#no lables
p8 + theme_bw() +
  labs(y = "", x = "", title = "") +
  theme(
    legend.key.size = unit(4, "lines"),
    legend.position = "none",
    plot.title = element_text(size = 24),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    panel.border = element_rect(colour = "black", size = 1),
    axis.ticks = element_line(size = 0.75, colour = "black"),
    axis.text.y = element_text(colour = "black", size = 26),
    axis.text.x = element_text(colour = "black", size = 24),
    axis.title.x = element_text(size = 40, vjust = -0.75)
  )

#fecal new plot
#bubble plot
fecal_core_new <- core_fecal %>%
  psmelt() %>%
  group_by(OTU, sample_type)

fecal_core_plot<- ggplot(fecal_core_new, aes(x= OTU, y = Abundance, fill = sample_type))+
  geom_point(aes(),position = position_jitterdodge(), size = 5, shape = 21, alpha = .9)+
  scale_fill_manual(values = fecal_core_palette )+
  #geom_linerange(mapping =aes(ymin=lower_ci,ymax=upper_ci),position = position_dodge(width = .5), linetype = 1)+
  labs(x = "", y = "" )+
  scale_y_continuous(breaks = seq(0, 1, by = 0.05))+
  theme_bw()+
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", size = 1),
        title = element_text(size = 28),
        axis.ticks = element_line(colour = "black", size = 0.75),
        axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(size = 24),
        axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1))
fecal_core_plot





#liver abscess
la_core <- subset_samples(core_ra_family, body_site== "Liver_Abscess")
La_ra <- transform(la_core, "compositional")
core_la <- core(la_core, detection = 0.001, prevalence = 0.90)
core_la@otu_table

