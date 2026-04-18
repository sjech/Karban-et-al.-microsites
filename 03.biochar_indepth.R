# The purpose of this code is to think about whether or not biochar affects the soil microbiome two years after pit installation. Here we are only comparing pits without biochar to pits with biochar.
# On 8/23/23, SDJ added code to think about some other ideas for the manuscript

### Load libraries
library(tidyverse)
library(phyloseq) # For all bioinformatics analyses post Dada2
packageVersion('phyloseq')
library(ggplot2)
library(vegan)
library(car) # For Levene Test and anova
#library(PMCMRplus) # For Nemenyi posthoc test
library(FSA)
#library("DESeq2")
library(indicspecies) # For indicator species
library(microbiome)
library(ggpubr) #for ggarrange
library(stringr)
library(RVAideMemoire)

# Functions
se <- function(x) {sd(x,na.rm=TRUE)/sqrt(length(x))}

# Color Palettes
grafify_muted_greenstart <- c("#117733", "#88ccee", "#882255", "#44aa99", "#999933", "#aa4499", "#dddddd","#cc6677", "#332288", "#ddcc77")


### Load the phyloseq data table which comes from the pits16Sphyloseq_1.10.22.R script (note the date should be 1.10.23 - I got the year wrong). 
list.files()
qPCR_masses <- read.csv("qPCR_mass_extracted.csv", header = TRUE)
# these two were when I was rarefying and pruning
# pits_16S_rar_filt_pruned_ps <- readRDS("pits_processed_rare_filt_pruned.rds")
# pits_16S_rar_filt_pruned_ps@sam_data
# these are more current
pits_16S_filt_ps <- readRDS("pits_processed_filt_ps.rds")
pits_16S_filt_ps@sam_data

# Let's compare the microbial communities from 2022 in controls and pits
samplesToKeep <- c("BP22a","BP22b","BP22c","P22a","P22b","P22c","C22a","C22b","C22c")
pitBiochar <- prune_samples(sample_names(pits_16S_filt_ps) %in% samplesToKeep , pits_16S_filt_ps)

pitBiochar@sam_data$bulk_soil <- factor(pitBiochar@sam_data$bulk_soil, levels = c("C","P","BP"))

# Is 16S rRNA gene abundance different for the three treatments?
sample_variables(pitBiochar)
qPCR <- as.data.frame(pitBiochar@sam_data$qPCR_16Scopies)
qPCR$treatment <- c("BP","BP","BP","C","C","C","P","P","P")
qPCR$Sample_ID <- sample_names(pitBiochar)
colnames(qPCR)[1] <- c("qPCR_count")
qPCR <- left_join(qPCR, qPCR_masses, by = c("Sample_ID" = "SampleID"))

# Statistical test 
qPCR.aov <- aov(qPCR_count ~ treatment, data = qPCR)
shapiro.test(qPCR.aov$residuals) # ok
leveneTest(qPCR_count ~ treatment, data = qPCR)
#plot(qPCR.aov)
summary(qPCR.aov) # not significantly different (F=0.602, p = 0.578)
mean(qPCR$qPCR_count)
se(qPCR$qPCR_count)


# Richness 
# You have to rarefy before calculating richness
# mbcyano_16S_ps # this is the phyloseq object to start with
min(sample_sums(pitBiochar)) # 19473 is the minimum
max(sample_sums(pitBiochar)) # 56395 is the maxiumum
samplesumsdf <- as.matrix(sample_sums(pitBiochar))
samplesumsdf <- as.data.frame(samplesumsdf)
samplesumsdf$samples <- rownames(samplesumsdf)
samplesumsdf <- arrange(samplesumsdf, -V1)

samplesumsdf %>% ggplot()+
  geom_point(mapping = aes(x = V1, y = reorder(samples,V1)))+
  geom_vline(mapping = aes(xintercept = 19000))+
  theme_classic()
# by rarefying to 19,000 we do not lose any samples

# Check read counts 
colSums(otu_table(pitBiochar)) # all read counts are different
sum(colSums(otu_table(pitBiochar))) # total number of reads in the dataset is 277,501
# remove ones less than 1% abundance in any sample
0.0001*277501 # 277.501 reads

# if using an abundance cutoff, you should apply it here too. 
# old: to keep consistent with the abundance cutoff I used when assigning GS names, I'd like to impose an abundance filter and include only the taxa that have an Abundance > 0.1% (raw reads = 20) in any sample 

pitBiochar_filt = filter_taxa(pitBiochar, function(x) max(x) >= 20, TRUE)
taxa_sums(pitBiochar_filt)
taxa_sums(pitBiochar)
min(taxa_sums(pitBiochar))
min(taxa_sums(pitBiochar_filt))
max(taxa_sums(pitBiochar))
max(taxa_sums(pitBiochar_filt))
mean(taxa_sums(pitBiochar))
mean(taxa_sums(pitBiochar_filt))

min(sample_sums(pitBiochar_filt)) # 19057 is the minimum
max(sample_sums(pitBiochar_filt)) # 54698 is the maxiumum
0.1*54698  #if rarefied to 10% of the max
0.9*19057 # if rarefied to 90% of the minimum, the level would be 17151
samplesumsdf <- as.matrix(sample_sums(pitBiochar_filt))
samplesumsdf <- as.data.frame(samplesumsdf)
samplesumsdf$samples <- rownames(samplesumsdf)
samplesumsdf <- arrange(samplesumsdf, -V1)
samplesumsdf %>% ggplot()+
  geom_point(mapping = aes(x = V1, y = reorder(samples,V1)))+
  geom_vline(mapping = aes(xintercept = 17151))+
  theme_classic()

# Rarefying required to assess richness
set.seed(1223)
pitBiochar_filt.rarefied <- rarefy_even_depth(pitBiochar_filt, sample.size = 17151, replace = FALSE)
sample_sums(pitBiochar_filt.rarefied) # they are all at 17151

# Now run all of the stuff needed to calculate richness
# melt out of the phyloseq object 
rare_psmelt <- psmelt(pitBiochar_filt.rarefied)

# remove rows with Abundance == 0
rare_psmelt_nozero <- rare_psmelt %>% filter(Abundance > 0) # this is a presence/absence filter which you have to do in order to use the summarise function to calculate richness

# ASV counts
ASVcount <- rare_psmelt_nozero %>%
  group_by(Sample, SampleID_norep, bulk_soil) %>%
  summarise(ASV_count = n())

# Richness Summaries 
DiversityIndices <- ASVcount %>%
  group_by(bulk_soil) %>%
  summarise(meanRichness_ASV = mean(ASV_count),
            sdRichness_ASV = sd(ASV_count),
            seRichness_ASV = se(ASV_count))

# ASV-level counts grouped by Phylum
ASVcount_byPhyla <- rare_psmelt_nozero %>%
  group_by(Sample, SampleID_norep, bulk_soil, Phylum) %>%
  summarise(ASV_count = n())

# check ricness at phyla level
# Bacter
richness_bacter <- ASVcount_byPhyla %>% filter(Phylum == "Bacteroidota")
kruskal.test(ASV_count ~ bulk_soil, data = richness_bacter) #significant
# Kruskal-Wallis chi-squared = 7.2, df = 2, p-value = 0.02732
dunnTest(ASV_count ~ bulk_soil, data = richness_bacter)
#  BP - C  2.683282 0.007290358 0.02187107
# Actino
richness_actino <- ASVcount_byPhyla %>% filter(Phylum == "Actinobacteriota")
kruskal.test(ASV_count ~ bulk_soil, data = richness_actino) #not significant
# Kruskal-Wallis chi-squared = 5.4678, df = 2, p-value = 0.06497
#Chloro
richness_chloro <- ASVcount_byPhyla %>% filter(Phylum == "Chloroflexi")
kruskal.test(ASV_count ~ bulk_soil, data = richness_chloro) #significant
# Kruskal-Wallis chi-squared = 7.2, df = 2, p-value = 0.02732
dunnTest(ASV_count ~ bulk_soil, data = richness_chloro)
#  C - P  2.683282 0.007290358 0.02187107
#Proteo
richness_proteo <- ASVcount_byPhyla %>% filter(Phylum == "Proteobacteria")
kruskal.test(ASV_count ~ bulk_soil, data = richness_proteo) #not significant
# Kruskal-Wallis chi-squared = 5.6, df = 2, p-value = 0.06081


DiversityIndices_byPhyla <- ASVcount_byPhyla %>%
  group_by(bulk_soil, Phylum) %>%
  summarise(meanRichness_ASV = mean(ASV_count),
            sdRichness_ASV = sd(ASV_count),
            seRichness_ASV = se(ASV_count))

DiversityIndices_byPhyla %>%
  ggplot()+
  geom_col(mapping = aes(x = Phylum, y = meanRichness_ASV, fill = bulk_soil), position = position_dodge())+
  geom_errorbar(mapping = aes(x = Phylum, ymin = meanRichness_ASV-seRichness_ASV, ymax = meanRichness_ASV + seRichness_ASV, color = bulk_soil), position = position_dodge())+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))

# Statistical Assessment for Richness
# ANOVA three levels: control vs. pits vs. pits+biochar
leveneTest(ASVcount$ASV_count ~ ASVcount$bulk_soil) # Variance homogeneous (p > 0.05)
# greater than 0.05 is good - data is homogeneous
shapiro.test(ASVcount$ASV_count) # Richness normally distributed (p > 0.05)
# p >0.05 is good - richness is normally distributed
m <- aov(ASVcount$ASV_count ~ ASVcount$bulk_soil)
shapiro.test(m$residuals) # all good 
plot(m) # QQ plot is no good - variance looks different
summary(m) #significant
# post-hoc test if there is a significant difference
TukeyHSD(m)
# BP-C 127.66667   33.04449 222.28884 0.0143437
# BP-P 178.33333   83.71116 272.95551 0.0028279

# rank based test to compare richness for the bulk soils - do this if the data cannot be analyzed with ANOVA
kruskal.test(ASV_count ~ bulk_soil, data = ASVcount) # significant
# Kruskal-Wallis chi-squared = 6.5434, df = 2, p-value = 0.03794
dunnTest(ASV_count ~ bulk_soil, data = ASVcount, method = "holm")
# BP - P 2.544836 0.01093291 0.03279872

# Graph the richness 
DiversityIndices %>%
  ggplot()+
  geom_col(mapping = aes(x = bulk_soil, y = meanRichness_ASV, fill = bulk_soil))+
  geom_errorbar(mapping = aes(x = bulk_soil, ymin = meanRichness_ASV-seRichness_ASV, ymax = meanRichness_ASV + seRichness_ASV))+ # width = 0.5
  theme_classic() + 
  labs(x = "Bulk Soil Type", y = "Richness (ASVlevel)") + 
  scale_fill_manual(values = c("orange4","darkslategray", "goldenrod")) + 
  scale_x_discrete(labels = c("Control","Pit", "Pit+Biochar"))+
  theme(axis.text.x = element_text(size=12), 
        axis.text.y = element_text(size=12), 
        axis.title.y = element_text(size=14), 
        axis.title.x = element_text(size=14),
        legend.position = "none")

#ggsave("output/richnessbars.pdf", width = 5, height = 5)






# Is beta (dissimilarities in communities) diversity different for pits+biochar vs. pits?
# NMDS with bray curtis on relative abundances
#pitBiochar.ord <- ordinate(pitBiochar, "NMDS", "bray") #stress = 0.038
# PCoA with bray curtis on relative abundances
pitBiochar.ord <- ordinate(pitBiochar, "PCoA", "bray") 


# centroids for NMDS
pitBiochar.ord$points # gives me the x and y coordinates for each 
scores(pitBiochar.ord)$sites # gives me the same thing
# centroids are the average value
pitBiochar.ord.scores <- as(scores(pitBiochar.ord)$sites,"matrix")
pitBiochar.ord.scores <- as.data.frame(pitBiochar.ord.scores) #rows are samples and columns are sample data
# add this to the mapping dataframe
mapping <- as(sample_data(pitBiochar), "matrix")
mapping <- merge(mapping, pitBiochar.ord.scores, by = 0)
# calculate the centroids
pitBiochar.ord.scores.centroids <- mapping %>% group_by(bulk_soil) %>% summarise(NMDS1_ave = mean(NMDS1),NMDS2_ave = mean(NMDS2)) # plot the centroids of the reference biocrusts onto every panel of the figure 

# centroids for the PCoA
mapping <- as(sample_data(pitBiochar), "matrix")
axis_values <- pitBiochar.ord$vectors[ ,1:2]
mapping <- merge(mapping, axis_values, by = 0)
# calculate the centroids
pitBiochar.ord.scores.centroids <- mapping %>% group_by(bulk_soil) %>% summarise(pcoa_axis1_ave = mean(Axis.1),pcoa_axis2_ave = mean(Axis.2)) # plot the centroids of the reference biocrusts onto every panel of the figure 

# Plotting
pitBiochar@sam_data$bulk_soil <- factor(pitBiochar@sam_data$bulk_soil, levels = c("C", "P","BP"))

ordination_DF <- plot_ordination(pitBiochar, pitBiochar.ord, justDF = TRUE)

ordination_DF %>%
  ggplot(aes(x = ordination_DF[,1], y = ordination_DF[,2], color = bulk_soil))+
  geom_point(size = 5)+ # alpha = 0.7
  scale_color_manual(values = c("orange4","darkslategray", "goldenrod"), labels = c("Control", "Pit", "Pit + Biochar")) + 
  labs(x = "Axis 1 [28.2%]", y = "Axis 2 [22.8%]", color = "Microsite\nTreatment")+
  theme_classic()+
  theme(strip.text = element_text(size=14),
        axis.text = element_text(size = 14),
        axis.title = element_text(face = "bold", size = 15),
        legend.text = element_text(size = 14),
        legend.title = element_text(face = "bold", size = 14),
        legend.position = "right") + 
  guides(color = guide_legend(ncol = 1))

#ggsave("output/PCoA_bray_relAbund.pdf", width = 7, height = 5)

# Statistical analysis for ordination
bray_dist <- phyloseq::distance(pitBiochar, method="bray", weighted= F)

# PERMANOVA - test for differences in centroid of different pCoA hull groups
adon <- adonis2(bray_dist ~ sample_data(pitBiochar)$bulk_soil, method = "bray") # centroid location not significantly different for controls vs. pits 
adon # significant 
pairwise.perm.manova(bray_dist, sample_data(pitBiochar)$bulk_soil, test = "Wilks", p.method = "holm", F = TRUE, R2 = TRUE) # not significant
#library(EcolUtils) # alternative option for pairwise contrasts
#adonis.pair(bray_dist, sample_data(pitBiochar)$bulk_soil) # not significant

# PERMDISP - multivariate version of Levenne Test. Difference in dispersion for each factor level
m1 <- betadisper(bray_dist, sample_data(pitBiochar)$bulk_soil)
anova(m1) # Dispersion not different (p = 0.5466)
# if permdisp were significant, then you would want to use a different test than PERMANOVA to compare the groups
#tukeyHSD(m1) # post hoc test
#scores(m1) # access the scores







# Taxonomic Analyses
# MULTIPATT doesn't give a value for every OTU. It only outputs important OTUs
set.seed(1223) # For reproducibility
OTU_intact <- as.data.frame(pitBiochar@otu_table)
OTU_noZero <- subset(OTU_intact, rowSums(OTU_intact) > 0)
# factor the bulk_soil levels
pitBiochar@sam_data$bulk_soil <- factor(pitBiochar@sam_data$bulk_soil, levels = c("P","C", "BP"))
mp <- multipatt(t(OTU_noZero), 
                sample_data(pitBiochar)$bulk_soil, 
                #func = "IndVal.g", 
                control = how(nperm=999))
summary(mp, alpha = 1) # Pit treatment doesn't show because there are no significant ones

# save the output
mp_results <- mp$sign # this one is better
# remove rows with significant p-value
# initially the list of ASVs is 1278 rows. Significant p-values results in 119 rows of data
mp_results <- mp_results %>% filter(p.value < 0.05)

# subset into dataframes based on bulk soil - only if there is a 1 in C and no where else. This is keeping the ones that are only indicators for a single desert
multipatt_C <- mp_results %>% filter(s.C == 1 & s.P == 0 & s.BP == 0) # 17 significant indicators for the control treatment 
multipatt_P <- mp_results %>% filter(s.C == 0 & s.P == 1 & s.BP == 0) # no significant indicators for the pit treatment 
multipatt_BP <- mp_results %>% filter(s.C == 0 & s.P == 0 & s.BP == 1) # 46 significant indicators for the pit + biochar treatment 
multipatt_microsite <- mp_results %>% filter(s.C == 0 & s.P == 1 & s.BP == 1) # 46 significant indicators for the pit + biochar treatment 

# make a dataframe 
multipatt_C$bulk_soil <- "C"
#multipatt_P$bulk_soil <- "P" # not needed because there are none
multipatt_BP$bulk_soil <- "BP"
multipatt_microsite$bulk_soil <- "microsite"
# transfer the row names to a column
multipatt_C$ASV <- row.names(multipatt_C)
multipatt_BP$ASV <- row.names(multipatt_BP)
multipatt_microsite$ASV <- row.names(multipatt_microsite)
# merge 
multipatt_results <- rbind(multipatt_C, multipatt_BP, multipatt_microsite)

# merge the taxonomy to the multipatt results dataframe
pitBiochar_taxonomy <- as.data.frame(pitBiochar@tax_table)
pitBiochar_taxonomy$ASV <- row.names(pitBiochar_taxonomy)
multipatt_results <- left_join(multipatt_results, pitBiochar_taxonomy, by = "ASV")

# export this list for the supplementary data
#write.csv(multipatt_results , "output/indicator_list_taxonomy.csv")

# make a visualization of this?
# mp_results_summary <- multipatt_results %>% group_by(bulk_soil, Phylum) %>% summarise(count = n())
# # multiple bacteroidota are enhanced with BP
# multipatt_results %>% filter(bulk_soil == "BP") %>% filter(Phylum == "Bacteroidota")
# multipatt_results %>% filter(bulk_soil == "BP") %>% filter(Phylum == "Firmicutes")
# multipatt_results %>% filter(bulk_soil == "BP") %>% filter(Phylum == "Proteobacteria")
# multipatt_results %>% filter(bulk_soil == "BP") %>% filter(Phylum == "Verrucomicrobiota")
# multipatt_results %>% filter(bulk_soil == "BP") %>% filter(Phylum == "Acidobacteriota")
# multipatt_results %>% filter(bulk_soil == "BP") %>% filter(Phylum == "Chloroflexi")
# multipatt_results %>% filter(bulk_soil == "BP") %>% filter(Phylum == "Crenarchaeota")
# 
# multipatt_results %>% filter(bulk_soil == "C") %>% filter(Phylum == "Actinobacteriota")
# multipatt_results %>% filter(bulk_soil == "C") %>% filter(Phylum == "Chloroflexi")
# multipatt_results %>% filter(bulk_soil == "C") %>% filter(Phylum == "Planctomycetota")

#modified columns to keep taxonomy information
multipatt_results$k_Kingdom <- paste("k_", multipatt_results$Kingdom)
multipatt_results$p_Phylum <- paste("p_", multipatt_results$Phylum)
multipatt_results$c_Class <- paste("c_", multipatt_results$Class)
multipatt_results$o_Order <- paste("o_", multipatt_results$Order)
multipatt_results$f_Family <- paste("f_", multipatt_results$Family)
multipatt_results$g_Genus <- paste("g_", multipatt_results$Genus)

# make a new column 
multipatt_results$taxonomy_k_ASV <- paste(multipatt_results$k_Kingdom, multipatt_results$p_Phylum, multipatt_results$c_Class, multipatt_results$o_Order, multipatt_results$f_Family, multipatt_results$g_Genus, multipatt_results$ASV, sep = ";")

# make a new column without phylym 
multipatt_results$taxonomy_c_ASV <- paste(multipatt_results$c_Class, multipatt_results$o_Order, multipatt_results$f_Family, multipatt_results$g_Genus, multipatt_results$ASV, sep = ";")

# save the list of indicator taxa
#write.csv(multipatt_results, "output/indicator_list_taxonomy.csv")









# Plotting Phyla Relative Abundance
#Get count of reads by phyla
table(phyloseq::tax_table(pitBiochar)[, "Phylum"])
# get the total read count per sample
reads_per_sample <- as.data.frame (sample_sums(pitBiochar)) # there are different values here when they should be all the same based on the rarefaction process!!!! 
#copy the rownames to a column
reads_per_sample <- tibble::rownames_to_column(reads_per_sample, "Sample")
# change the name of the second column
colnames(reads_per_sample) <- c("Sample", "totalReads")
nrow(reads_per_sample) # 6 samples are included

# Think about the phyla 
get_taxa_unique(pitBiochar, "Phylum") # 28 unique Phyla
# I only want to plot the top 9 and the rest go into an Other category

# first melt the data out of phyloseq format and into long dataframe
pits_psmelt <- psmelt(pitBiochar)

# Get a list of the most abundant Phyla
temp <- pits_psmelt %>%
  group_by(Phylum) %>%
  summarise(totalAbund = sum(Abundance)) %>%
  arrange(-totalAbund)
temp <- temp[1:9,]

top9Phyla <- temp$Phylum # list that holds the Top9Phyla names

#subset the psmelt dataframe for only the Top9Phyla
phylaSubset <- pits_psmelt %>%
  filter(Phylum %in% top9Phyla)
# check that we got the right ones
unique(phylaSubset$Phylum) # yes

# calculate the reads per Phylum, grouped by Sample
phylaSubset_abund <- phylaSubset %>%
  group_by(Sample, Phylum) %>%
  summarise(abund = sum(Abundance))

#make a new column were I divide each of the values by the correct one in the reads_per_sample list
# first merge the reads_per_sample data onto the dataframe
phylaSubset_abund2 <- inner_join(phylaSubset_abund, reads_per_sample, by = "Sample")
# divide the abund column by the total reads column
phylaSubset_abund2$relAbund <- phylaSubset_abund2$abund / phylaSubset_abund2$totalReads
# check that it does not add to 1 because there should be some in the Other category
total_relAbund <- phylaSubset_abund2 %>%
  group_by(Sample) %>%
  summarise(sum = sum(relAbund))
# add a new column which contains the relative abundance value for the Other category
total_relAbund$other <- (1 - total_relAbund$sum)
# delete the sum column
total_relAbund <- total_relAbund %>%
  select(Sample, other)
#add column which is "Other" repeated
total_relAbund$Phylum <- "Other"
#change column header "other" to relAbund
colnames(total_relAbund)[which(names(total_relAbund) == "other")] <- "relAbund"
# select columns to keep in the dataframe we want
phylaSubset_abund2 <- phylaSubset_abund2 %>%
  select(Sample, relAbund, Phylum)
# rbind the other values to the phylaSubset_abund2 dataframe
phylaSubset_abund3 <- rbind(phylaSubset_abund2, total_relAbund)


# plot it!
# phylaSubset_abund3 %>%
#   ggplot(aes(y = relAbund, x = Sample, fill = Phylum)) +
#   geom_bar(stat = "identity", colour = "black", linewidth = 0.25) +
#   labs(x = "Sample", y = "Relative Abundance", fill = "Phyla") +
#   scale_fill_brewer(palette = "Paired") +
#   theme_classic()+
#   theme(axis.title = element_text(face = "bold", size = 16),
#         axis.text.x = element_text(size = 14, angle=90),
#         legend.text = element_text(face = "italic"))

# plot it! new colors and saved for ggarrange plot
#phylaSubset_abund3$bulk_soil <- factor(phylaSubset_abund3$bulk_soil, levels = c("C","P","B"))

phylaSubset_abund3 %>%
  ggplot(aes(y = relAbund, x = Sample, fill = Phylum)) +
  geom_bar(stat = "identity", colour = "black", linewidth = 0.25) +
  labs(x = "Sample", y = "Relative Abundance", fill = "Phyla") +
  scale_fill_manual(values = c("brown4", "mediumpurple4","bisque3","cadetblue3","darkslategray","goldenrod","lightsteelblue4","coral2", "cornflowerblue", "cornsilk"))+
  scale_x_discrete(labels = c("","Control","","","Pit","","","Pit+Biochar",""))+
  theme_classic()+
  theme(axis.title = element_text(face = "bold", size = 15),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(face = "italic", size = 14),
        axis.ticks.x = element_blank(),
        legend.title=element_text(face = "bold",size=14))

#ggsave("output/relAbund_barplot.pdf", width = 7, height = 5)





#Compare Phyla Relative Abundances
phylaSubset_abund3$bulk_soil <- substr(phylaSubset_abund3$Sample, 1, 1)  
phylaSubset_abund3$Phy_bulk_soil <- paste(phylaSubset_abund3$Phylum, phylaSubset_abund3$bulk_soil)
# CLiff's Tutorial uses a Kruskal-Wallis test (multiple phyla at once)
kruskal.test(relAbund ~ Phy_bulk_soil, data=phylaSubset_abund3)
# As the p-value is less than the significance level 0.05, we can conclude that there are significant differences between the treatment groups.
pairwise.wilcox.test(phylaSubset_abund3$relAbund, phylaSubset_abund3$Phy_bulk_soil,p.adjust.method = "BH")
# none seem different with this method even though the overall kruskall test is significant



# Compare the top 9 phyla between groups
# These may be wrong because the data is relative abundance and so all the values are connected to one another. I also think I need a p-value correction if I am running multiple tests in a row???
#
acidobacteria <- phylaSubset_abund3 %>%
  filter(Phylum == "Acidobacteriota")
#acidobacteria.aov <- aov(acidobacteria$relAbund ~ acidobacteria$bulk_soil) 
#shapiro.test(acidobacteria.aov$residuals) # Residuals normally distributed (p > 0.05)
#qqnorm(acidobacteria.aov$residuals)
#qqline(acidobacteria.aov$residuals, col = "red")
#bartlett.test(relAbund ~ bulk_soil, data=acidobacteria) # samples have equal variance
#plot(acidobacteria.aov) # Console for plots
#summary(acidobacteria.aov) #significant
#TukeyHSD(acidobacteria.aov) # controls are significantly different (higher)
# Kruskal-test if assumptions are violated or simply use it since it is compositional data
kruskal.test(relAbund ~ bulk_soil, data = acidobacteria) # not significant

#
actinobacteria <- phylaSubset_abund3 %>%
  filter(Phylum == "Actinobacteriota")
#actinobacteria.aov <- aov(actinobacteria$relAbund ~ actinobacteria$bulk_soil) 
#shapiro.test(actinobacteria.aov$residuals) # Residuals normally distributed (p < 0.05)
#qqnorm(actinobacteria.aov$residuals)
#qqline(actinobacteria.aov$residuals, col = "red")
#bartlett.test(relAbund ~ bulk_soil, data=actinobacteria) # samples have equal variance
#plot(actinobacteria.aov) # Console for plots
#summary(actinobacteria.aov) #significant
#TukeyHSD(actinobacteria.aov) # controls are significantly different (higher)
# Kruskal-test instead
kruskal.test(relAbund ~ bulk_soil, data = actinobacteria) # not significant

#
bacteroidota <- phylaSubset_abund3 %>%
  filter(Phylum == "Bacteroidota")
#bacteroidota.aov <- aov(bacteroidota$relAbund ~ bacteroidota$bulk_soil) 
#shapiro.test(bacteroidota.aov$residuals) # Residuals normally distributed (p < 0.05)
#qqnorm(bacteroidota.aov$residuals)
#qqline(bacteroidota.aov$residuals, col = "red") # looks bad
#bartlett.test(relAbund ~ bulk_soil, data=bacteroidota) # samples do not have equal variance
#plot(bacteroidota.aov) # Console for plots
#summary(bacteroidota.aov) #significant but not valid
#TukeyHSD(bacteroidota.aov) # controls are significantly different (lower)
# Kruskal-test instead
kruskal.test(relAbund ~ bulk_soil, data = bacteroidota) # not significant
#
chloroflexi <- phylaSubset_abund3 %>%
  filter(Phylum == "Chloroflexi")
#chloroflexi.aov <- aov(chloroflexi$relAbund ~ chloroflexi$bulk_soil) 
#shapiro.test(chloroflexi.aov$residuals) # Residuals normally distributed (p < 0.05)
#qqnorm(chloroflexi.aov$residuals)
#qqline(chloroflexi.aov$residuals, col = "red") 
#bartlett.test(relAbund ~ bulk_soil, data=chloroflexi) # samples have equal variance
#plot(chloroflexi.aov) # Console for plots
#summary(chloroflexi.aov) #significant
#TukeyHSD(chloroflexi.aov) # controls are significantly different (higher)
# Kruskal-test instead
kruskal.test(relAbund ~ bulk_soil, data = chloroflexi) # not significant
#
crenarchaeota <- phylaSubset_abund3 %>%
  filter(Phylum == "Crenarchaeota")
#crenarchaeota.aov <- aov(crenarchaeota$relAbund ~ crenarchaeota$bulk_soil) 
#shapiro.test(crenarchaeota.aov$residuals) # Residuals normally distributed (p < 0.05)
#qqnorm(crenarchaeota.aov$residuals)
#qqline(crenarchaeota.aov$residuals, col = "red") 
#bartlett.test(relAbund ~ bulk_soil, data=crenarchaeota) # samples have equal variance
#plot(crenarchaeota.aov) # Console for plots
#summary(crenarchaeota.aov) #significant
#TukeyHSD(crenarchaeota.aov) # controls are significantly different (higher)
# Kruskal-test instead
kruskal.test(relAbund ~ bulk_soil, data = crenarchaeota) # significant
dunnTest(relAbund ~ bulk_soil, data = crenarchaeota, method = "bonferroni")

#
cyanobacteria <- phylaSubset_abund3 %>%
  filter(Phylum == "Cyanobacteria")
#cyanobacteria.aov <- aov(cyanobacteria$relAbund ~ cyanobacteria$bulk_soil) 
#shapiro.test(cyanobacteria.aov$residuals) # Residuals normally distributed (p < 0.05)
#qqnorm(cyanobacteria.aov$residuals)
#qqline(cyanobacteria.aov$residuals, col = "red") 
#bartlett.test(relAbund ~ bulk_soil, data=cyanobacteria) # samples have equal variance
#plot(cyanobacteria.aov) # Console for plots
#summary(cyanobacteria.aov) # not significant
#TukeyHSD(cyanobacteria.aov) 
# Kruskal-test instead
kruskal.test(relAbund ~ bulk_soil, data = cyanobacteria) # not significant
#
planctomycetota <- phylaSubset_abund3 %>%
  filter(Phylum == "Planctomycetota")
#planctomycetota.aov <- aov(planctomycetota$relAbund ~ planctomycetota$bulk_soil) 
#shapiro.test(planctomycetota.aov$residuals) # Residuals normally distributed (p < 0.05)
#qqnorm(planctomycetota.aov$residuals)
#qqline(planctomycetota.aov$residuals, col = "red") 
#bartlett.test(relAbund ~ bulk_soil, data=planctomycetota) # samples have equal variance
#plot(planctomycetota.aov) # Console for plots
#summary(planctomycetota.aov) #significant
#TukeyHSD(planctomycetota.aov) # controls are significantly different (higher)
# Kruskal-test instead
kruskal.test(relAbund ~ bulk_soil, data = planctomycetota) # significant
dunnTest(relAbund ~ bulk_soil, data = planctomycetota, method = "bonferroni")

#
proteobacteria <- phylaSubset_abund3 %>%
  filter(Phylum == "Proteobacteria")
#proteobacteria.aov <- aov(proteobacteria$relAbund ~ proteobacteria$bulk_soil) 
#shapiro.test(proteobacteria.aov$residuals) # Residuals normally distributed (p < 0.05)
#qqnorm(proteobacteria.aov$residuals)
#qqline(proteobacteria.aov$residuals, col = "red") 
#bartlett.test(relAbund ~ bulk_soil, data=proteobacteria) # samples do not have equal variance
#plot(proteobacteria.aov) # Console for plots
#summary(proteobacteria.aov) # significant but not valid
#TukeyHSD(proteobacteria.aov) # pits are significantly different (higher)
# Kruskal-test instead
kruskal.test(relAbund ~ bulk_soil, data = proteobacteria) # significant
dunnTest(relAbund ~ bulk_soil, data = proteobacteria, method = "bonferroni")

#
verrucomicrobiota <- phylaSubset_abund3 %>%
  filter(Phylum == "Verrucomicrobiota")
#verrucomicrobiota.aov <- aov(verrucomicrobiota$relAbund ~ verrucomicrobiota$bulk_soil) 
#shapiro.test(verrucomicrobiota.aov$residuals) # Residuals normally distributed (p < 0.05)
#qqnorm(verrucomicrobiota.aov$residuals)
#qqline(verrucomicrobiota.aov$residuals, col = "red") 
#bartlett.test(relAbund ~ bulk_soil, data=verrucomicrobiota) # samples have equal variance
#plot(verrucomicrobiota.aov) # Console for plots
#summary(verrucomicrobiota.aov) # significant
#TukeyHSD(verrucomicrobiota.aov) # biochar is significantly different
# Kruskal-test instead
kruskal.test(relAbund ~ bulk_soil, data = verrucomicrobiota) # not significant

# Plot it 
# change the levels of Phyla
gotcha <- phylaSubset_abund3 %>%
  group_by(Phylum) %>%
  summarise(meanrelAbund = mean(relAbund))
gotcha <- gotcha %>%
  arrange(-meanrelAbund)
top10phyla <- gotcha$Phylum

phylaSubset_abund3$Phylum <- factor(phylaSubset_abund3$Phylum, levels = c("Actinobacteriota", "Bacteroidota" ,"Proteobacteria","Chloroflexi","Crenarchaeota","Other","Acidobacteriota","Cyanobacteria","Planctomycetota","Verrucomicrobiota"))

phylaSubset_abund3 <- phylaSubset_abund3 %>%
  mutate(bulk_soil = case_when(
    startsWith(Sample, "B") ~ "BP",
    startsWith(Sample, "C") ~ "C",
    startsWith(Sample, "P") ~ "P",
  ))

phylaSubset_abund3$bulk_soil <- factor(phylaSubset_abund3$bulk_soil, levels = c("C","P","BP"))

# Save this plot
ggplot(data = phylaSubset_abund3, aes(x = bulk_soil, y = relAbund)) +
  geom_boxplot(outlier.shape  = NA, aes(fill = bulk_soil), alpha = 0.7) +
  geom_jitter(height = 0, width = .2, alpha = 0.5, aes(color = bulk_soil)) +
  labs(x = "", y = "Relative Abundance") +
  facet_wrap(~ Phylum, scales = "free")+
  scale_fill_manual(values = c("orange4","darkslategray", "goldenrod"), labels = c("Control","Pit", "Pit+Biochar")) +
  scale_color_manual(values = c("orange4","darkslategray", "goldenrod"), guide = "none") +
  #scale_x_discrete(labels = c("Control","Pit", "Pit+Biochar")) +
  labs(fill = "Treatment")+
  ylim(0,NA)+
  theme_classic() +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(size=12), 
        axis.title.y = element_text(size=14, face = "bold"), 
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())

#ggsave("output/phyla_relAbund_boxplots.pdf",width = 10, height = 6)


# note that some of the differences are between control and microsite...so different statistical test: 
phylaSubset_abund3$microsite <- "microsite"
phylaSubset_abund3$microsite[phylaSubset_abund3$microsite == 'Orange St'] <- 'Portola Pkwy'
phylaSubset_abund3 <- phylaSubset_abund3 %>%
  mutate(microsite = case_when(
    # When "self manag" is found, fill Plan with "S"
    grepl(pattern = "C", x = bulk_soil) ~ "control",
    # Same, with "P"
    grepl(pattern = "BP", x = bulk_soil) ~ "microsite",
    grepl(pattern = "P", x = bulk_soil) ~ "microsite"
  ))

# Now try t.tests
acidobacteria <- phylaSubset_abund3 %>%
  filter(Phylum == "Acidobacteriota")
wilcox.test(acidobacteria$relAbund ~ acidobacteria$microsite)
# W = 18, p-value = 0.02381

actinobacteria <- phylaSubset_abund3 %>%
  filter(Phylum == "Actinobacteriota")
wilcox.test(actinobacteria$relAbund ~ actinobacteria$microsite)
# W = 18, p-value = 0.02381

bacteroidota <- phylaSubset_abund3 %>%
  filter(Phylum == "Bacteroidota")
wilcox.test(bacteroidota$relAbund ~ bacteroidota$microsite)
# W = 0, p-value = 0.02381

chloroflexi <- phylaSubset_abund3 %>%
  filter(Phylum == "Chloroflexi")
wilcox.test(chloroflexi$relAbund ~ chloroflexi$microsite)
# W = 18, p-value = 0.02381

crenarchaeota <- phylaSubset_abund3 %>%
  filter(Phylum == "Crenarchaeota")
wilcox.test(crenarchaeota$relAbund ~ crenarchaeota$microsite)
# W = 18, p-value = 0.02381

cyanobacteria <- phylaSubset_abund3 %>%
  filter(Phylum == "Cyanobacteria")
wilcox.test(cyanobacteria$relAbund ~ cyanobacteria$microsite)
# W = 0, p-value = 0.02381

planctomycetota <- phylaSubset_abund3 %>%
  filter(Phylum == "Planctomycetota")
wilcox.test(planctomycetota$relAbund ~ planctomycetota$microsite)
# W = 18, p-value = 0.02381

proteobacteria <- phylaSubset_abund3 %>%
  filter(Phylum == "Proteobacteria")
wilcox.test(proteobacteria$relAbund ~ proteobacteria$microsite)
# W = 0, p-value = 0.02381

verrucomicrobiota <- phylaSubset_abund3 %>%
  filter(Phylum == "Verrucomicrobiota")
wilcox.test(verrucomicrobiota$relAbund ~ verrucomicrobiota$microsite)
# W = 2, p-value = 0.09524

# Linear models may be better if kruskal wallis is not satistfactory
# fixed effect is treatment type (bulk_soil)
# fixed effect biochar
#phylaSubset_abund3$microsite <- "m"
#phylaSubset_abund3 <- phylaSubset_abund3 %>%
#  mutate(microsite = case_when(
#  grepl(pattern = "C", x = bulk_soil) ~ "c",
#  TRUE ~ microsite
#  ))
#phylaSubset_abund3$biochar <- "no_biochar"
#phylaSubset_abund3$biochar[phylaSubset_abund3$bulk_soil == "BP"] <- "biochar"

#models - would need to check each one for parametric assumptions. 
# actino.lm <- lm(relAbund ~ bulk_soil, data = actinobacteria)
# summary(actino.lm)
# plot(actino.lm, 1) # idk
# hist(actino.lm$residuals) # not normal
# plot(actino.lm, 2) # looks good
# bacter.lm <- lm(relAbund ~ bulk_soil, data = bacteroidota)
# summary(bacter.lm)
# proteo.lm <- lm(relAbund ~ bulk_soil, data = proteobacteria)
# summary(proteo.lm)
# chloro.lm <- lm(relAbund ~ bulk_soil, data = chloroflexi)
# summary(chloro.lm)
# cren.lm <- lm(relAbund ~ bulk_soil, data = crenarchaeota)
# summary(cren.lm)
# acido.lm <- lm(relAbund ~ bulk_soil, data = acidobacteria)
# summary(acido.lm)
# cyano.lm <- lm(relAbund ~ bulk_soil, data = cyanobacteria)
# summary(cyano.lm)
# plancto.lm <- lm(relAbund ~ bulk_soil, data = planctomycetota)
# summary(plancto.lm)
# verru.lm <- lm(relAbund ~ bulk_soil, data = verrucomicrobiota)
# summary(verru.lm)







### Plot one bar for each sample type instead of three
# first calculate the averages
phylaSubset_abund3_aggregate <- phylaSubset_abund3 %>% 
  group_by(bulk_soil, Phylum) %>%
  summarise(meanrelAbund = mean(relAbund),
            sdrelAbund = sd(relAbund))

# plot: Absolute abundance phylum level for supplementary materials
phylaSubset_abund3_aggregate$Phylum <- factor(phylaSubset_abund3_aggregate$Phylum, levels = c("Actinobacteriota", "Bacteroidota", "Proteobacteria", "Chloroflexi", "Crenarchaeota", "Acidobacteriota", "Cyanobacteria", "Planctomycetota", "Verrucomicrobiota", "Other"))


phylaSubset_abund3_aggregate %>%
  ggplot(aes(y = meanrelAbund, x = bulk_soil, fill = Phylum)) +
  geom_bar(stat = "identity", colour = "black", linewidth = 0.25) +
  labs(x = "", y = "Relative Abundance", fill = "Phylum") +
  scale_fill_manual(values = c("brown4", "mediumpurple4","bisque3","cadetblue3","darkslategray","goldenrod","lightsteelblue4","coral2", "cornflowerblue", "cornsilk"))+
  scale_x_discrete(labels = c("Control","Pit","Pit+Biochar"))+
  theme_classic()+
  theme(axis.title = element_text(face = "bold", size = 15),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(face = "italic", size = 14),
        axis.ticks.x = element_blank(),
        legend.title=element_text(face = "bold",size=14))

#ggsave("output/relAbund_barplot.pdf", width = 7, height = 5)










### Absolute Abundances ### 
# Note: when writing about this, just say 16S rRNA copies to avoid confusion around # of 16S rRNA copies per cell
### Add in the qPCR data ###
phylaSubset_abund3 <- left_join(phylaSubset_abund3, qPCR, by = c("Sample" = "Sample_ID"))
phylaSubset_abund3$absolute_abundance <- ((phylaSubset_abund3$relAbund * phylaSubset_abund3$qPCR_count)) / phylaSubset_abund3$Mass_extracted_g
# the math I did to get number of 16S rRNA copies x 10^6 x cm^-2 of soil
# I found the average weight of each sample = 11 g
# we took 7 cores to get that amount, so each core weighs 1.57 g
# each core had a diameter of 1 cm, so the area of one core is 0.785
# we only extracted 0.25g of soil, so when I convert, I get 0.125 cm^2 that was extracted for total DNA
# this means I should divide the absolute abundance column by 0.125 cm^2 to get the number of copies per area
#phylaSubset_abund3$abso_abund_cm2 <- phylaSubset_abund3$absolute_abundance / 0.125
# per gram of soil
phylaSubset_abund3$abso_abund <- phylaSubset_abund3$absolute_abundance / (10^6)

# All together but average (not every replicate)
# first calculate the averages
phylaSubset_abund3_aggregate <- phylaSubset_abund3 %>% 
  group_by(bulk_soil, Phylum) %>%
  summarise(meanAbso_abund = mean(abso_abund),
            sdAbso_abund = sd(abso_abund))

# plot: Absolute abundance phylum level for supplementary materials
phylaSubset_abund3_aggregate$Phylum <- factor(phylaSubset_abund3_aggregate$Phylum, levels = c("Other", "Verrucomicrobiota", "Planctomycetota", "Cyanobacteria", "Acidobacteriota", "Crenarchaeota", "Chloroflexi", "Proteobacteria", "Bacteroidota", "Actinobacteriota"))

new_labels <- c("BP" = "Pit + Biochar", "P" = "Pit","C" = "Degraded Soil")

phylaSubset_abund3_aggregate %>%
  ggplot(aes(y = meanAbso_abund, x = bulk_soil, fill = Phylum)) +
  geom_bar(stat = "identity", alpha = 0.8, width = 0.85) +
  scale_fill_manual(values = grafify_muted_greenstart, breaks=c("Actinobacteriota", "Bacteroidota" ,"Proteobacteria","Chloroflexi","Crenarchaeota","Acidobacteriota","Cyanobacteria","Planctomycetota","Verrucomicrobiota", "Other")) +
  scale_y_continuous(expand = c(0, 0.1), breaks=seq(0, 50, 10)) +
  theme_classic()+
  coord_flip()+
  facet_grid(rows = vars(bulk_soil), scales = "free", space = "free", switch = "y",labeller = labeller(bulk_soil = new_labels))+
  labs(x = "", y = expression(paste("No. 16S rRNA gene copies x ", 10^6, " x "  ,g^-1,)), fill = "")+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size = 12),
        #strip.background = element_rect(fill="white", colour=NA),
        strip.background = element_blank(),
        #strip.text = element_blank(),
        strip.text = element_text(face="bold", size=12),
        strip.text.y.left = element_text(angle = 0),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 10),
        plot.margin = margin(0.5,1,0,0, "cm"))+
  guides(fill = guide_legend(nrow = 3, byrow = TRUE))

#ggsave("output/phylaabso_abund_stackedbar.pdf", width = 6.5, height = 3.5)

# calculate more information for the manuscript
absoAbund_minmax <- phylaSubset_abund3 %>%
  #filter(Phylum == "Cyanobacteria") %>%
  #filter(Phylum %in% c("Cyanobacteria","Bacteroidota","Chloroflexi","Acidobacteriota")) %>%
  group_by(bulk_soil, Phylum) %>%
  summarise(meanAbso = mean(abso_abund),
            sdAbso = sd(abso_abund),
            meanRel = mean(relAbund),
            sdRel = sd(relAbund))


## I STOPPED HERE on 9/12/23




# I want to know the taxonomic information for the ASVs that were detected in the pits+biochar and not in the controls or pits
# get a list of the ASVs in pits+biochar
melt_pitBiochar <- psmelt(pitBiochar)
BP <- melt_pitBiochar %>%
  filter(bulk_soil == "BP" & Abundance > 0)
BP_ASVs <- data.frame(unique(BP$OTU))
# get a list of the ASVs in controls
C <- melt_pitBiochar %>%
  filter(bulk_soil == "C" & Abundance > 0)
C_ASVs <- data.frame(unique(C$OTU))
# filter the Pit+Biochar treatment for the ASVs that are not found in the controls
BP_uniqueASVs <- BP_ASVs %>%
  filter(!unique.BP.OTU. %in% C_ASVs$unique.C.OTU.) # there are 196 ASVs
colnames(BP_uniqueASVs) <- "OTU"
# now take a look at their taxonomy
temp <- inner_join(BP, BP_uniqueASVs, by = "OTU")
#get rid of extra columns
temp <- temp %>%
  select("OTU", "Sample","Abundance","Kingdom","Phylum","Class","Order","Family","Genus")
# graph these at Phylum level
temp2 <- temp %>%
  group_by(Sample,Phylum) %>%
  summarise(abund = sum(Abundance),
            n=n())
# this plot shows me the total number of reads OR the number of ASVs associated with each phylum which was unique to pit+biochar when comparing that treatment to controls. 
temp2 %>%
  ggplot()+
  geom_col(mapping = aes(x=Phylum, y = n, fill = Sample), position = position_dodge(preserve = "single"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90)) # I see that the majority of reads that differ between pits+biochar and controls are from the Bacteroidota and Proteobacteria. 

### Now do the same thing comparing pits+biochar to pits only 
# get a list of the ASVs in pits
P <- melt_pitBiochar %>%
  filter(bulk_soil == "P" & Abundance > 0)
P_ASVs <- data.frame(unique(P$OTU))
# filter the Pit+Biochar treatment for the ASVs that are not found in the pits
BP_uniqueASVs2 <- BP_ASVs %>%
  filter(!unique.BP.OTU. %in% P_ASVs$unique.P.OTU.) # there are 196 ASVs
colnames(BP_uniqueASVs2) <- "OTU"
# now take a look at their taxonomy
temp <- inner_join(BP, BP_uniqueASVs2, by = "OTU")
#get rid of extra columns
temp <- temp %>%
  select("OTU", "Sample","Abundance","Kingdom","Phylum","Class","Order","Family","Genus")
# graph these at Phylum level
temp2 <- temp %>%
  group_by(Sample,Phylum) %>%
  summarise(abund = sum(Abundance),
            n=n())
# this plot shows me the total number of reads OR the number of ASVs associated with each phylum which was unique to pit+biochar when comparing that treatment to controls. 
temp2 %>%
  ggplot()+
  geom_col(mapping = aes(x=Phylum, y = n, fill = Sample), position = position_dodge(preserve = "single"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90)) # I see that the majority of reads that differ between pits+biochar and pits are from the phyla Acidobacteriota, Actinobacteriota, Crenarchaeota, Proteobacteria. 

## Now I am interested in the genera associated with Bacteroidota in each treatment category
bacteroidota_abundances <- melt_pitBiochar %>%
  filter(Phylum == "Bacteroidota")
unique(bacteroidota_abundances$Genus) # there are 32 unique genera,
unique(bacteroidota_abundances$OTU) #  there are 232 unique ASVs

# number of unique ASVs by sample 
temp <- bacteroidota_abundances %>%
  group_by(Sample,Genus) %>%
  summarise(totalAbund = sum(Abundance),
            n = n())

temp$bulk_soil <- substr(temp$Sample, 1,1)

output <- temp %>%
  group_by(bulk_soil, Genus) %>%
  summarise(meantotalAbund = mean(totalAbund),
            sdtotalAbund = sd(totalAbund))

#simple graph
output <- output %>%
  filter(meantotalAbund > 0)
output %>%
  ggplot()+
  geom_point(mapping = aes( x=Genus, y = meantotalAbund, color = bulk_soil))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))

#boxplot graph of the genera of interest: Flaviaesturariibacter, Cnuella, Rhodocytophaga, Segetibacter, Adhaeribacter, Flavisolibacter, Pedobacter, Nibribacter
bacteroidota_subset <- temp %>% filter(Genus %in% c("Flaviaesturariibacter", "Cnuella", "Rhodocytophaga", "Segetibacter", "Adhaeribacter", "Flavisolibacter", "Pedobacter", "Nibribacter"))

bacteroidota_subset$bulk_soil <- factor(bacteroidota_subset$bulk_soil, levels = c("C","P","B"))

ggplot(data = bacteroidota_subset, aes(x = bulk_soil, y = totalAbund)) +
  geom_boxplot(outlier.shape  = NA, aes(fill = bulk_soil), alpha = 0.7) +
  geom_jitter(height = 0, width = .2, alpha = 0.5, aes(color = bulk_soil)) +
  labs(x = "", y = "Relative Abundance") +
  facet_wrap(~ Genus, scales = "free")+
  scale_fill_manual(values = c("orange4","darkslategray", "goldenrod"), labels = c("Control","Pit", "Pit+Biochar")) +
  scale_color_manual(values = c("orange4","darkslategray", "goldenrod"), guide = "none") +
  #scale_x_discrete(labels = c("Control","Pit", "Pit+Biochar")) +
  labs(fill = "Treatment")+
  ylim(0,NA)+
  theme_classic() +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(size=12), 
        axis.title.y = element_text(size=14, face = "bold"), 
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())

#ggsave("output/bacteroidota_relAbund_boxplots.pdf",width = 10, height = 6)



## Now I am interested in the genera associated with Actinobacteriota in each treatment category
actino_abundances <- melt_pitBiochar %>%
  filter(Phylum == "Actinobacteriota")
unique(actino_abundances$Genus) # there are 57 unique genera,
unique(actino_abundances$OTU) #  there are 477 unique ASVs
# number of unique ASVs by sample 
temp <- actino_abundances %>%
  group_by(Sample,Genus) %>%
  summarise(totalAbund = sum(Abundance),
            n = n())
temp$bulk_soil <- substr(temp$Sample, 1,1)
output <- temp %>%
  group_by(bulk_soil, Genus) %>%
  summarise(meantotalAbund = mean(totalAbund))
output <- output %>%
  filter(meantotalAbund > 0)
output %>%
  ggplot()+
  geom_point(mapping = aes( x=Genus, y = meantotalAbund, color = bulk_soil))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))

## Now I am interested in the abundances of Verrucomicrobiota because they were unique to biochar # But there were not many genera identified. The most dominant (>75%) were from Chthoniobacter. A smaller amount were from Candidatus Udaeobacter. 

# Evenness with estimateEvenness(index = "simpson")
# Dominance - usually negatively correlated with diversity - estimateDominance(index = "relative") which gives the relative abundance of the most dominant species in the sample 














# Mantel Test for soil characteristics and whole community dissimilarity measurements
# First, I only measured these variables for controls and pits, so make a new subset 
samplesToKeep2 <- c("P22a","P22b","P22c","C22a","C22b","C22c")
pitBiochar2 <- prune_samples(sample_names(pits_16S_filt_ps) %in% samplesToKeep2 , pits_16S_filt_ps)
# factor the bulk soil levels
pitBiochar2@sam_data$bulk_soil <- factor(pitBiochar2@sam_data$bulk_soil, levels = c("C","P"))
# make new bray curtis dissimilarity matrix
bray_dist2 <- phyloseq::distance(pitBiochar2, method="bray", weighted= F)

# Mantel tests
ph <- pitBiochar2@sam_data$pH
dist.ph <- dist(ph, method = "euclidean")
abund_ph <- mantel(bray_dist2, dist.ph, method = "spearman", permutations = 999, na.rm = TRUE)
abund_ph #Mantel statistic r: 0.3071, Significance: 0.077778 

bd <- pitBiochar2@sam_data$bulk_density_g_cm3
dist.bd <- dist(bd, method = "euclidean")
abund_bd <- mantel(bray_dist2, dist.bd, method = "spearman", permutations = 999, na.rm = TRUE)
abund_bd #Mantel statistic r: 0.1441, Significance: 0.23333 


sand <- pitBiochar2@sam_data$sand_percent
dist.sand <- dist(sand, method = "euclidean")
abund_sand <- mantel(bray_dist2, dist.sand, method = "spearman", permutations = 999, na.rm = TRUE)
abund_sand #Mantel statistic r: 0.8393, Significance: 0.016667 


silt <- pitBiochar2@sam_data$silt_percent
dist.silt <- dist(silt, method = "euclidean")
abund_silt <- mantel(bray_dist2, dist.silt, method = "spearman", permutations = 999, na.rm = TRUE)
abund_silt #Mantel statistic r: 0.6893, Significance: 0.0625 

clay <- pitBiochar2@sam_data$clay_percent
dist.clay <- dist(clay, method = "euclidean")
abund_clay <- mantel(bray_dist2, dist.clay, method = "spearman", permutations = 999, na.rm = TRUE)
abund_clay #Mantel statistic r: 0.3929, Significance: 0.13472

water <- pitBiochar2@sam_data$MarchApril_meanVWC
dist.water <- dist(water, method = "euclidean")
abund_water <- mantel(bray_dist2, dist.water, method = "spearman", permutations = 999, na.rm = TRUE)
abund_water #Mantel statistic r: 0.8504, Significance: 0.1






#### Tracking ASVs ####
pits_psmelt <- left_join(pits_psmelt, reads_per_sample, by = "Sample")
pits_psmelt$RelAbund <- pits_psmelt$Abundance / pits_psmelt$totalReads
pits_ASVs <- pits_psmelt %>% group_by(OTU) %>% summarise(sumrelAbund = sum(RelAbund))
pits_ASVs <- pits_ASVs %>% arrange(desc(sumrelAbund))
sum(pits_ASVs$sumrelAbund)
# 9 because each sample has 1 relative abundance and we have 9 samples 
pits_top10ASVs <- pits_ASVs$OTU[1:10]
pits_top10ASVs
# The top 20 ASVs represent ~30% of the reads in the 9 samples 
# The top 10 ASVs represent ~20% of the reads in the 9 samples
pits_psmelt_top10ASV <- pits_psmelt %>% filter(OTU %in% pits_top10ASVs)

# Add taxonomy information
pits_psmelt_top10ASV$taxa <- paste(pits_psmelt_top10ASV$Phylum, ";", pits_psmelt_top10ASV$Class, ";", pits_psmelt_top10ASV$Order, ";", pits_psmelt_top10ASV$Family, ";", pits_psmelt_top10ASV$Genus, ";", pits_psmelt_top10ASV$OTU)


# Keep individual samples (replicates)
pits_psmelt_top10ASV$bulk_soil <- factor(pits_psmelt_top10ASV$bulk_soil, levels = c("C", "P", "BP"))

# pits_psmelt_top10ASV %>% 
#   filter(OTU == "ASV_3") %>%
#   ggplot()+
#   geom_point(mapping = aes(x = bulk_soil, y = RelAbund))+
#   geom_line(mapping = aes(x = bulk_soil, y = RelAbund, group = Sample))+
#   theme_classic()+
#   facet_grid(cols = vars(OTU))




# summarise for plotting
pits_psmelt_top10ASV_summary <- pits_psmelt_top10ASV %>% 
  group_by(bulk_soil, OTU, taxa) %>%
  summarise(meanRelAbund = mean(RelAbund),
            sdRelAbund = sd(RelAbund))

# factor for plotting
pits_psmelt_top10ASV_summary$bulk_soil <- factor(pits_psmelt_top10ASV_summary$bulk_soil, levels = c("C", "P", "BP"))

pits_psmelt_top10ASV_summary %>% 
  ggplot()+
  geom_point(mapping = aes(x = bulk_soil, y = meanRelAbund, group = OTU, color = taxa))+
  #geom_errorbar(mapping = aes(x = bulk_soil, ymin = meanRelAbund-sdRelAbund, ymax = meanRelAbund+sdRelAbund, group = OTU, color = taxa))+
  geom_line(mapping = aes(x = bulk_soil, y = meanRelAbund, group = OTU, color = taxa))+
  theme_classic()



#### Genera
pits_Genera <- pits_psmelt %>% group_by(Genus) %>% summarise(sumrelAbund = sum(RelAbund))
pits_Genera <- pits_Genera %>% arrange(desc(sumrelAbund))
pits_Genera <- pits_Genera %>% filter(!Genus %in% c("NA"))
sum(pits_Genera$sumrelAbund)

# 5.8 because I deleted a lot of NA Genera which would make this sum to 9
pits_top50Genera <- pits_Genera$Genus[1:50]
pits_top50Genera

pits_psmelt_top50Genera <- pits_psmelt %>% filter(Genus %in% pits_top50Genera)

# Add taxonomy information
pits_psmelt_top50Genera$taxa <- paste(pits_psmelt_top50Genera$Phylum, ";", pits_psmelt_top50Genera$Class, ";", pits_psmelt_top50Genera$Order, ";", pits_psmelt_top50Genera$Family, ";", pits_psmelt_top50Genera$Genus)

pits_psmelt_top50Genera$taxa <- gsub(" ;", ";", pits_psmelt_top50Genera$taxa)

# Keep individual samples (replicates)
pits_psmelt_top50Genera$bulk_soil <- factor(pits_psmelt_top50Genera$bulk_soil, levels = c("C", "P", "BP"))

# summarise for plotting
pits_psmelt_top50Genera_summary <- pits_psmelt_top50Genera %>% 
  group_by(bulk_soil, Genus, taxa, Phylum) %>%
  summarise(meanRelAbund = mean(RelAbund),
            sdRelAbund = sd(RelAbund))

# factor for plotting
pits_psmelt_top50Genera_summary$bulk_soil <- factor(pits_psmelt_top50Genera_summary$bulk_soil, levels = c("C", "P", "BP"))

# new column to highlight the correct cyanobacteria name 
pits_psmelt_top50Genera_summary$CydrasilName <- pits_psmelt_top50Genera_summary$taxa
pits_psmelt_top50Genera_summary$CydrasilName[pits_psmelt_top50Genera_summary$CydrasilName == "Cyanobacteria; Cyanobacteriia; Cyanobacteriales; Phormidiaceae; Tychonema CCAP 1459-11B"] <- "Microcoleus vaginatus"

unique(pits_psmelt_top50Genera_summary$CydrasilName)

# basic version 
pits_psmelt_top50Genera_summary %>% 
  ggplot()+
  geom_point(mapping = aes(x = bulk_soil, y = meanRelAbund, group = Genus, color = taxa))+
  #geom_errorbar(mapping = aes(x = bulk_soil, ymin = meanRelAbund-sdRelAbund, ymax = meanRelAbund+sdRelAbund, group = OTU, color = taxa))+
  geom_line(mapping = aes(x = bulk_soil, y = meanRelAbund, group = Genus, color = taxa))+
  theme_classic()+
  theme(legend.position = "bottom")+
  guides(color=guide_legend(ncol=2))

#ggsave("output/ASV_microsites.pdf", width=15, height=15)

# A cleaned up version for submitting for publication (08/30/2025)
pits_psmelt_top50Genera_summary %>% 
  ggplot()+
  geom_point(mapping = aes(x = bulk_soil, y = meanRelAbund, group = Genus, color = CydrasilName))+
  #geom_errorbar(mapping = aes(x = bulk_soil, ymin = meanRelAbund-sdRelAbund, ymax = meanRelAbund+sdRelAbund, group = OTU, color = taxa))+
  geom_line(mapping = aes(x = bulk_soil, y = meanRelAbund, group = Genus, color = CydrasilName))+
  scale_x_discrete(labels = c("Control", "Pit", "Pit+Biochar"))+
  scale_color_manual(values = c("Microcoleus vaginatus" = "blue", "other" = "black")) + # Set desired colors
  theme_classic()+
  labs(x = "", y = "Mean Relative Abundance", color = "Genus")+
  theme(legend.position = "bottom")+
  guides(color=guide_legend(ncol=2))

ggsave("output/ASV_microsites_clean.pdf", width=6, height=5)


# subset this for cyanobacteria
pits_psmelt_top50Genera_cyano <- pits_psmelt_top50Genera %>% filter(Phylum == "Cyanobacteria")
unique(pits_psmelt_top50Genera_cyano$OTU)
# there are three cyanobacteria ASVs that made the top 50 taxa and they are all M. vaginatus

subset(pits_psmelt, OTU == "ASV_3")[1,] # M. vaginatus
subset(pits_psmelt, OTU == "ASV_121")[1,] # M. vaginatus
subset(pits_psmelt, OTU == "ASV_325")[1,] # M. vaginatus
