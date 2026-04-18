# The purpose of this code is to think about whether or not the soil microbiome at the degraded site differs from reference soils. Part of this comparison includes thinking about pitting as digging into the soil profile and exposing deep microbes to surface conditions, which may account for some of the differences we see between controls and pits. 

### Load libraries
library(tidyverse)
library(phyloseq) # For all bioinformatics analyses post Dada2
packageVersion('phyloseq') # ‘1.52.0’
library(ggplot2)
library(vegan)
library(car) # For Levene Test and anova
library(PMCMRplus) # For Nemenyi posthoc test
#library("DESeq2")
library(indicspecies) # For indicator species
library(microbiome)

# Color palettes
grafify_muted_greenstart <- c("#117733", "#88ccee", "#882255", "#44aa99", "#999933", "#aa4499", "#dddddd","#cc6677", "#332288", "#ddcc77")

natural_colors <- c("brown4", "mediumpurple4","bisque3","cadetblue3","darkslategray","goldenrod","lightsteelblue4","coral2", "cornflowerblue", "cornsilk")

natural_colors2 <- c("cornsilk", "cornflowerblue", "coral2", "lightsteelblue4", "goldenrod", "darkslategray", "cadetblue3", "bisque3", "mediumpurple4", "brown4")

### Load the phyloseq data table which comes from the pits16Sphyloseq_1.10.22.R script (note the date should be 1.10.23 - I got the year wrong). 
list.files()
pits_16S_rar_filt_pruned_ps <- readRDS("pits_processed_filt_ps.rds")
#pits_16S_rar_filt_pruned_ps@sam_data

qPCR_masses <- read.csv("qPCR_mass_extracted.csv", header = TRUE)
qPCR <- as.data.frame(pits_16S_rar_filt_pruned_ps@sam_data$qPCR_16Scopies)
qPCR$treatment <- pits_16S_rar_filt_pruned_ps@sam_data$bulk_soil
qPCR$Sample_ID <- sample_names(pits_16S_rar_filt_pruned_ps)
colnames(qPCR)[1] <- c("qPCR_count")
qPCR <- left_join(qPCR, qPCR_masses, by = c("Sample_ID" = "SampleID"))
# remove rows with NA
qPCR <- drop_na(qPCR)

# Let's compare the microbial communities from 2022 in healthy and control plots
bulk <- subset_samples(pits_16S_rar_filt_pruned_ps, rhizosphere_bulk=="b" & year == 2022)
bulk@sam_data
#HC <- subset_samples(bulk, bulk_soil != "BP")
#HC@sam_data
# removing deep samples too
HC <- subset_samples(bulk, !SampleID_norep %in% c("BP22","CD22","HbD22", "P22"))
HC@sam_data



#### Plotting Phyla Relative Abundance ####

#Get count of reads by phyla
table(phyloseq::tax_table(HC)[, "Phylum"]) #29 Phyla included
# get the total read count per sample
reads_per_sample <- as.data.frame (sample_sums(HC)) # they are all different because we did not rarefy in this dataset
#copy the rownames to a column
reads_per_sample <- tibble::rownames_to_column(reads_per_sample, "Sample")
# change the name of the second column
colnames(reads_per_sample) <- c("Sample", "totalReads")
nrow(reads_per_sample) # 6 samples are included

# Think about the phyla 
#get_taxa_unique(initial_biocrust, "Phylum") # 34 unique Phyla
# I only want to plot the top 9 and the rest go into an Other category

# first melt the data out of phyloseq format and into long dataframe
pits_psmelt <- psmelt(HC)

# Get a list of the most abundant Phyla by comboLabel
phyla_summary <- pits_psmelt %>%
  group_by(bulk_soil, Phylum) %>%
  summarise(totalAbund = sum(Abundance)) %>%
  arrange(-totalAbund)

# Get a list of the most abundant Phyla
temp <- pits_psmelt %>%
  group_by(Phylum) %>%
  summarise(totalAbund = sum(Abundance)) %>%
  arrange(-totalAbund)

temp <- temp[1:9,]

top9Phyla <- temp$Phylum # list that holds the Top9Phyla names
top9Phyla

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
  dplyr::select(Sample, other)
#add column which is "Other" repeated
total_relAbund$Phylum <- "Other"
#change column header "other" to relAbund
colnames(total_relAbund)[which(names(total_relAbund) == "other")] <- "relAbund"
# select columns to keep in the dataframe we want
phylaSubset_abund2 <- phylaSubset_abund2 %>%
  dplyr::select(Sample, relAbund, Phylum)
# rbind the other values to the phylaSubset_abund2 dataframe
phylaSubset_abund3 <- rbind(phylaSubset_abund2, total_relAbund)

# Now check that they sum to 1
total_relAbund2 <- phylaSubset_abund3 %>%
  group_by(Sample) %>%
  summarise(sum = sum(relAbund)) # YES!!!!
head(total_relAbund2)

# plot it!
# Plot Relative Abundance
temp2 <- unique(pits_psmelt %>% dplyr::select(Sample, bulk_soil))
phylaSubset_abund3 <- left_join(phylaSubset_abund3, temp2, by = "Sample")

phylaSubset_abund3$Phylum <- factor(phylaSubset_abund3$Phylum, levels = c("Other", "Armatimonadota", "Planctomycetota", "Cyanobacteria", "Acidobacteriota", "Bacteroidota", "Crenarchaeota", "Proteobacteria", "Chloroflexi", "Actinobacteriota"))

# Plot Relative Abundance - basic
phylaSubset_abund3 %>%
  ggplot(aes(y = relAbund, x = Sample, fill = Phylum)) +
  geom_bar(stat = "identity", alpha = 0.8, width = 0.85) +
  scale_fill_manual(values = grafify_muted_greenstart, breaks=c("Actinobacteriota","Chloroflexi",  "Proteobacteria","Crenarchaeota", "Bacteroidota","Acidobacteriota", "Cyanobacteria", "Planctomycetota","Armatimonadota", "Other")) +
  scale_y_continuous(expand = c(0, 0.01), breaks=seq(0, 1, .2)) +
  theme_classic()+
  #coord_flip()+
  #facet_grid(rows = vars(bulk_soil), scales = "free", space = "free", switch = "y", labeller = labeller(bulk_soil = new_labels))+
  labs(x = "", y = "Relative Abundance", fill = "")+
  theme(axis.text.x = element_text(size = 12),
        #axis.text.y = element_blank(),
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
# note Hb22a has significantly more cyanobacteria reads than the other samples and I should consider removing it for statistical analyses? I am going to keep it in because I only have three samples of each and because it is real. When I remove this sample from the dataset, the results do not change (statistical comparisons at the phylum level)



# All together but average (not every replicate)
# first calculate the averages
phylaSubset_abund3_agregate <- phylaSubset_abund3 %>% 
  #filter(Sample != "Hb22a") %>%
  group_by(bulk_soil, Phylum) %>%
  summarise(meanRel_abund = mean(relAbund),
            sdRel_abund = sd(relAbund))

#export this dataframe for use in other scripts
#write.csv(phylaSubset_abund3_agregate, "output_v1/phylaSubsetagg_initialbiocrust.csv", row.names = FALSE)


# plot: Relative abundance, phylum level
new_labels <- c("H" = "Reference", "C" = "Degraded")

phylaSubset_abund3_agregate$bulk_soil <- factor(phylaSubset_abund3_agregate$bulk_soil, levels = c("H", "C"))

phylaSubset_abund3_agregate$Phylum <- factor(phylaSubset_abund3_agregate$Phylum, levels = c("Other", "Armatimonadota", "Planctomycetota", "Cyanobacteria", "Acidobacteriota", "Bacteroidota", "Crenarchaeota", "Proteobacteria", "Chloroflexi", "Actinobacteriota"))

phylaSubset_abund3_agregate %>%
  ggplot(aes(y = meanRel_abund, x = bulk_soil, fill = Phylum)) +
  geom_bar(stat = "identity", width = 0.85, color = "black", size = 0.25) +
  scale_fill_manual(values = natural_colors, breaks=c("Actinobacteriota","Chloroflexi",  "Proteobacteria","Crenarchaeota", "Bacteroidota","Acidobacteriota", "Cyanobacteria", "Planctomycetota","Armatimonadota", "Other")) +
  scale_y_continuous(expand = c(0, 0.01), breaks=seq(0, 1, .2)) +
  theme_classic()+
  coord_flip()+
  facet_grid(rows = vars(bulk_soil), scales = "free", space = "free", switch = "y", labeller = labeller(bulk_soil = new_labels))+
  labs(x = "", y = "Relative Abundance", fill = "")+
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

#ggsave("output/HD_relabund_stackedbar_naturalcolors.pdf", width = 6.5, height = 3.5)


# Same graph but vertical to match the other one in the paper (08/30/2025)
phylaSubset_abund3_agregate$Phylum <- factor(phylaSubset_abund3_agregate$Phylum, levels = c("Actinobacteriota","Chloroflexi",  "Proteobacteria","Crenarchaeota", "Bacteroidota","Acidobacteriota", "Cyanobacteria", "Planctomycetota","Armatimonadota", "Other"))

phylaSubset_abund3_agregate %>%
  ggplot(aes(y = meanRel_abund, x = bulk_soil, fill = Phylum)) +
  geom_bar(stat = "identity", color = "black", size = 0.25) +
  scale_fill_manual(values = natural_colors, breaks=c("Actinobacteriota","Chloroflexi",  "Proteobacteria","Crenarchaeota", "Bacteroidota","Acidobacteriota", "Cyanobacteria", "Planctomycetota","Armatimonadota", "Other")) +
  scale_x_discrete(labels = c("Reference","Degraded"))+
  scale_y_continuous(expand = c(0, 0.01), breaks=seq(0, 1, 0.25)) +
  theme_classic()+
  #coord_flip()+
  #facet_grid(rows = vars(bulk_soil), scales = "free", space = "free", switch = "y", labeller = labeller(bulk_soil = new_labels))+
  labs(x = "", y = "Relative Abundance", fill = "")+
  theme(axis.title = element_text(face = "bold", size = 15),
              axis.text.x = element_text(size = 14),
              axis.text.y = element_text(size = 14),
              legend.text = element_text(face = "italic", size = 14),
              axis.ticks.x = element_blank(),
              legend.title=element_text(face = "bold",size=14))
  guides(fill = guide_legend(nrow = 3, byrow = TRUE))

ggsave("output/HD_relabund_stackedbar_naturalcolors.pdf", width = 6, height = 5)





# # For supplementary figure, take a different subset of the data
# temp_12phyla <- c("Actinobacteriota", "Bacteroidota", "Acidobacteriota", "Chloroflexi", "Proteobacteria", "Thaumarchaeota", "Planctomycetota", "Cyanobacteria", "Crenarchaeota", "Verrucomicrobia", "Armatimonadota")
# 
# #subset the psmelt dataframe for only the Top12Phyla
# phyla12Subset <- pits_psmelt %>%
#   filter(Phylum %in% temp_12phyla)
# 
# # check that we got the right ones
# unique(phyla12Subset$Phylum) # yes
# 
# # calculate the reads per Phylum, grouped by Sample
# phyla12Subset_abund <- phyla12Subset %>%
#   group_by(Sample, Phylum) %>%
#   summarise(abund = sum(Abundance))
# 
# #make a new column were I divide each of the values by the correct one in the reads_per_sample list
# # first merge the reads_per_sample data onto the dataframe
# phyla12Subset_abund2 <- inner_join(phyla12Subset_abund, reads_per_sample, by = "Sample")
# # divide the abund column by the total reads column
# phyla12Subset_abund2$relAbund <- phyla12Subset_abund2$abund / phyla12Subset_abund2$totalReads
# # check that it does not add to 1 because there should be some in the Other category
# total_12relAbund <- phyla12Subset_abund2 %>%
#   group_by(Sample) %>%
#   summarise(sum = sum(relAbund))
# # add a new column which contains the relative abundance value for the Other category
# total_12relAbund$other <- (1 - total_12relAbund$sum)
# # delete the sum column
# total_12relAbund <- total_12relAbund %>%
#   dplyr::select(Sample, other)
# #add column which is "Other" repeated
# total_12relAbund$Phylum <- "Other"
# #change column header "other" to relAbund
# colnames(total_12relAbund)[which(names(total_12relAbund) == "other")] <- "relAbund"
# # select columns to keep in the dataframe we want
# phyla12Subset_abund2 <- phyla12Subset_abund2 %>%
#   dplyr::select(Sample, relAbund, Phylum)
# # rbind the other values to the phylaSubset_abund2 dataframe
# phyla12Subset_abund3 <- rbind(phyla12Subset_abund2, total_12relAbund)
# 
# # Now check that they sum to 1
# total_12relAbund2 <- phyla12Subset_abund3 %>%
#   group_by(Sample) %>%
#   summarise(sum = sum(relAbund)) # YES!!!!
# head(total_12relAbund2)
# 
# temp_12phyla <- unique(pits_psmelt %>% dplyr::select(Sample, bulk_soil))
# phyla12Subset_abund3 <- left_join(phyla12Subset_abund3, temp_12phyla, by = "Sample")
# 
# # All together but average (not every replicate)
# # first calculate the averages
# phyla12Subset_abund3_aggregate <- phyla12Subset_abund3 %>% 
#   #filter(Sample != "Hb22a") %>%
#   group_by(bulk_soil, Phylum) %>%
#   summarise(meanRel_abund = mean(relAbund),
#             sdRel_abund = sd(relAbund))
# 
# phyla12Subset_abund3_aggregate$seqRound <- 2
# 
# #export this dataframe for use in other scripts
# #write.csv(phyla12Subset_abund3_aggregate, "output/suppFig_12phyla_round2_aggregate.csv", row.names = FALSE)



# Check out the firmicutes in the bulk soils
firmicutes <- pits_psmelt %>% filter(Phylum == "Firmicutes")
unique(firmicutes$Genus) # 25 genera present
firmicutes_summary <- firmicutes %>% group_by(Genus, bulk_soil) %>% summarise(sumAbund = sum(Abundance))






# Where did the Thaumarchaeota and Verrucomicrobia go? 
#unique(pits_psmelt$Phylum)
# there aren't any in the dataset...not sure why

# Statistical comparison at phylum level 
# Calculate median and standard deviation for a table in the manuscript
phylaSubset_medians <- phylaSubset_abund3 %>% 
  #filter(Sample != "Hb22a") %>%
  group_by(Phylum, bulk_soil) %>%
  summarise(medianRelAbund = round(median(relAbund), 3),
            sdRelAbund = round(sd(relAbund),2))

# Actinobacteria
actino <- phylaSubset_abund3 %>% filter(Phylum == "Actinobacteriota")
t.test(relAbund ~ bulk_soil, alternative = "two.sided", data = actino)
# t = 1.0691, df = 2.0537, p-value = 0.3945
# check assumptions
with(actino, shapiro.test(relAbund[bulk_soil == "H"]))# p = 0.1703 (normal)
with(actino, shapiro.test(relAbund[bulk_soil == "C"]))# p = 0.0892 (normal)# p = 0.3295 (normal)
var.test(relAbund ~ bulk_soil, data = actino) # p = 0.0265, variances are statistically different
# cannot use t-test, use wilcoxon instead
wilcox.test(relAbund ~ bulk_soil, data = actino)
# W = 7, p-value = 0.4
phylaSubset_abund3 %>% 
  #filter(Sample != "Hb22a") %>%
  filter(Phylum == "Actinobacteriota") %>%
  ggplot() + 
  geom_boxplot(mapping = aes(x = bulk_soil, y = relAbund))+
  geom_jitter(mapping = aes(x = bulk_soil, y = relAbund))+
  theme_classic() # note the low value for one healthy sample



# Acidobacteriota
acido <- phylaSubset_abund3 %>% filter(Phylum == "Acidobacteriota")
t.test(relAbund ~ bulk_soil, alternative = "two.sided", data = acido)
# t = 1.1772, df = 2.2073, p-value = 0.3504
# check assumptions
with(acido, shapiro.test(relAbund[bulk_soil == "H"]))# p = 0.1625 (normal)
with(acido, shapiro.test(relAbund[bulk_soil == "C"]))# p = 0.3556 (normal)# p = 0.3295 (normal)
var.test(relAbund ~ bulk_soil, data = acido) # p = 0.01524, variances are not statistically different
# t.test ok to use
phylaSubset_abund3 %>% 
  filter(Phylum == "Acidobacteriota") %>%
  ggplot() + 
  geom_boxplot(mapping = aes(x = bulk_soil, y = relAbund))+
  geom_jitter(mapping = aes(x = bulk_soil, y = relAbund))+
  theme_classic()

# Chloroflexi
chloro <- phylaSubset_abund3 %>% filter(Phylum == "Chloroflexi")
t.test(relAbund ~ bulk_soil, alternative = "two.sided", data = chloro)
# t = 0.069211, df = 2.1712, p-value = 0.9507
# check assumptions
with(chloro, shapiro.test(relAbund[bulk_soil == "H"]))# p = 0.572 (normal)
with(chloro, shapiro.test(relAbund[bulk_soil == "C"]))# p = 0.1344 (normal)
var.test(relAbund ~ bulk_soil, data = chloro) # p = 0.08225, variances are not statistically different
# t.test ok to use
phylaSubset_abund3 %>% 
  filter(Phylum == "Chloroflexi") %>%
  ggplot() + 
  geom_boxplot(mapping = aes(x = bulk_soil, y = relAbund))+
  geom_jitter(mapping = aes(x = bulk_soil, y = relAbund))+
  theme_classic()


# Proteobacteria
protob <- phylaSubset_abund3 %>% filter(Phylum == "Proteobacteria")
t.test(relAbund ~ bulk_soil, alternative = "two.sided", data = protob)
# t = 0.091998, df = 2.4342, p-value = 0.9337
# check assumptions
with(protob, shapiro.test(relAbund[bulk_soil == "H"]))# p = 0.6015 (normal)
with(protob, shapiro.test(relAbund[bulk_soil == "C"]))# p = 0.9777 (normal)
var.test(relAbund ~ bulk_soil, data = protob) # p = 0.198, variances are not statistically different
# t.test ok to use
phylaSubset_abund3 %>% 
  filter(Phylum == "Proteobacteria") %>%
  ggplot() + 
  geom_boxplot(mapping = aes(x = bulk_soil, y = relAbund))+
  geom_jitter(mapping = aes(x = bulk_soil, y = relAbund))+
  theme_classic()


# Crenarchaeota
cren <- phylaSubset_abund3 %>% filter(Phylum == "Crenarchaeota")
t.test(relAbund ~ bulk_soil, alternative = "two.sided", data = cren)
# t = 2.9922, df = 3.3123, p-value = 0.05117
# check assumptions
with(cren, shapiro.test(relAbund[bulk_soil == "H"]))# p = 0.4156 (normal)
with(cren, shapiro.test(relAbund[bulk_soil == "C"]))# p = 0.1829 (normal)
var.test(relAbund ~ bulk_soil, data = cren) # p = 0.5443, variances are not statistically different
# t.test ok to use
phylaSubset_abund3 %>% 
  filter(Phylum == "Crenarchaeota") %>%
  ggplot() + 
  geom_boxplot(mapping = aes(x = bulk_soil, y = relAbund))+
  geom_jitter(mapping = aes(x = bulk_soil, y = relAbund))+
  theme_classic()

# Bacteroidota
bacter <- phylaSubset_abund3 %>% filter(Phylum == "Bacteroidota")
t.test(relAbund ~ bulk_soil, alternative = "two.sided", data = bacter)
# t = -1.5246, df = 3.7525, p-value = 0.2066
# check assumptions
with(bacter, shapiro.test(relAbund[bulk_soil == "H"]))# p = 0.6319 (normal)
with(bacter, shapiro.test(relAbund[bulk_soil == "C"]))# p = 0.01282 (normal)
var.test(relAbund ~ bulk_soil, data = bacter) # p = 0.7432, variances are not statistically different
# t.test not ok to use, wilcoxon test
wilcox.test(relAbund ~ bulk_soil, data = bacter)
#W = 1, p-value = 0.2
phylaSubset_abund3 %>% 
  filter(Phylum == "Bacteroidota") %>%
  ggplot() + 
  geom_boxplot(mapping = aes(x = bulk_soil, y = relAbund))+
  geom_jitter(mapping = aes(x = bulk_soil, y = relAbund))+
  theme_classic()


# Cyanobacteria
cyano <- phylaSubset_abund3 %>% filter(Phylum == "Cyanobacteria")
t.test(relAbund ~ bulk_soil, alternative = "two.sided", data = cyano)
# t = -1.2532, df = 2.0399, p-value = 0.3347
# check assumptions
with(cyano, shapiro.test(relAbund[bulk_soil == "H"]))# p = 0.2604 (normal)
with(cyano, shapiro.test(relAbund[bulk_soil == "C"]))# p = 0.3667 (normal)
var.test(relAbund ~ bulk_soil, data = cyano) # p = 0.01976, variances are not statistically different
# t.test not ok to use, wilcoxon
wilcox.test(relAbund ~ bulk_soil, data = cyano)
# W = 1, p-value = 0.2
phylaSubset_abund3 %>% 
  filter(Phylum == "Cyanobacteria") %>%
  ggplot() + 
  geom_boxplot(mapping = aes(x = bulk_soil, y = relAbund))+
  geom_jitter(mapping = aes(x = bulk_soil, y = relAbund))+
  theme_classic()


# Planctomycetota
plancto <- phylaSubset_abund3 %>% filter(Phylum == "Planctomycetota")
t.test(relAbund ~ bulk_soil, alternative = "two.sided", data = plancto)
# t = 1.5756, df = 3.925, p-value = 0.1916
# check assumptions
with(plancto, shapiro.test(relAbund[bulk_soil == "H"]))# p = 0.4619 (normal)
with(plancto, shapiro.test(relAbund[bulk_soil == "C"]))# p = 0.5702 (normal)
var.test(relAbund ~ bulk_soil, data = plancto) # p = 0.8617, variances are not statistically different
# t.test ok to use
phylaSubset_abund3 %>% 
  filter(Phylum == "Planctomycetota") %>%
  ggplot() + 
  geom_boxplot(mapping = aes(x = bulk_soil, y = relAbund))+
  geom_jitter(mapping = aes(x = bulk_soil, y = relAbund))+
  theme_classic()

# Armatimonadota
arma <- phylaSubset_abund3 %>% filter(Phylum == "Armatimonadota")
t.test(relAbund ~ bulk_soil, alternative = "two.sided", data = arma)
# t = 0.11206, df = 2.081, p-value = 0.9207
# check assumptions
with(arma, shapiro.test(relAbund[bulk_soil == "H"]))# p = 0.4985 (normal)
with(arma, shapiro.test(relAbund[bulk_soil == "C"]))# p = 0.0737 (normal)
var.test(relAbund ~ bulk_soil, data = arma) # p = 0.03973, variances are not statistically different
# t.test not ok to use, wilcoxon test
wilcox.test(relAbund ~ bulk_soil, data = arma)
# W = 3, p-value = 0.7
phylaSubset_abund3 %>% 
  filter(Phylum == "Armatimonadota") %>%
  ggplot() + 
  geom_boxplot(mapping = aes(x = bulk_soil, y = relAbund))+
  geom_jitter(mapping = aes(x = bulk_soil, y = relAbund))+
  theme_classic()

### Statistical Comparison with Ranks 
abund <- otu_table(HC)
# First remove rows that sum to zero
abund <- abund[rowSums(abund[]) > 0, ]
abund_ranks <- apply(abund, 2, rank, ties.method = "average") # calculation of the ranks at ASV level with the suggested method for calculating ties 
# Thinking about log transformation as a normalization (proof that you do not need to log transform when using ranks)
#pslog <- transform_sample_counts(HC, function(x) log(1+x)) # log transformation to deal with normality
#abund_1 <- otu_table(pslog)
#abund_logranks <- apply(abund_1, 2, rank, ties.method = "average")

#  In the case where many bacteria are absent or present at trace amounts, an artificially large difference in rank could occur for minimally abundant taxa. To avoid this, all those microbes with rank below some threshold are set to be tied at 1. The ranks for the other microbes are shifted down, so there is no large gap between ranks.
# I am not sure how to calculate an appropriate threshold
hist(abund_ranks[,1])
# max is > 1000
# min is 311 (the value given to the ties average)
abund_ranks <- abund_ranks - 311
abund_ranks[abund_ranks < 1] <- 1

# So, now run the t-tests on the ranked data 
# First, add the rank data to the dataframe ... these are at the ASV level so I would need to recalculate at the right level for the statistical test 
abund_ranks <- as.data.frame(abund_ranks)
abund_ranks <- rownames_to_column(abund_ranks)
colnames(abund_ranks)[1] <- "ASV"
datalabels <- pits_psmelt %>% select("OTU", "SampleID_norep","Kingdom":"Genus")
abund_ranks <- left_join(abund_ranks, datalabels, by = c("ASV"="OTU"))

# IDK WHAT COMES NEXT 






## Absolute Abundances
# Reporting as No of 16S rRNA copies and applying across all phyla
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
phylaSubset_abund3$abso_abund <- phylaSubset_abund3$absolute_abundance / (10^7)


# Not aggregated 
#phylaSubset_abund3$Phylum <- factor(phylaSubset_abund3$Phylum, levels = c("Other", "Armatimonadota", "Planctomycetota", "Cyanobacteria", "Acidobacteriota", "Bacteroidota", "Crenarchaeota", "Proteobacteria", "Chloroflexi", "Actinobacteriota"))

# Plot Relative Abundance - basic
phylaSubset_abund3 %>%
  ggplot(aes(y = abso_abund, x = Sample, fill = Phylum)) +
  geom_bar(stat = "identity", alpha = 0.8, width = 0.85) +
  scale_fill_manual(values = grafify_muted_greenstart, breaks=c("Actinobacteriota","Chloroflexi",  "Proteobacteria","Crenarchaeota", "Bacteroidota","Acidobacteriota", "Cyanobacteria", "Planctomycetota","Armatimonadota", "Other")) +
  scale_y_continuous(expand = c(0, 0.01), breaks=seq(0, 10, 1)) +
  theme_classic()+
  coord_flip()+
  #facet_grid(rows = vars(bulk_soil), scales = "free", space = "free", switch = "y", labeller = labeller(bulk_soil = new_labels))+
  labs(x = "", y = expression(paste("No. 16S rRNA gene copies x ", 10^7, " x "  ,g^-1,)), fill = "")+
  theme(axis.text.x = element_text(size = 12),
        #axis.text.y = element_blank(),
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
# note Hb22a has significantly more reads than the other samples and I should consider removing it for statistical analyses... I am going to keep it in because I only have three samples of each and because it is real. When I remove this sample from the dataset, the results do not change (statistical comparisons at the phylum level)???



# All together but average (not every replicate)
# first calculate the averages
phylaSubset_abund3_aggregate <- phylaSubset_abund3 %>% 
  group_by(bulk_soil, Phylum) %>%
  summarise(meanAbso_abund = mean(abso_abund),
            sdAbso_abund = sd(abso_abund))

# plot: Absolute abundance phylum level for supplementary materials
phylaSubset_abund3_aggregate$Phylum <- factor(phylaSubset_abund3_aggregate$Phylum, levels = c("Other", "Armatimonadota", "Planctomycetota", "Cyanobacteria", "Acidobacteriota", "Bacteroidota", "Crenarchaeota", "Proteobacteria", "Chloroflexi", "Actinobacteriota"))

phylaSubset_abund3_aggregate %>%
  ggplot(aes(y = meanAbso_abund, x = bulk_soil, fill = Phylum)) +
  geom_bar(stat = "identity", alpha = 0.8, width = 0.85) +
  scale_fill_manual(values = grafify_muted_greenstart, breaks=c("Actinobacteriota","Chloroflexi",  "Proteobacteria","Crenarchaeota", "Bacteroidota","Acidobacteriota", "Cyanobacteria", "Planctomycetota","Armatimonadota", "Other")) +
  scale_y_continuous(expand = c(0, 0.1), breaks=seq(0, 6, 1)) +
  theme_classic()+
  coord_flip()+
  facet_grid(rows = vars(bulk_soil), scales = "free", space = "free", switch = "y",labeller = labeller(bulk_soil = new_labels))+
  labs(x = "", y = expression(paste("No. 16S rRNA gene copies x ", 10^7, " x "  ,g^-1,)), fill = "")+
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

#ggsave("output/HD_phyla_abso_abund_stackedbar.pdf", width = 6.5, height = 3.5)

# I do not know that abosolute abundances are any more helpful in this story than relative abundances, so for now, I am going to not use them at all. 










# Is richness different in these bulk soils?
HC.rarefied = rarefy_even_depth(HC, rngseed=123, sample.size=0.9*min(sample_sums(HC)), replace=F)
phyloseq::sample_sums(HC.rarefied) # all rarefied to 20455 reads per sample

pAlpha = plot_richness(HC.rarefied,
                       x = "SampleID_norep", # if you include your group here, then you can compare
                       #shape = "",
                       #color = "rhizosphere_bulk",
                       measures = c("Observed", "Shannon"),
                       #can also calculate c("Chao1", "InvSimpson")
                       title = "Alpha Diversity")

pAlpha + 
  #geom_point()+
  geom_boxplot(outlier.shape = NA) + 
  #geom_jitter(alpha = 0.2)+ 
  theme_classic() + 
  labs(x = "Bulk Soil Type", y = "Alpha Diversity") + 
 # scale_x_discrete(labels = c("Control","Pit")) + 
  theme(axis.text.x = element_text(size=12), 
        axis.text.y = element_text(size=12), 
        axis.title.y = element_text(size=14), 
        axis.title.x = element_text(size=14))

# estimate richness
rich <- estimate_richness(HC.rarefied)

# Statistical Assessment for Richness
#sample_data(HC)$SampleID_norep <- factor(sample_data(HC)$SampleID_norep, levels = c("C22","CD22","Hb22","HbD22","P22"))
sample_data(HC.rarefied)$SampleID_norep <- factor(sample_data(HC.rarefied)$SampleID_norep, levels = c("C22","Hb22"))
# Five Levels - 
leveneTest(rich$Observed ~ sample_data(HC.rarefied)$SampleID_norep) # Variance homogeneous (p > 0.05)
# greater than 0.05 is good - data is homogeneous
shapiro.test(rich$Observed) # Richness not normally distributed (p < 0.05)
# p >0.05 is good - richness is normally distributed
#m <- aov(rich$Observed ~ sample_data(HC)$SampleID_norep)
#shapiro.test(m$residuals)
# Residuals normally distributed (p > 0.05)
# Other diagnostics - learn more here (https://data.library.virginia.edu/diagnostic-plots/)
#plot(m) # Click in the console and hit Return to see each diagnostic plot
#summary(m) # not significant
# If levenne test or shapiro test < 0.05, then use: 
#kruskalTest(rich$Observed ~ sample_data(HC)$SampleID_norep)
#ans <- kwAllPairsNemenyiTest(rich$Observed ~ sample_data(HC)$SampleID_norep)
#summary(ans) # none are significant
t.test(rich$Observed ~ sample_data(HC.rarefied)$SampleID_norep, alternative = "two.sided")
# t = 1.1113, df = 2.772, p-value = 0.3535 not significant
with(rich, shapiro.test(Observed[sample_data(HC.rarefied)$SampleID_norep == "C22"]))# p = 0.6167 (normal)
with(rich, shapiro.test(Observed[sample_data(HC.rarefied)$SampleID_norep == "Hb22"]))# p = 0.6068 (normal)
var.test(Observed ~ sample_data(HC.rarefied)$SampleID_norep, data = rich) # 0.3344
# t-test is fine to use 

# SDJ calculate the median richness for the two samples 
# Observed
median(rich$Observed[1:3]) #500
sd(rich$Observed[1:3]) #54
median(rich$Observed[4:6]) #469
sd(rich$Observed[4:6]) #24
# Shannon
median(rich$Shannon[1:3]) #5.6
sd(rich$Shannon[1:3]) #0.08
median(rich$Shannon[4:6]) #5.6
sd(rich$Shannon[4:6]) #0.53


# Is beta diversity different for these soils?
# PCoA plot using Bray Curtis distance
bray_dist <- phyloseq::distance(HC, method="bray", weighted= F)
ordination <- ordinate(HC, method="PCoA", distance=bray_dist)
#ordination.nmds <- ordinate(year1, method="NMDS", distance="bray") # insufficient data here

plot_ordination(HC, ordination, color="SampleID_norep") +
  theme(aspect.ratio=1) + 
  theme_classic()

# PERMANOVA - test for differences in centroid of different pCoA hull groups
adon <- adonis2(bray_dist ~ sample_data(HC)$SampleID_norep, method = "bray") # centroid location, can include an interaction term here if needed
adon # not significant


# PERMDISP - multivariate version of Levenne Test. Difference in dispersion for each factor level
m1 <- betadisper(bray_dist, sample_data(HC)$SampleID_norep)
anova(m1) # Dispersion is not different (p = 0.4363)
#scores(m1)




# Taxonomic Analyses
# SIMPER (list how much each ASV contributes to dissimilarity among groups)
OTU <- as.data.frame(HC@otu_table)
# # get rid of rows (ASVs) that are all zero
OTU_noZero <- OTU[rowSums(OTU[]) > 0, ]
sim <- simper(t(OTU_noZero), 
              sample_data(HC)$SampleID_norep)
s <- summary(sim)
# Let's look at the top 20 contributing to dissimilarity between the two groups
head(s$C22_Hb22, n = 20)

# average is the proportion contribution, cumsum is cumulative, ava and avb are mean sequence abundances per group.

# MULTIPATT (list ASVs associated with each group)
# difference between this and SIMPER is in what it is reporting
# MULTIPATT doesn't give a value for every OTU. It only outputs important OTUs
set.seed(1223) # For reproducibility
mp <- multipatt(t(OTU_noZero), 
                sample_data(HC)$SampleID_norep, 
                func = "IndVal.g", 
                control = how(nperm=999))
summary(mp) # no OTUs are important for differentiating the two soil types


