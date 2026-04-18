# Alpha and Beta Diversity Analysis
# Following starter code from Claire Winfrey
# Sierra updated code from 12.15.22 on 1.10.23 and nothing really changed except the mapping file now has more columns
# SDJ re-ran this code for Claire's samples on 2/7/23

#### Setup ####
# We'll use phyloseq for the first time 
library(tidyverse)
library(phyloseq) # For all bioinformatics analyses post Dada2
packageVersion('phyloseq')
library(ggplot2)
library(vegan)
library(car) # For Levene Test and anova
library(PMCMRplus) # For Nemenyi posthoc test
#library("DESeq2")
library(indicspecies) # For indicator species
library(microbiome)

#### Functions ####
# 1. Function to get the taxonomy table out of phyloseq:
# (Inspired by similar function at https://jacobrprice.github.io/2017/08/26/phyloseq-to-vegan-and-back.html)
taxtable_outta_ps <- function(physeq){ #input is a phyloseq object
  taxTable <- tax_table(physeq)
  return(as.data.frame(taxTable))
}

### Load the ASV table (called OTU in phyloseq)
list.files()
otumat <- read.table(file = "seqtab_wTax_mctoolsr.txt",header=T)
dim(otumat) # 52 samples and 7096 ASVs
otumat <- otumat %>%
  tibble::column_to_rownames("X.ASV_ID") # push the ASV_ID column into the rownames position

### Load the taxonomy table - this comes from two columns in the dada2 output
# Find the ASV names
asv <- rownames(otumat)
# Find the taxonomy info
tax <- otumat$taxonomy
# Combine them into a dataframe
taxmat <- as.data.frame(cbind(asv, tax)) 
# change the column names to be "ASV_ID" and taxonomy
colnames(taxmat) <- c("ASV_ID", "taxonomy")
# split the column into each taxonomy
taxmat <- separate(taxmat, 
                   col = taxonomy, 
                   into= c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";")
# Format tax_sep for phyloseq by making ASV_ID the rownames
taxmat <- taxmat %>%
  tibble::column_to_rownames("ASV_ID")
# Delete the species column because it is ASV_ID
taxmat <- taxmat[, 1:6]
# remove the taxonomy column from the OTU table 
otumat <- otumat[, -51]

#### Load the data mapping (metadata) file
mapping <- read.table("pits_mapping_12202022.txt", header = TRUE) # all metadata
colnames(mapping)
mapping <- remove_rownames(mapping)
sampledata <- column_to_rownames(mapping, var = "SampleID") # rownames should be the sample IDs

# Transform to phyloseq objects
# the sample dataframe includes 4 extra samples because we took more replicates of pH measurements than we took samples for sequencing. So, I add delete those columns
sampledata <- sampledata[1:50, ]
# check that the names match
all.equal(sort(colnames(otumat)), sort(rownames(sampledata))) 
# Transform otu table and taxonomy tables into matrices
otumat <- as.matrix(otumat)
taxmat <- as.matrix(taxmat)
#View(otu_mat)
OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
samples = sample_data(sampledata)
pits_16S_raw.ps <- phyloseq(OTU, TAX, samples)
pits_16S_raw.ps

#### Check on the naming
colnames(otu_table(pits_16S_raw.ps)) # names associated with the otu table
sample_names(pits_16S_raw.ps) #names associated with the mapping table
sample_variables(pits_16S_raw.ps)
sample_sums(pits_16S_raw.ps) #returns the ASV count for each sample

#### Get Phylum counts ####
rank_names(pits_16S_raw.ps)
table(tax_table(pits_16S_raw.ps)[, "Phylum"], exclude = NULL)
# probably want to remove the NA phyla and the phyla with only 1 feature
pits_16S_rmNAPhyla.ps <- subset_taxa(pits_16S_raw.ps, !is.na(Phylum))

## If you need to select particular samples to analyze (like only the rhizosphere samples)
#pits_16S_raw.ps <- subset_samples(pits_16S_raw.ps, rhizosphere_bulk =="Yes")
#pits_16S_raw.ps

#first visualization to make sure things look right...this takes time so only if you need it 
#plot_bar(pits_16S_raw.ps, fill = "Phylum")

#### Prevalence ####
# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(pits_16S_rmNAPhyla.ps),
               MARGIN = ifelse(taxa_are_rows(pits_16S_rmNAPhyla.ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(pits_16S_rmNAPhyla.ps),
                    tax_table(pits_16S_rmNAPhyla.ps))
head(arrange(prevdf,-TotalAbundance))
# note, there are 50 samples, but two are blanks. So the prevalence column gives me a count of the number of samples which have that ASV. We do not really have a reason to limit the prevalence of ASVs in the analysis. 
# note, the Total Abundance column is "the sum of all reads observed for each ASV". When people say that they removed rare species which had < 1% abundance...they are talking about an abundance cutoff. The phyloseq manual uses a prevalence cutoff of 5%
prevalenceThreshold = 0.05 * nsamples(pits_16S_rmNAPhyla.ps)
prevalenceThreshold # so they have to be in at least 2.5 samples to be included in the analysis

#Are there phyla that are comprised of mostly low-prevalence features? Compute the total and average prevalence of the features in each phylum.
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
# now prune those with a prevalence < 5%
keepTaxa = rownames(prevdf)[(prevdf$Prevalence >= prevalenceThreshold)]
pits_16S_rmNAPhyla_prevalent.ps = prune_taxa(keepTaxa, pits_16S_rmNAPhyla.ps)

# check your work
table(tax_table(pits_16S_rmNAPhyla_prevalent.ps)[, "Phylum"], exclude = NULL)
# this table is much smaller, I think it tells us the number of ASVs in each phylum

##################################################################
# FILTER OUT MITOCHONDRIA, CHLOROPLASTS, and KINGDOM EUKARYOTA
# chloroplast is an order and mitochondria is a family in this version of SILVA
##################################################################
# Get taxonomy table out of phyloseq:
pits_taxTable <- psmelt(pits_16S_rmNAPhyla_prevalent.ps)
#View(SRS_taxTable)

# First check what will be removed
tax_filt <- pits_taxTable %>%
  filter(Order != "Chloroplast") %>%
  filter(Family != "Mitochondria") %>%
  filter(Kingdom != "Eukaryota") %>%
  filter(Kingdom != "NA") %>% #Remove ASVs where kingdom is unknown ("NA" in first column)
  filter(Phylum != "NA") #Remove ASVs where phylum is unknown ("NA" in first column)

# View(tax_noeuksorNAs)
dim(pits_taxTable)[1] - dim(tax_filt)[1]
# There are 106,000 ASVs in the dataset
# There are 103,600 ASVs staying in the dataset
# We will filter out a total of 2400 taxa
# This means that at 5% prevalence, we removed 2.26% of ASVS and kept 97.7% of ASVs

# Now actually remove them from the pyloseq object. First find and save each set as their own object. We will then remove all the ones we do not want below. 
# Chloroplasts (order)
chloros <- pits_taxTable %>%
  filter(Order == "Chloroplast")
dim(chloros) #200 chloroplast ASVs
chloro_names <- rownames(chloros)

# Mitochondria (family)
mitos <- pits_taxTable %>%
  filter(Family == "Mitochondria")
dim(mitos) #850 mitochondria ASVs
mito_names <- rownames(mitos)

# Eukaryota (kingdom)
kingdomEuks <- pits_taxTable %>%
  filter(Kingdom == "Eukaryota")
dim(kingdomEuks) #0
euks_names <- rownames(kingdomEuks)

# NA's (kingdom) - these were removed above as NA's
# kingdomNAs <- pits_taxTable %>%
#   filter(Kingdom == "NA")
# dim(kingdomNAs) #0
# kNAs_names <- rownames(kingdomNAs)
# NA's (phyla)
PhylumNAs <- pits_taxTable %>%
  filter(Phylum == "NA")
dim(PhylumNAs) #1350
pNAs_names <- rownames(PhylumNAs)

remove_ASVs <- c(chloro_names, mito_names, euks_names,pNAs_names)
length(remove_ASVs) # removing a total of 2400 ASVs

# Removing ASVs in phyloseq comes from: Joey711's code: https://github.com/joey711/phyloseq/issues/652
# Remove the ASVs that identified:
all_Taxa <- taxa_names(pits_16S_rmNAPhyla_prevalent.ps) #get all tax names in original, uncleaned dataset
ASVstoKeep <- all_Taxa[!(all_Taxa %in% remove_ASVs)]
length(ASVstoKeep) #2120
pits_16S_filt_ps <- prune_taxa(ASVstoKeep, pits_16S_rmNAPhyla_prevalent.ps) # new phyloseq object with just the taxa we want!

##################################################################
# Checking the data
# Code here comes from Cliff Tutorial
##################################################################
# again calculate the number of ASV per sample
sort(phyloseq::sample_sums(pits_16S_filt_ps)) 
# minimum reads in blanks = 110
# maximum reads in blanks = 178
# minimum reads in a sample = 4,720
# maximum reads in a sample = 56,395
mean(phyloseq::sample_sums(pits_16S_filt_ps))
# mean reads per sample = 24,065.9

# How are number of reads distributed across the samples?
seqcounts <- as.data.frame(sort(colSums(otu_table(pits_16S_filt_ps)))) %>%
  rename("seqs" = "sort(colSums(otu_table(pits_16S_filt_ps)))") %>%
  rownames_to_column(var = "sampleID")

# Now we have a dataframe with two columns, seqs and sampleID which we can plot
ggplot(seqcounts, aes(reorder(sampleID, seqs, mean), seqs)) + # Dataframe and variables
  geom_bar(stat = "identity") + # Type of graph
  labs(y = "# Reads", x = "Sample") + # Axes labels
  coord_flip() + # Flip axes
  geom_hline(yintercept = 4500, color = "blue") + # this seems backwards because of the coord_flip() command above
  theme_classic() + 
  theme(axis.text.y = element_text(size = 2)) # rarefaction at 9750 looks ok to me because it makes the blanks drop out but keeps as many possible reads as we can. 

# Rarefaction - decided to skip rarefaction for now
# Plot rarefaction curves - this isn't working
# rarecurve(t(otu_table(pits_16S_filt_ps)), sample = 9750, step=50, cex = 0.6)

# Look at the barplot or rarecurve to decide where you want to rarefy to. If you go too high, you will drop samples. If you go too low, you lose a lot of data unnecessarily
# set.seed(19)
# pits_16S_rar_filt_ps <- rarefy_even_depth(pits_16S_filt_ps, sample.size = 9750, replace=FALSE, trimOTUs=TRUE) # the rngseed warning can be ignored because I use "set.seed()" before running the rarefy() function. Phyloseq warns against rarefaction, so I may want to reconsider this step
# I drop both blank samples during this step
# 350 OTUs were removed because they are no longer present in any sample after random sub-sampling
# check the number of reads per sample
#colSums(otu_table(pits_16S_rar_filt_ps)) # did not rarefy this time
colSums(otu_table(pits_16S_filt_ps)) # all read counts are different
sum(colSums(otu_table(pits_16S_filt_ps))) # total number of reads in the dataset is 1,203,295

####################################################################
# Remove Rare Taxa = phyla that have less than 1% abundance
# prevalence - the number of samples in which a taxa appears, and 
# total counts - the total number (or proportion) of observations of a taxa across all samples.
# Samples were filtered based on 5% prevalence above, so I am not removing rare taxa anymore.
####################################################################
# choose the abundance level across the whole dataset which is to be the cut off for dropping rare taxa
# minTotRelAbun = 1e-5 # this represents 0.005%
# # get total number of reads of each ASV in the dataset (add the counts for all 6370 ASVs across 50 samples)
# x = taxa_sums(pits_16S_rar_filt_ps)
# # Make a dataframe of ASV_IDs which are to be kept because they exceed the threshold
# keepTaxa = taxa_names(pits_16S_rar_filt_ps)[which((x / sum(x)) > minTotRelAbun)]
# # prune out the unwanted taxa
# pits_16S_rar_filt_pruned_ps = prune_taxa(keepTaxa, pits_16S_rar_filt_ps)
# # Check how many ASVs were removed because they had too low of abundances
# ntaxa(pits_16S_rar_filt_ps)
# ntaxa(pits_16S_rar_filt_pruned_ps)
# I threw out 6353-4255 = 2098 taxa

#second visualization to make sure things look right
#plot_bar(pits_16S_filt_ps, fill = "Phylum")
# Alternatively
#ps.phylum = tax_glom(pits_16S_filt_ps, taxrank="Phylum", NArm=FALSE)
#plot_bar(ps.phylum, fill="Phylum")

####################################################################
# Save these Pre-Processing Steps so that you do not have to do them again
#save(pits_16S_filt_ps, file = "pits_processed_filt_ps.RData")
#saveRDS(pits_16S_filt_ps, file = "pits_processed_filt_ps.rds")
####################################################################

#### Check the ASVs that are in the blanks ###
# subset for the blank samples
samplesToKeep <- c("blank1","blank2")
blanks <- prune_samples(sample_names(pits_16S_filt_ps) %in% samplesToKeep , pits_16S_filt_ps)
colSums(otu_table(blanks)) # number of reads in these samples
table(tax_table(blanks)[, "Phylum"], exclude = NULL) # 28 phyla represented, but these have zero abundance in these samples. So actually...
melt_blanks <- psmelt(blanks)
# remove ASVs from the table with Abundance = 0
temp <- melt_blanks %>%
  filter(Abundance > 0)
# now there are only three ASVs showing up ... ~100 reads from Bradyrhizobium in each blank and a Bacteroidia in blank 1


