# The purpose of this code is to compare the soil properties for control plots vs. healthy soils

# Load libraries
library(tidyverse)
library(ggplot2)

# Functions
se <- function(x) {sd(x,na.rm=TRUE)/sqrt(length(x))}

# Load data
soils <- read.csv("soil_properties_control_ref.csv", header = TRUE)

# Compare reference and bulk soils - visualizations
# pH
soils %>%
  ggplot() +
  geom_boxplot(mapping = aes(x = Sample_Type, y=pH), outlier.shape = NA) +
  geom_jitter(mapping = aes(x = Sample_Type, y=pH), alpha = 0.5) +
  theme_classic()

# Bulk Density
soils %>%
  ggplot() +
  geom_boxplot(mapping = aes(x = Sample_Type, y=Bulk_density), outlier.shape = NA) +
  geom_jitter(mapping = aes(x = Sample_Type, y=Bulk_density), alpha = 0.5) +
  theme_classic()

# Sand
soils %>%
  ggplot() +
  geom_boxplot(mapping = aes(x = Sample_Type, y=sand), outlier.shape = NA) +
  geom_jitter(mapping = aes(x = Sample_Type, y=sand), alpha = 0.5) +
  theme_classic()

# Silt
soils %>%
  ggplot() +
  geom_boxplot(mapping = aes(x = Sample_Type, y=silt), outlier.shape = NA) +
  geom_jitter(mapping = aes(x = Sample_Type, y=silt), alpha = 0.5) +
  theme_classic()

# Clay
soils %>%
  ggplot() +
  geom_boxplot(mapping = aes(x = Sample_Type, y=clay), outlier.shape = NA) +
  geom_jitter(mapping = aes(x = Sample_Type, y=clay), alpha = 0.5) +
  theme_classic()

# Organics
soils %>%
  ggplot() +
  geom_boxplot(mapping = aes(x = Sample_Type, y=organics), outlier.shape = NA) +
  geom_jitter(mapping = aes(x = Sample_Type, y=organics), alpha = 0.5) +
  theme_classic()

# calculate the mean, standard error, and range
summary_table <- soils %>%
  group_by(Sample_Type) %>%
  summarise(mean_pH = mean(pH),
            se_pH = se(pH),
            min_pH = min(pH),
            max_pH = max(pH),
            mean_BD = mean(Bulk_density),
            se_BD = se(Bulk_density),
            min_BD = min(Bulk_density),
            max_BD = max(Bulk_density),
            mean_sand = mean(sand),
            se_sand = se(sand),
            min_sand = min(sand),
            max_sand = max(sand),
            mean_silt = mean(silt),
            se_silt = se(silt),
            min_silt = min(silt),
            max_silt = max(silt),
            mean_clay = mean(clay),
            se_clay = se(clay),
            min_clay = min(clay),
            max_clay = max(clay),
            mean_organics = mean(organics),
            se_organics = se(organics),
            min_organics = min(organics),
            max_organics = max(organics))
# Statistics - unequal group sizes, so Welch's t-test if normally distributed, but I have low sample numbers so Wilcoxon test might be required
shapiro.test(soils$pH) # not normally distributed
shapiro.test(soils$Bulk_density)
shapiro.test(soils$sand)
shapiro.test(soils$silt)
shapiro.test(soils$clay)
shapiro.test(soils$organics)

#pH
#parametric test
with(soils, t.test(pH[Sample_Type == "degraded"], pH[Sample_Type == "reference"], alternative = "two.sided"))
# t = -1.7931, df = 4.3677, p-value = 0.1414
# non-parametric test
with(soils, wilcox.test(pH[Sample_Type == "degraded"], pH[Sample_Type == "reference"], alternative = "two.sided", exact = FALSE))
# W = 2, p-value = 0.136
# the p-value is a normal approximation (not an exact p-value)

#Bulk Density
#parametric test
with(soils, t.test(Bulk_density[Sample_Type == "degraded"], Bulk_density[Sample_Type == "reference"], alternative = "two.sided"))
# t = -2.2306, df = 3.0551, p-value = 0.1103

# sand
with(soils, t.test(sand[Sample_Type == "degraded"], sand[Sample_Type == "reference"], alternative = "two.sided"))
# t = -1.8876, df = 3.1832, p-value = 0.1502

# silt
with(soils, t.test(silt[Sample_Type == "degraded"], silt[Sample_Type == "reference"], alternative = "two.sided"))
# t = 0.72887, df = 2.5992, p-value = 0.5262

# clay
with(soils, t.test(clay[Sample_Type == "degraded"], clay[Sample_Type == "reference"], alternative = "two.sided"))
# t = 4.13, df = 5.8259, p-value = 0.00655

# organics
with(soils, t.test(organics[Sample_Type == "degraded"], organics[Sample_Type == "reference"], alternative = "two.sided"))
# t = 2.7976, df = 2.8407, p-value = 0.07237

