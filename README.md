# Supporting Information

This repository contains all code for bioinformatic analysis for the manuscript, Karban, C. C., Gugel, S. D., & Barger, N. N. (In press). Microsite creation increases native plant density and biomass in a semi-arid grassland restoration. Ecological Applications.


# What can I do with this code?

In publishing this repository, our hope is that this code is useful to other members of the scientific community. This repository is released under a Creative Commons BY (CC-BY) license, which means that all code published here can be shared and adapted for any purposes so long as appropriate credit and citation of the original paper is given.


# How do I run this code?

1. Go to NCBI and download the raw sequence files. 
2. Download and install R for your operating system.
3. Download and install RStudio for your operating system.
4. Download a zip file of this repository and decompress it in a directory of your choosing on your computer.
5. Navigate to the directory and open the Rstudio Project file to load this project's files.
6. Open the script(s) you would like to run. Scripts are numbered in the order they should be executed e.g, 01, 02, 03. You do not have to start with dada2 script (01.processing_bioinformatics.Rmd) or the phyloseq script (02.pits16Sphyloseq_1.10.22.R). The resulting file is included as "pits_processed_filt_ps.rds" so that you can used them in the data analysis scripts (03-05 instead). 
7. Ensure that you have all of the required libraries installed by inspecting the Setup chunks. In these scripts, we note the CRAN/GitHub version/release that was used. If any libraries fail to install, note the name of the library and attempt to manually install its most recent version via CRAN or GitHub.
8. Contact Sierra Gugel for more information or questions: sierra.jech@nau.edu


# Scripts
- 01.processing_bioinformatics.Rmd
- 02.pits16Sphyloseq_1.10.22.R
- 03.soilProperties.R
- 04.biochar_indepth.R
- 05.healthy_degraded_comparison.R

# Output/Input
- pits_processed_filt_ps.rds
- qPCR_mass_extracted.csv
