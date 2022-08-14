###############################################################################
###############################################################################
# cull-probes.R
# Written November 2020

# Runs like `Rscript cull-probes.R --options` in terminal 
# (from containing folder)

# Eliminates sites in methylation arrays if...
# 1. Detection p-value is not low enough
# 2. Polymorphic target
# 3. Cross hybridized (CpG site or otherwise)

# Keeps track of why a probe was eliminates (for further analysis/visualization)

# Inputs: none

# Outputs:
# 1. M-filtered.csv: M values that are subset (filtered) to only keep those with 
#    "good" probes/sites
# 2. beta-filtered.csv: same but with beta
# 3. probe-fate.csv: data frame that records which tests each site passed

# Takes ~45 seconds to run on local machine

###############################################################################
###############################################################################


# Setup -------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(minfi)
  library(IlluminaHumanMethylationEPICmanifest)
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  library(data.table)
  library(argparse)
})


MAF_cutoff <- 0.05
detectP_cutoff <- 0.01

# File names ------------------------------------------------------------------

poly.file <- paste0(args$probe_fol, "pub-2016McCartney-PolymorphicTargets.txt")
cross.file <- paste0(args$probe_fol, "pub-2016McCartney-CrossHybridCpG.txt")
crossother.file <- paste0(args$probe_fol, "pub-2016McCartney-CrossHybridNonCpG.txt")

detectP.file <- paste0(args$probe_fol, "detectP.csv")
annotation.file <- paste0(args$probe_fol, "annotation.csv")

M.file <- "dataDerived/M.csv"
M2.file <- "dataDerived/M.filtered.csv"

beta.file <- "../dataDerived/beta.csv"
beta2.file <- "../dataDerived/beta.filtered.csv"


# Read in data ------------------------------------------------------------

poly.df <- fread(poly.file)
cross.df <- fread(cross.file, col.names = "CpG_id")
crossother.df <- fread(crossother.file, col.names = "CpG_id")

detectP.df <- fread(detectP.file)
annotation.df <- fread(annotation.file) %>%
  mutate(probes_on_sex_chr = ifelse(chr %in% c("chrX", "chrY"), T, F)) %>% 
  mutate(CpG_SNP = ifelse(CpG_maf > 0.05, T, F)) %>%
  mutate(SBE_SNP = ifelse(SBE_maf > 0.05, T, F))





# Data processing ---------------------------------------------------------
# At each CpG, count the number of patients whose detection p values are above
# cutoff -- if there are one or more, call it a bad detection
detectP.df$is_bad_detection <- (apply(select(detectP.df, -"CpG_id") >= 
                                        args$detectP_cutoff, 1, sum) > 0)


# We'll consider the max minor allele frequency among the different ethnic groups
# Standard in the field is to drop those above 0.05
poly.df$max_MAF <- poly.df %>%
  select((contains("AF"))) %>%
  apply(MARGIN = 1, FUN = max)


# Pull vectors
SNPs <- annotation.df %>%
  filter(CpG_SNP | SBE_SNP) %>%
  pull(CpG_id)

# Vector of positions on sex chromosome
sex_chr <- annotation.df %>% 
  filter(probes_on_sex_chr) %>% 
  pull(CpG_id)

# Vector of 
polymorphisms <- poly.df %>%
  filter(max_MAF >= args$MAF_cutoff) %>%
  pull(IlmnID)

# Cobble all the information into a dataframe to eventually write out
probe_fate.df <- detectP.df %>%
  select(c(CpG_id, is_bad_detection)) %>%
  mutate(has_SnP = ifelse(CpG_id %in% SNPs, T, F)) %>%
  mutate(on_sex_chr = ifelse(CpG_id %in% sex_chr, T, F)) %>%
  mutate(has_CpG_crosshybrid = ifelse(CpG_id %in% cross.df$CpG_id, T, F)) %>%
  mutate(has_nonCpG_crosshybrid = ifelse(CpG_id %in% crossother.df$CpG_id, T, F)) %>%
  mutate(is_polymorphic = ifelse(CpG_id %in% polymorphisms, T, F)) %>%
  mutate(should_remove = select(., -CpG_id) %>% rowSums()) %>%
  mutate(should_remove = ifelse(should_remove > 0, T, F))



# Write out with a bit of in-processing -----------------------------------



fwrite(probe_fate.df, file = ofile)

good_probes <- probe_fate.df %>% 
  filter(!should_remove) %>% 
  pull(CpG_id)

# Free up RAM
rm(probe_fate.df, annotation.df, cross.df, crossother.df)

# Create a filtered data set of Methylation beta values
fread(beta.file) %>%
  filter(CpG_id %in% good_probes) %>%
  fwrite(beta2.file)

# Ibid for M-values
fread(M.file) %>%
  filter(CpG_id %in% good_probes) %>%
  fwrite(M2.file)

#END