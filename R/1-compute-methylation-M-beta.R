###############################################################################
###############################################################################
# compute-methylation-M-beta.R
# Written October 2020, cleaned up August 2022

# runs like `Rscript compute-methylation-M-beta.R` 

# Does
# 1. Computes the beta and M methylation signals given idat files
# 2. Computes detection p-values for each probe

# Outputs:
# 1. 'M.csv': methylation signal (log version)
# 2. 'beta.csv': methylation signal ("linear" version)
# 3. 'detect-P.csv': detection p-value for each probe (low is good)
# 4. 'annotation.csv': includes chromosome, position, corresponding genes, etc.

# On Ubuntu 20, R 3.6.3, 16 GB RAM, takes 7-9 mins
###############################################################################
###############################################################################

start <- Sys.time()

# Libraries ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(IlluminaHumanMethylationEPICmanifest)
  library(limma)
})


# File names --------------------------------------------------------------

ss.ifile <- "../dataRaw/master.csv"
beta.ofile <- "../dataDerived/beta.csv"
M.ofile <- "../dataDerived/M.csv"
detectP.ofile <- "../dataDerived/detectP.csv"
annotation.ofile <- "../dataDerived/annotation.csv"


# Sample sheet with demographic, phenotypic information
ss.df <- read_csv(ss.ifile)

# Collect base names
rel.names <- file.path("../dataRaw/850Karrays/", ss.df$patient_folder, ss.df$patient_id)

# Pull in red-green channel set
rgset <- read.metharray(rel.names)
gc()

annotation <- getAnnotation(rgset)
gc()



# Pre-processing -----------------------------------------------------------

# Detection P-values (computed for each probe for each patient)
detectP.df <- as.data.frame(detectionP(rgset))
detectP.df$CpG_id <- row.names(detectP.df)

# Save name "./probes/detect-P.csv" or something similar
fwrite(detectP.df, detectP.ofile)

# Memory management
rm(detectP.df); gc()

# One pre-processing step
gen.ratio <- preprocessQuantile(rgset)
rm(rgset)

# Measure the methylation signal (beta version)
beta.df <- as.data.frame(getBeta(gen.ratio))
beta.df$CpG_id <- row.names(beta.df) 
fwrite(beta.df, beta.ofile)
rm(beta.df); gc()


# M methylatation signal
M.df <- as.data.frame(getM(gen.ratio))
M.df$CpG_id <- row.names(M.df)
fwrite(M.df, M.ofile)
rm(M.df); gc()

# This returns a data frame that can be subset
notes.df <- as.data.frame(getAnnotation(gen.ratio))
notes.df$CpG_id <- row.names(notes.df) # redundant, but better name
fwrite(notes.df, annotation.ofile)


# Cleanup -----------------------------------------------------------------

print("Script took ")
print(Sys.time() - start)

rm(list = ls())
#END