###############################################################################
###############################################################################
# compute-methylation.R
# Written October 2020

# runs like `Rscript compute-methylation.R` from container folder in terminal
# Or `Rscript compute-methylation.R --ifolder path/to/idats/` ibid
# Or source("compute-methylation.R") from R studio

# Does
# 1. Computes the beta and M methylation signals given idat files
# 2. Computes detection p-values for each probe

# Required inputs: none
# (unless default paths are wrong in argument parsing step)

# Outputs:
# 1. 'M.csv': methylation signal (log version)
# 2. 'beta.csv': methylation signal ("linear" version)
# 3. 'detect-P.csv': detection p-value for each probe (low is good)
# 4. 'annotation.csv': includes chromosome, position, corresponding genes, etc.

# IlluminaHumanMethylationEPICanno.ilm10b4.hg19

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
  library(argparse)
})

# Parse arguments ---------------------------------------------------------

# Instantiate argument parser a la Python version
parser <- ArgumentParser(description = 'Compute methylation signal given .idat files')

parser$add_argument('--in_fol', 
                    default = '../../data/850Karrays/', 
                    help = 'Relative path to where .idat files are stored')

parser$add_argument('--sheets_fol', 
                    default = '../../data/samplesheets/', 
                    help = 'Samplesheet with info on patient characteristics')

parser$add_argument('--out_fol', default = '../../data/', 
                    help = 'Root directory in which to save outputs')

parser$add_argument('--probe_fol', default = '../../data/probes/', 
                    help = 'Root directory in which to save probe information')

# Make arguments accessible using object$var syntax (object.var in Python)
args <- parser$parse_args()
  


# File names --------------------------------------------------------------

ss.file <- paste0(args$sheets_fol, "/master.csv")
wb.file <- paste0(args$sheets_fol, "/master-WB.csv")

beta.file <- paste0(args$out_fol, "/beta.csv")
M.file <- paste0(args$out_fol, "/M.csv")

detectP.file <- paste0(args$probe_fol, "detectP.csv")
annotation.file <- paste0(args$probe_fol, "annotation.csv")

# Set up output paths, load data ------------------------------------------


create_new_dir <- function(suffix){
  # Takes the global output root folder, appends the 'suffix' input to the
  # end, makes the full directory (if needed), and returns the file name.
  new.folder <- paste0(args$ofolder, suffix)
  
  if (!exists(new.folder)) {
    dir.create(new.folder)
  }
  
  return(new.folder)
}



# Sample sheet with demographic, phenotypic information
ss.df <- read_csv(ss.file) 

# Collect base names
rel.names <- paste(args$in_fol, ss.df$patient_folder, ss.df$patient_id, sep = "/")

# Pull in red-green channel set
rgset <- read.metharray(rel.names)


annotation <- getAnnotation(rgset)



# Pre-processing -----------------------------------------------------------

# Estimate white blood cell counts
wb.cells <- estimateCellCounts(rgset)

# Join with ss.df to create a new master sample sheet
wb.df <- as.data.frame(wb.cells) %>%
  mutate(patient_id = row.names(wb.cells)) %>%
  full_join(ss.df, by = "patient_id")

# Save amplesheet with white blood counts
write_csv(wb.df, wb.file)

# Detection P-values (computed for each probe for each patient)
detectP.df <- as.data.frame(detectionP(rgset))
detectP.df$CpG_id <- row.names(detectP.df)

# Save name "./probes/detect-P.csv" or something similar
fwrite(detectP.df, detectP.file)

# Memory management
rm(detectP.df)

# One pre-processing step
gen.ratio <- preprocessQuantile(rgset)
rm(rgset)

# Measure the methylation signal (beta version)
beta.df <- as.data.frame(getBeta(gen.ratio))
beta.df$CpG_id <- row.names(beta.df) 

fwrite(beta.df, beta.file)

rm(beta.df)

# M methylatation signal
M.df <- as.data.frame(getM(gen.ratio))
M.df$CpG_id <- row.names(M.df)

fwrite(M.df, M.file)
rm(M.df)

# This returns a data frame that can be subset
notes.df <- as.data.frame(getAnnotation(gen.ratio))
notes.df$CpG_id <- row.names(notes.df) # redundant, but better name

fwrite(notes.df, annotation.file)






# Cleanup -----------------------------------------------------------------

print("compute_mehtylation.R took ")
print(Sys.time() - start)


rm(list = ls())
