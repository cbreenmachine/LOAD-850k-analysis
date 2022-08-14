###############################################################################
###############################################################################
# make-bed-files.R
# Written May 2021

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
parser$add_argument('--indir', default = '../../data/850Karrays/', help = 'Relative path to where .idat files are stored')
parser$add_argument('--samplesheets', default = '../../data/samplesheets/', help = 'Samplesheet with info on patient characteristics')
parser$add_argument('--outdir', default = '../../data/BED/', help = 'Root directory in which to save outputs')
parser$add_argument('--probeinfo', default = '../../data/probes/probe-fate.csv', help = 'Good/bad probes.')
parser$add_argument('--offset', default = 0, help = 'Index correction')
args <- parser$parse_args()

# Samplesheets and Probes ------------------------------------------------------
DHMRI.df <- read.csv(paste0(args$samplesheets, "/DHMRI_samplesheet_pilot.csv"), stringsAsFactors = F) 
ss.df <- read.csv(paste0(args$samplesheets, "/master.csv"), stringsAsFactors = F) %>%
  dplyr::filter(adrcnum %in% DHMRI.df$adrcnum)
good_cpgs <- read.csv(args$probeinfo) %>%
                filter(should_remove == FALSE) %>%
                pull(CpG_id)

# Set up output paths, load data ------------------------------------------
rel.names <- paste(args$indir, ss.df$patient_folder, ss.df$patient_id, sep = "/")
rgset <- read.metharray(rel.names)

# Pre-processing -----------------------------------------------------------

# One pre-processing step
gen.ratio <- preprocessQuantile(rgset)
rm(rgset)

annot.df <- getAnnotation(gen.ratio) %>% 
              as.data.frame() %>%
              rownames_to_column(var="CpG_id") %>% 
             transmute(chrom = chr, chromStart = pos + args$offset, chromEnd = pos + 1 , CpG_id, strand) %>%
              filter(CpG_id %in% good_cpgs)

for (i in 1:nrow(ss.df)) {
  ix <- which(colnames(gen.ratio[,i]) == ss.df$patient_id)
  adrc_id <- ss.df$adrcnum[ix]
  
  if (adrc_id %in% DHMRI.df$adrcnum){
    out_name <- DHMRI.df[adrc_id == DHMRI.df$adrcnum, "Alisch_ID"]
    beta.df <- as.data.frame(getBeta(gen.ratio[,i])) %>% rownames_to_column(var="CpG_id") %>% dplyr::rename(Beta = 2)
    M.df <- as.data.frame(getM(gen.ratio[,i])) %>% rownames_to_column(var="CpG_id") %>% dplyr::rename(M = 2)
    annot.df %>%
      merge(beta.df, by = "CpG_id") %>% merge(M.df, by = "CpG_id") %>%
      dplyr::select(c(chrom, chromStart, chromEnd, M, Beta, CpG_id, strand)) %>%
      write.table(file=paste0(args$outdir, "s", out_name, ".bed"), row.names = FALSE, quote=FALSE, col.names = TRUE)
  }
}
# END
