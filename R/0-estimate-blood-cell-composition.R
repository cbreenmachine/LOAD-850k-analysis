library(minfi)
library(tidyverse)
library(limma)
library(methylCC)

ss.df <- read_csv("../dataRaw/master.csv")
rel.names <- file.path("../dataRaw/850Karrays/", ss.df$patient_folder, ss.df$patient_id)

# Pull in red-green channel set
rgset <- read.metharray(rel.names)

gc()

# Estimate blood cell composition using Hansen's method
minfi.counts <- minfi::estimateCellCounts(rgset)
gc()

methylCC.counts <- cell_counts(methylCC::estimatecc(object = rgset))
gc()


# Munge and put together
df.1 <- data.frame(minfi.counts) %>% 
  rownames_to_column("patient_id") %>% 
  pivot_longer(cols = -patient_id, names_to = "cell_type", values_to = "minfi_prop")

df.2 <- data.frame(methylCC.counts) %>% 
  rownames_to_column("patient_id") %>% 
  pivot_longer(cols = -patient_id, names_to = "cell_type", values_to = "methylCC_prop")


df <- full_join(df.1, df.2, by = c("patient_id", "cell_type"))


write_csv(df, "../dataDerived/blood-cell-composition.csv")
# 
# p <- df %>% 
#   ggplot(aes(x = minfi_prop, y = methylCC_prop, color = cell_type)) +
#   geom_point(size = 1.5) +
#   theme_minimal() +
#   ggtitle("84 samples") +
#   labs(caption = "")
