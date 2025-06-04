rm(list = ls())
# GETTING STARTED ---------------------------------------------------------
library(readxl)

# LOADING DATA ------------------------------------------------------------
traits_full <- read_xlsx('data/0_traits/birdTraits-merged-all.xlsx')

# select only relevant species 
traits_subset <- traits_full[ rowSums(traits_full[, c("not_breeding", "introduced", "rare", "domestic")] == 1, na.rm = TRUE) == 0, ]

# underscores to spaces
traits_subset$latin_DOF_underscores <- gsub(' ','_',traits_subset$latin_DOF)

# grab columns of interest 
traits_whack <- traits_subset %>%
  select(latin_DOF_underscores,Migration_AVONET,foraging_guild_consensus)


# WRITE OUT ---------------------------------------------------------------
write.csv(traits_whack,'data/1_preprocessing/Tr_aits/traits-guild_migration.csv')


