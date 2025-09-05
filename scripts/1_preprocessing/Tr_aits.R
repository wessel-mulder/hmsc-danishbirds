rm(list = ls())
# GETTING STARTED ---------------------------------------------------------
library(readxl)
library(dplyr)

# LOADING DATA ------------------------------------------------------------
traits_full <- read_xlsx('data/0_traits/birdTraits-merged-all.xlsx',
                         sheet = 'birdTraits-merged')

# select only relevant species 
traits_subset <- traits_full[ rowSums(traits_full[, c("not_breeding", "introduced", "rare", "domestic")] == 1, na.rm = TRUE) == 0, ]

# underscores to spaces
traits_subset$latin_DOF_underscores <- gsub(' ','_',traits_subset$latin_DOF)

# grab columns of interest 
traits_whack <- traits_subset %>%
  select(latin_DOF_underscores,Migration_a3_DOF,foraging_guild_consensus)

# check for errors 
unique(traits_whack$foraging_guild_consensus)

# remove the subspecies 
traits_whack <- traits_whack[!(traits_whack$latin_DOF_underscores %in% c('Motacilla_alba_yarrellii','Motacilla_flava_flavissima')), ]
traits_whack$latin_DOF_underscores[traits_whack$latin_DOF_underscores == "Acanthis_flammea_cabaret"] <- "Acanthis_flammea"

# turn migratory to character
traits_whack$Migration_a3_DOF[traits_whack$Migration_a3_DOF==1] <- 'sedentary'
traits_whack$Migration_a3_DOF[traits_whack$Migration_a3_DOF==1.5] <- 'sedentary and short-distance'
traits_whack$Migration_a3_DOF[traits_whack$Migration_a3_DOF==2] <- 'short-distance'
traits_whack$Migration_a3_DOF[traits_whack$Migration_a3_DOF==2.5] <- 'short-and long-distance'
traits_whack$Migration_a3_DOF[traits_whack$Migration_a3_DOF==3] <- 'long-distance'

# WRITE OUT ---------------------------------------------------------------
write.csv(traits_whack,'data/1_preprocessing/Tr_aits/traits-guild_migration.csv')


