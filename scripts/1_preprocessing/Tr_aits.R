rm(list = ls())
# GETTING STARTED ---------------------------------------------------------
library(readxl)
library(dplyr)

# LOADING DATA ------------------------------------------------------------
traits_full <- read_xlsx('data/0_traits/birdTraits-merged-all.xlsx',
                         sheet = 'birdTraits-merged')
head(traits_full)
# select only relevant species 
traits_subset <- traits_full[ rowSums(traits_full[, c("not_breeding", "introduced", "rare_quant", "domestic")] == 1, na.rm = TRUE) == 0, ]
# remove the subspecies 
traits_subset <- traits_subset[!(traits_subset$Species2_AVONET %in% c('Motacilla alba yarrellii','Motacilla flava flavissima')), ]
traits_subset$Species2_AVONET[traits_subset$Species2_AVONET == "Acanthis flammea cabaret"] <- "Acanthis flammea"

# underscores to spaces
traits_subset$latin_DOF_underscores <- gsub(' ','_',traits_subset$latin_DOF)


# MAKE TRAITS FILE --------------------------------------------------------
# grab columns of interest 
traits_whack <- traits_subset %>%
  select(latin_DOF_underscores,Migration_a3_DOF,foraging_guild_consensus)

# check for errors 
unique(traits_whack$foraging_guild_consensus)


# turn migratory to character
traits_whack$Migration_a3_DOF[traits_whack$Migration_a3_DOF==1] <- 'sedentary'
traits_whack$Migration_a3_DOF[traits_whack$Migration_a3_DOF==1.5] <- 'sedentary and short-distance'
traits_whack$Migration_a3_DOF[traits_whack$Migration_a3_DOF==2] <- 'short-distance'
traits_whack$Migration_a3_DOF[traits_whack$Migration_a3_DOF==2.5] <- 'short-and long-distance'
traits_whack$Migration_a3_DOF[traits_whack$Migration_a3_DOF==3] <- 'long-distance'

write.csv(traits_whack,'data/1_preprocessing/Tr_aits/traits-guild_migration.csv')


# MAKE PHYLO FILE ---------------------------------------------------------
# grab columns of interest 
phylo_whack <- traits_subset %>%
  select(latin_DOF_underscores,Species2_AVONET,Family2_AVONET,Order2_AVONET)

phylo_whack$genus <- sapply(strsplit(phylo_whack$Species2_AVONET,' '), function(x) x[1])
write.csv(phylo_whack,'data/1_preprocessing/Taxonomy/taxonomy.csv')


