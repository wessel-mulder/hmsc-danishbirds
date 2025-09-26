input <- '.'
library(ggplot2)
library(dplyr)
# --> LOAD TRAITS 
Tr <- readxl::read_xlsx(file.path(input,"data/0_traits/birdTraits-merged-all.xlsx"))
species_list <- read.csv(file.path(input,'data/1_preprocessing/Tr_aits/traits-guild_migration.csv'))

species <- gsub('_',' ',species_list$latin_DOF_underscores)
Tr$latin_DOF[Tr$latin_DOF == 'Acanthis flammea cabaret'] <- 'Acanthis flammea'
Tr_species <- Tr[Tr$latin_DOF %in% species,]

Tr_sub <- Tr_species %>%
  select(english_DOF,latin_DOF,
         Migration_a3_DOF_description,Migration_source,
         foraging_guild_consensus,foraging_guild_source,
         )

setdiff(species,Tr_sub$latin_DOF)
setdiff(Tr_sub$latin_DOF,species)

write.csv(Tr_sub,
          file.path(input,'data/1_preprocessing/Tr_aits/guilds-strategies-sources.csv'))
openxlsx::write.xlsx(Tr_sub,
                     file.path(input,'data/1_preprocessing/Tr_aits/guilds-strategies-sources.xlsx'))

          