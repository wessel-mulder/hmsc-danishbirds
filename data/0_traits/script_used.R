rm(list = ls())
# GETTING STARTED ---------------------------------------------------------
library(readxl)


# LOAD DATA ---------------------------------------------------------------
# load traits 
traits <- read_xlsx('data/0_traits/birdTraits-merged-all.xlsx')

# subset the weird species out 
traits_subset <- traits[ rowSums(traits[, c("not_breeding", "introduced", "rare", "domestic")] == 1, na.rm = TRUE) == 0, ]

# select traits im interested 
traits_subset_select <- traits_subset[,3:ncol(traits_subset)] %>%
  select(-Species2_AVONET,-Family2_AVONET,-Order2_AVONET)

# load taxonomy 
taxon <- read.csv('data/0_taxonomy/IOC_Names_File_Plus-15.1_wide.csv')
taxon <- taxon %>% 
  select(-X)

names(taxon) <- c('order_IOC','family_IOC','genus_IOC','species_IOC','ssp_IOC')

# MERGE BASED ON NUMBER OF NAMES  ---------------------------------------------------
# Split DOF into two groups
taxon_2 <- taxon[,1:4]
taxon_2 <- distinct(taxon_2)

traits_2 <- traits_subset_select %>%
  filter(str_count(latin_DOF, "\\S+") == 2) %>%
  mutate(species_IOC = latin_DOF)

traits_3 <- traits_subset_select %>%
  filter(str_count(latin_DOF, "\\S+") == 3) %>%
  mutate(ssp_IOC = latin_DOF)

# Join both
# Join taxonomy info from df into english_DOF subsets
joined_2 <- traits_2 %>%
  left_join(taxon_2, by = "species_IOC")
names(joined_2)

joined_3 <- traits_3 %>%
  left_join(taxon, by = "ssp_IOC")
names(joined_3)

# Combine back into one dataframe
final_joined <- bind_rows(joined_2, joined_3)

# reorder
final_joined <- final_joined %>%
  select(
    1, 2,                    # keep col1, col2
    (ncol(final_joined)-4):ncol(final_joined),  # last 5 columns
    3:(ncol(final_joined)-4)  # columns 3 up to just before last 4 columns
  )

write.csv(final_joined,
          'data/0_traits/birdTraits-merged-subset.csv',
          row.names = F)
metaGenerator('data/0_traits/')
