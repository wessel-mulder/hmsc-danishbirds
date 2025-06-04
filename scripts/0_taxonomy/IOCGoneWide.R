rm(list=ls())

# GETTING STARTED ---------------------------------------------------------
library(readxl)
library(dplyr)
library(stringr)
library(supportScripts)

# LOADING DATA ------------------------------------------------------------
# t <- read.tree('/Users/bhr597/Downloads/example-hmsc-anthony/Cumulative-thermal-stress-predicts-functional-change/Data/phy/phy_cut.tre')

# read IOC 2024
ioc2024 <- read_xlsx('../data/taxonomy/taxonomy-IOC/IOC_Names_File_Plus-15.1_red.xlsx')
ioc <- ioc2024 %>%
  select(Rank,`Scientific Name`)

# ranks of interest
unique(ioc$Rank)
ranks <- c('ORDER','Family','Genus','Species','ssp')
ioc_relevant <- ioc[ioc$Rank %in% ranks,]
names(ioc_relevant) <- c('rank','name')

# fill order 
library(dplyr)
library(stringr)
library(tidyr)

# Map order
ioc_with_order <- ioc_relevant %>%
  mutate(
    order = if_else(rank == "ORDER", str_remove(name, "ORDER\\s+"), NA_character_),
    order_group = cumsum(rank == "ORDER")
  ) %>%
  group_by(order_group) %>%
  fill(order, .direction = "down") %>%
  ungroup()

# Map family
ioc_with_family <- ioc_with_order %>%
  mutate(
    family = if_else(rank == "Family", str_remove(name, "Family\\s+"), NA_character_),
    family_group = cumsum(rank %in% c("ORDER", "Family"))
  ) %>%
  group_by(family_group) %>%
  fill(family, .direction = "down") %>%
  ungroup()

# map genus
ioc_with_genus <- ioc_with_family %>%
  mutate(
    genus = if_else(rank == "Genus", str_remove(name, "Genus\\s+"), NA_character_),
    genus_group = cumsum(rank %in% c("ORDER", "Family", "Genus"))
  ) %>%
  group_by(genus_group) %>%
  fill(genus, .direction = "down") %>%
  ungroup()

# map species 
ioc_with_species <- ioc_with_genus %>%
  mutate(
    species = if_else(rank == "Species", str_remove(name, "Species\\s+"), NA_character_),
    species_group = cumsum(rank %in% c("ORDER", "Family", "Genus", "Species"))
  ) %>%
  group_by(species_group) %>%
  fill(species, .direction = "down") %>%
  ungroup()

# map ssp 
ioc_with_ssp <- ioc_with_species %>%
  mutate(ssp = if_else(rank == "ssp", name, NA_character_))

# remove NAs
ioc_clean <- ioc_with_ssp %>%
  filter(!is.na(family), !is.na(genus), !is.na(species))

# selected columns
ioc_clean_selected <- ioc_clean %>%
  select(order,family,genus,species,ssp)

# remove caps
ioc_clean_nocaps <- ioc_clean_selected %>%
  mutate(order = str_to_title(order))

ioc_clean_ssp_spelled <- ioc_clean_nocaps %>%
  mutate(
    ssp = if_else(
      !is.na(ssp),
      paste(genus, word(species, -1), word(ssp, -1)),
      NA_character_
    )
  )

# write out 
write.csv(ioc_clean_ssp_spelled,
          '../data/taxonomy/taxonomy-IOC/IOC_Names_File_Plus-15.1_wide.csv')
write.csv(ioc_clean_ssp_spelled,
          'data/0_taxonomy/IOC_Names_File_Plus-15.1_wide.csv')
metaGenerator('../data/taxonomy/taxonomy-IOC/')
