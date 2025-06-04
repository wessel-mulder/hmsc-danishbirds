rm(list = ls())
# GETTING STARTED ---------------------------------------------------------
library('terra')
library('dplyr')

# LOADING OCCURRENCES  ----------------------------------------------------
data <- read.csv('../data/distributions/dof_atlas/atlasdata_DOF_27-01-2025.csv',
                 sep = ';')

data <- data[data$sandsynlighed %in% c('sikker','sandsynlig'),] # remove the mulig records 
data

# survey dataframe
survey <- read.csv('data/1_preprocessing/design/studyDesign.csv')

# EXTRACT CODES AND RELEVANT INFO ---------------------------------------------------------
translate_letters <- function(code) {
  # get first two letters 
  first <- unlist(strsplit(substr(as.character(code),1,2),""))
  letter_one <- toupper(letters[as.integer(first[1])])
  letter_two <- toupper(letters[as.integer(first[2])])
  
  # get the numbers 
  last <- substr(as.character(code),3,4)
  name <- paste0(letter_one,letter_two,last)
  return(name)
}

data$site <- vapply(data$kvadratnr,translate_letters,character(1))
data$survey <- paste0(data$site,'_',data$atlas)
data$latin_nospace_DOF <- gsub(' ','_',data$latin)
data$presence <- 1

data <- data[,c('survey','latin_nospace_DOF','presence')]

# call new data
df <- data.frame(matrix(nrow = nrow(survey),
                        ncol = 0))
df$survey <- survey$survey

# CALL NEW DATAFRAME ------------------------------------------------------

for(species in unique(data$latin_nospace_DOF)){
  print(species)
  
  if(species == unique(data$latin_nospace_DOF)[1]){
    df2 <- df
  }
  
  df2 <- left_join(df2,data[data$latin_nospace_DOF == species,c('survey','presence')],by='survey')
  colnames(df2)[which(names(df2) == "presence")] <- species
  
  if (species == unique(data$latin_nospace_DOF)[length(unique(data$latin_nospace_DOF))]) {
    df2[is.na(df2)] <- 0
  }
}

# LOADING TRAITS 
traits <- read.csv('data/1_preprocessing/Tr_aits/traits-guild_migration.csv')
species_list <- traits$latin_DOF_underscores #199 individuals 

# SELECT ONLY SPECIES IN THE SPECIES LIST 
occurrences_subset <- df2 %>%
  select(all_of(c('survey',species_list))) 

# 200 colums  correct, because 1 one of them is the survey 

head(occurrences_subset,5)
write.csv(occurrences_subset,
          'data/1_preprocessing/Y_occurrences/Y_occurrences.csv',
          row.names = F)




