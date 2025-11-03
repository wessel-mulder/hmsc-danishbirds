rm(list = ls())


# SANITY CHECK ALLSPECIES ATLAS1 ------------------------------------------
#### GETTING STARTED 
library(Hmsc)
library(jsonify)
library(dplyr)
library(tibble)
library(ggplot2)
library(fields)
#### LOADING MODEL
between <- 'mods-complexity-v2'
dirs <- list.dirs(file.path('./tmp_rds',between),recursive=F)
input <- dirs[grepl('threeenv_allspecies_atlas_1',dirs)]
  
# load model j
m <- readRDS(file.path(input,'m_object.rds'))
# load params 
params <- readRDS(file.path(input,'params.rds'))
nChains <- params$nChains
nSamples <- params$nSamples
thin <- params$thin
transient <- params$transient

chainList = vector("list", nChains)
for(cInd in 1:nChains){
  chain_file_path = file.path(input, sprintf("post_chain%.2d_file.rds", cInd-1))
  #print(chain_file_path)
  if(file.exists(chain_file_path)) {
    chainList[[cInd]] = from_json(readRDS(file = chain_file_path)[[1]])[[1]]
  }
}

filteredList <- chainList
fitSepTF = importPosteriorFromHPC(m, filteredList, nSamples, thin, transient)
mpost <-convertToCodaObject(fitSepTF)

### LOAD ENVIRONMENT
XData <- fitSepTF$XData
Y <- as.data.frame(fitSepTF$Y)
studyDesign <- fitSepTF$studyDesign
coords <- fitSepTF$ranLevels$site$s
merged <- XData %>%
  rownames_to_column('id') %>%
  left_join(studyDesign %>% rownames_to_column("id"), join_by(id)) %>%
  left_join(coords %>% rownames_to_column("site"), join_by(site)) %>%
  column_to_rownames('id')
# plot 
ggplot(merged,
       aes(x=X,y=Y,col=hh))+
  geom_point()

merged_occ <- Y %>%
  rownames_to_column('id') %>%
  left_join(studyDesign %>% rownames_to_column("id"), join_by(id)) %>%
  left_join(coords %>% rownames_to_column("site"), join_by(site)) %>%
  column_to_rownames('id')

ggplot(merged_occ,
       aes(x=X,y=Y,col=Picus_viridis))+
  geom_point()



# MAKING PREDICTIONS IN ANOTHER ATLAS ------------------------------------------
#### GETTING STARTED 
library(Hmsc)
library(jsonify)
library(dplyr)
library(tibble)
library(ggplot2)
library(fields)
#### LOADING MODEL
between <- 'mods-complexity-v2'
dirs <- list.dirs(file.path('./tmp_rds',between),recursive=F)
input <- dirs[grepl('threeenv_allspecies_atlas_1',dirs)]

# load model j
m <- readRDS(file.path(input,'m_object.rds'))
# load params 
params <- readRDS(file.path(input,'params.rds'))
nChains <- params$nChains
nSamples <- params$nSamples
thin <- params$thin
transient <- params$transient

chainList = vector("list", nChains)
for(cInd in 1:nChains){
  chain_file_path = file.path(input, sprintf("post_chain%.2d_file.rds", cInd-1))
  #print(chain_file_path)
  if(file.exists(chain_file_path)) {
    chainList[[cInd]] = from_json(readRDS(file = chain_file_path)[[1]])[[1]]
  }
}

filteredList <- chainList
fitSepTF = importPosteriorFromHPC(m, filteredList, nSamples, thin, transient)
mpost <-convertToCodaObject(fitSepTF)

### figure out which atlases are represented 
XData 