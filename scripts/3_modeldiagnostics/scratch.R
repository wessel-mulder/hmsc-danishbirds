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

# GETTING IDEA OF SPECIES TO LOOK AT --------------------------------------
library(dplyr)
library(ggplot2)
dir <- file.path('.','tmp_rds','mods-complexity-v2')
dirs <- list.dirs(dir,
                  recursive = F)
dirs <- dirs[grepl('.*_(1|2|3)$',dirs)]
dir<-dirs[3]
outputs <- lapply(dirs,function(dir){
  # read objects 
  MF <- readRDS(file.path(dir,'model-outputs','model-fit.rds'))
  VP <- readRDS(file.path(dir,'model-outputs','VP.rds'))
  m <- readRDS(file.path(dir,'m_object.rds'))
  X <- m$XData
  atlas <- unique(sapply(strsplit(rownames(X),'_'),'[',2))
  
  #transpose VP
  df <-VP$vals # get vals
  df <- as.data.frame(t(df)) # transpose into dataframe 
  
  # merge VP with species info 
  Tjur<- data.frame(tjur = MF$TjurR2)
  df_rel <- Tjur$tjur*df

  spNames <- m$spNames
  row.names(Tjur) <- spNames
  Tr <- m$TrData
  species_fit_char <- merge_by_rownames(Tr,Tjur)
  species_fit_char <- merge_by_rownames(species_fit_char,df_rel)
  species_fit_char$atlas <- atlas
  species_fit_char
})
merged <- do.call(rbind,outputs)
str(merged)
colnames(merged) <- c('Migration','Guild','TjurR2',
                      'Temperature','Precipitation','Landscape heterogeneity',
                      "Random effect - site", 'Atlas')

# count species per guild
guild_counts <- merged %>%
  group_by(Guild) %>%
  summarise(n_species = n()/3) %>%
  mutate(Guild_counts = paste0(Guild,' (n=',n_species,')')) %>%
  ungroup()

str(merged)
merged$Guild_counts <- merged$Guild
merged$Guild_counts <- guild_counts$Guild_counts[match(merged$Guild, guild_counts$Guild)]

vars2inspect <- c('Temperature','Precipitation','Landscape heterogeneity')
lapply(vars2inspect,function(var2inspect){
v2inspect <- sym(var2inspect)

guild_means <- merged %>%
  group_by(Guild_counts, Atlas) %>%
  summarise(!!v2inspect := mean(!!v2inspect, na.rm = TRUE), .groups = "drop")
colors <- c('darkblue','darkgreen','darkred')

ggplot(merged,
       aes(x = reorder(Guild_counts, !!v2inspect, FUN = mean, na.rm = TRUE),
           y=!!v2inspect,
           col = Atlas, fill = Atlas))+
  # violins for guilds with >1 species
  #geom_violin(drop = FALSE, alpha = 0.6) +
  geom_point(alpha = 0.5,size=0.5)+

  # Points for each species, colored by Atlas
  # Points for mean Tjur² per guild per Atlas
  geom_point(data = guild_means,
             aes(y = !!v2inspect, fill = Atlas,color=Atlas),
             size = 1, shape = 21, stroke = 1) +
  
  # set colors 
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  
  #scale axes 
  scale_y_continuous(limits=(c(0,0.175)))+
  
  # adjustements 
  labs(title=var2inspect)+
  xlab(label='Guilds (number of species)')+
  ylab(label='Relative variable importance')+
  
  coord_flip()+
  theme_minimal()+
  theme(
    axis.line = element_blank(),             # remove axis lines
    #panel.grid.major = element_blank(),      # remove major gridlines
    #panel.grid.minor = element_blank()       # remove minor gridlines
  ) 
})

# PLOTTING ENVIRONMENTAL VARS  --------------------------------------------
### load model with all 3 atlases 

library(ggplot2)
library(ggExtra)
library(cowplot)
library(gridExtra)
library(grid)

m <- readRDS('tmp_rds/mods-complexity-v2/2025-10-24_17-32-54_threeenv_allspecies_atlas_123/m_object.rds')
X <- m$XData
Y <- m$Y
design <- m$studyDesign

merge <- merge(design,X,by='row.names')
rownames(merge) <- merge$Row.names
merge <- merge[,!colnames(merge) %in% c('Row.names')]
merge <- merge(merge,Y,by='row.names')
rownames(merge) <- merge$Row.names
merge <- merge[,!colnames(merge) %in% c('Row.names')]

species <- 'Caprimulgus_europaeus'
species <- 'Astur_gentilis'
species <- 'Haematopus_ostralegus'
species_eng <- 'Oystercatcher'
subs <- c('atlas','tmean_year','prec_year','hh',species)
subs <- merge[,colnames(merge) %in% subs]

colnames(subs) <- c('Atlas',
                    'Temperature','Precipitation','Landscape heterogeneity',
                    'PresenceAbsence')
occs_a1 <- sum(subs$PresenceAbsence[subs$Atlas==1])
occs_a2 <- sum(subs$PresenceAbsence[subs$Atlas==2])
occs_a3 <- sum(subs$PresenceAbsence[subs$Atlas==3])
subs$PresenceAbsence <- as.character(subs$PresenceAbsence)

combinations <- list(
  c('Temperature','Precipitation'),
  c('Temperature','Landscape heterogeneity'),
  c('Precipitation','Landscape heterogeneity')
)
combination <- combinations[[2]]
# make some plots 
lapply(combinations,function(combination){
  grid.newpage()
  print(combination)
  v1 <- sym(combination[1])
  v2 <- sym(combination[2])

cols <- c('darkred','darkgreen','darkblue')
alphas <- c(0.1,1)
pmain <- ggplot(data=subs,
       aes(x=!!v1,y=!!v2,col=Atlas,alpha=PresenceAbsence))+
  geom_point()+
  scale_color_manual(values=cols)+
  scale_alpha_manual(values=alphas)+
  labs(title = abbreviate_genus(species))+
  theme_minimal()
### VERSION WITH COLORS 
# Marginal densities along x axis
xdens <- axis_canvas(pmain, axis = "x") +
  geom_density(data = subs, 
               aes(x = !!v1, 
                   #fill = Atlas,
                   col=Atlas,
                   linetype = PresenceAbsence),
               size = 0.2,
               linewidth = 0.5)+
  #scale_fill_manual(values=cols)+
  scale_color_manual(values=cols)+
  scale_linetype_manual(values=c(3,1))
# Marginal densities along x axis
ydens <- axis_canvas(pmain, axis = "y",coord_flip=T) +
  geom_density(data = subs, 
               aes(x = !!v2, 
                   #fill = Atlas,
                   col=Atlas,
                   linetype = PresenceAbsence),
               size = 0.2,
               alpha=1,
               linewidth=0.5)+
  #scale_fill_manual(values=cols)+
  scale_color_manual(values=cols)+
  scale_linetype_manual(values=c(3,1))+
coord_flip()

### GET TOGETHER 
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
ggdraw(p2)

### ANNOTATE 
# grid.text(label='Presences in:',x=0.87,y=0.725,just='left',
#           gp=gpar(col='black'))+
# grid.text(label=paste0('Atlas 1 - ',occs_a1),x=0.87,y=0.7,just='left',
#           gp=gpar(col='darkred'))+
# grid.text(label=paste0('Atlas 2 - ',occs_a2),x=0.87,y=0.675,just='left',
#           gp=gpar(col='darkgreen'))+
# grid.text(label=paste0('Atlas 3 - ',occs_a3),x=0.87,y=0.65,just='left',
#           gp=gpar(col='darkblue'))


})


ps

p+ggMarginal(x=c(8,9,10),
           y=c(600,700,900),
           'density')


# COMPARE MODEL WITH 1/2/3 ATLAS  ---------------------------------------------
dir <- file.path('.','tmp_rds','mods-complexity-v2')
dirs <- list.dirs(dir,
                  recursive = F)
dirs <- dirs[grepl('.*_(1|2|3)$',dirs)]
dir<-dirs[1]
outputs <- lapply(dirs,function(dir){
m <- readRDS(file.path(dir,'m_object.rds'))
X <- m$XData
Y <- m$Y
design <- m$studyDesign
design$atlas <- sapply(strsplit(rownames(X),'_'),'[',2)

merge <- merge(design,X,by='row.names')
rownames(merge) <- merge$Row.names
merge <- merge[,!colnames(merge) %in% c('Row.names')]
merge <- merge(merge,Y,by='row.names')
rownames(merge) <- merge$Row.names
merge <- merge[,!colnames(merge) %in% c('Row.names')]


species <- 'Picus_viridis'
species_eng <- 'Green woodpecker'
subs <- c('atlas','tmean_year','prec_year','hh',species)
subs <- merge[,colnames(merge) %in% subs]
colnames(subs) <- c('Atlas',
                    'Temperature','Precipitation','Habitat heterogeneity',
                    'PresenceAbsence')
subs
})
subs <- do.call(rbind,outputs)

occs_a1 <- sum(subs$PresenceAbsence[subs$Atlas==1])
occs_a2 <- sum(subs$PresenceAbsence[subs$Atlas==2])
occs_a3 <- sum(subs$PresenceAbsence[subs$Atlas==3])
subs$PresenceAbsence <- as.character(subs$PresenceAbsence)
head(subs)
ggplot(data=subs,
       aes(x=Temperature,y=Precipitation,col=Atlas,alpha=PresenceAbsence))+
  geom_point()+
  scale_color_manual(values=c('darkred','darkgreen','darkblue'))+
  scale_alpha_manual(values=c(0.1,1))+
  annotate('text',
           label=paste0('Number of presences:'),
           x = 9.78,y=920,
           cex=5,col='black')+
  annotate('text',
           label=paste0('Atlas 1: ',occs_a1),
           x = 9.78,y=890,
           cex=5,col='darkred')+
  annotate('text',
           label=paste0('Atlas 2: ',occs_a2),
           x = 9.78,y=870,
           cex=5,col='darkgreen')+
  annotate('text',
           label=paste0('Atlas 3: ',occs_a3),
           x = 9.78,y=850,
           cex=5,col='darkblue')+
  labs(title = paste0(abbreviate_genus(species),' - ',species_eng))+
  theme_minimal()

# PLOT DISTRIBUTIONS WITH LAND USE CATEGORIES -------------------------------------------------------------------------
# GETTING STARTED 
input <- '.'

### ENVIRONMENT
X <- read.csv(file.path(input,'data/1_preprocessing/X_environmental/X_Environmental.csv'),row.names=1)
X <- X[sort(row.names(X)),]

# get ocean thresholds
grids_thresholds <- st_read(file.path(input,'data/1_preprocessing/atlas-grids/grids-ocean-thresholds/grids_ocean_thresholds.shp'))
thresholds <- grids_thresholds$kvdrtkd[grids_thresholds$pct_lnd>=25]

X <- X[sub("_[123]$", "", rownames(X)) %in% thresholds,]
# any NAs?
table(is.na(X))
if(!env_vars=='all'){X <- X[,env_vars,drop=F]}
table(is.na(X)) # still som NAs
X <- na.omit(X)
table(is.na(X)) # fixed

# grab tmean_year
sites_actual <- row.names(X)

### OCCURRENCES 
Y <- read.csv(file.path(input,'data/1_preprocessing/Y_occurrences/Y_occurrences.csv'),row.names=1)

# remove sites without data 
Y <- Y[row.names(Y) %in% sites_actual,]
X_sub <- X[grep('LULC_*',colnames(X))]
XY <- merge_by_rownames(X_sub,Y)
head(XY)
species <- 'Haematopus_ostralegus'
species_eng <- 'Oystercatcher'
subs <- c('atlas','tmean_year','prec_year','hh',species)
subs <- merge[,colnames(merge) %in% subs]

colnames(subs) <- c('Atlas',
                    'Temperature','Precipitation','Landscape heterogeneity',
                    'PresenceAbsence')
occs_a1 <- sum(subs$PresenceAbsence[subs$Atlas==1])
occs_a2 <- sum(subs$PresenceAbsence[subs$Atlas==2])
occs_a3 <- sum(subs$PresenceAbsence[subs$Atlas==3])
subs$PresenceAbsence <- as.character(subs$PresenceAbsence)

combinations <- list(
  c('Temperature','Precipitation'),
  c('Temperature','Landscape heterogeneity'),
  c('Precipitation','Landscape heterogeneity')
)
combination <- combinations[[2]]
# make some plots 
lapply(combinations,function(combination){
  grid.newpage()
  print(combination)
  v1 <- sym(combination[1])
  v2 <- sym(combination[2])
  
  cols <- c('darkred','darkgreen','darkblue')
  alphas <- c(0.1,1)
  pmain <- ggplot(data=subs,
                  aes(x=!!v1,y=!!v2,col=Atlas,alpha=PresenceAbsence))+
    geom_point()+
    scale_color_manual(values=cols)+
    scale_alpha_manual(values=alphas)+
    labs(title = abbreviate_genus(species))+
    theme_minimal()
  ### VERSION WITH COLORS 
  # Marginal densities along x axis
  xdens <- axis_canvas(pmain, axis = "x") +
    geom_density(data = subs, 
                 aes(x = !!v1, 
                     #fill = Atlas,
                     col=Atlas,
                     linetype = PresenceAbsence),
                 size = 0.2,
                 linewidth = 0.5)+
    #scale_fill_manual(values=cols)+
    scale_color_manual(values=cols)+
    scale_linetype_manual(values=c(3,1))
  # Marginal densities along x axis
  ydens <- axis_canvas(pmain, axis = "y",coord_flip=T) +
    geom_density(data = subs, 
                 aes(x = !!v2, 
                     #fill = Atlas,
                     col=Atlas,
                     linetype = PresenceAbsence),
                 size = 0.2,
                 alpha=1,
                 linewidth=0.5)+
    #scale_fill_manual(values=cols)+
    scale_color_manual(values=cols)+
    scale_linetype_manual(values=c(3,1))+
    coord_flip()
  
  ### GET TOGETHER 
  p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
  p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
  ggdraw(p2)
  
  ### ANNOTATE 
  # grid.text(label='Presences in:',x=0.87,y=0.725,just='left',
  #           gp=gpar(col='black'))+
  # grid.text(label=paste0('Atlas 1 - ',occs_a1),x=0.87,y=0.7,just='left',
  #           gp=gpar(col='darkred'))+
  # grid.text(label=paste0('Atlas 2 - ',occs_a2),x=0.87,y=0.675,just='left',
  #           gp=gpar(col='darkgreen'))+
  # grid.text(label=paste0('Atlas 3 - ',occs_a3),x=0.87,y=0.65,just='left',
  #           gp=gpar(col='darkblue'))
  
  
})


# HELPER FUNCTIONS  -------------------------------------------------------
abbreviate_genus <- function(x) {
  sapply(x, function(nm) {
    parts <- unlist(strsplit(nm, "_"))  # split by underscore
    paste0(substr(parts[1], 1, 1), ". ", parts[2])
  })
}

merge_by_rownames <- function(x,y) {
  merge <- merge(x,y,by='row.names')
  rownames(merge) <- merge$Row.names
  merge[,!colnames(merge) %in% c('Row.names')]
}
