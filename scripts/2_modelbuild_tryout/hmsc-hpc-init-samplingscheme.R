rm(list = ls())

# Define MCMC settings
nChains <- 4
thin <- 10
nSamples <- 2000
transient <- thin*nSamples
verbose <- 100
params <- list(
  nChains = nChains,
  thin = thin,
  nSamples = nSamples,
  transient = transient,
  verbose = verbose
)


# GETTING STARTED ---------------------------------------------------------
if (interactive() && Sys.getenv("RSTUDIO") == "1") {
  message("Running in RStudio")
  library(Hmsc)
  library(jsonify)
  library(knitr)
  library(corrplot)
  library(ape)
  library(dplyr)
  input <- '.'
  python <- file.path(getwd(), 'hmsc-hpc-main',"hmsc-venv", "bin", "python3.11")
  flagInit = 0
  flagFitR = 1
} else {
  message("Running from terminal or non-interactive environment")
  library(RColorBrewer,lib="~/Rlibs")
  library(farver,lib="~/Rlibs")
  library(scales,lib="~/Rlibs")
  library(jsonify,lib="~/Rlibs")
  library(ape,lib="~/Rlibs")
  library(dplyr,lib="~/Rlibs")
  library(Hmsc,lib="~/Rlibs")
  input <- '~/home/projects/hmsc-danishbirds'
  python <- '~/opt/software/anaconda3/2021.05/bin/python'
  flagInit = 1
  flagFitR = 0
}

# check if accurately installed
summary(TD)

# LOADING DATA -----------------------------------------------------

# --> LOAD ENVIRONMENT
X <- read.csv(file.path(input,'data/1_preprocessing/X_environmental/X_Environmental.csv'),row.names=1)
X <- X[sort(row.names(X)),]

# keep only sites with data  
X <- na.omit(X)
sites_actual <- row.names(X)

# --> LOAD OCCURRENCES 
Y <- read.csv(file.path(input,'data/1_preprocessing/Y_occurrences/Y_occurrences.csv'),row.names=1)
# merge subspecies for now 
# Merge Motacilla_alba_yarrellii into Motacilla_alba
Y$Motacilla_alba <- pmax(Y$Motacilla_alba, Y$Motacilla_alba_yarrellii, na.rm = TRUE)
Y$Motacilla_alba_yarrellii <- NULL  # Remove the subspecies column

# Merge Motacilla_flava_flavissima into Motacilla_flava
Y$Motacilla_flava <- pmax(Y$Motacilla_flava, Y$Motacilla_flava_flavissima, na.rm = TRUE)
Y$Motacilla_flava_flavissima <- NULL  # Remove the subspecies column

# rename acanthis to species name
names(Y)[names(Y) == "Acanthis_flammea_cabaret"] <- "Acanthis_flammea"
Y <- Y[sort(row.names(Y)), sort(colnames(Y))]

# remove sites without data 
Y <- Y[row.names(Y) %in% sites_actual,]

# --> LOAD TRAITS 
Tr <- read.csv(file.path(input,"data/1_preprocessing/Tr_aits/traits-guild_migration.csv"),row.names = 2)[,c(2,3)]
Tr <- Tr[sort(row.names(Tr)),]
Tr['Motacilla_alba_yarrellii',] == Tr['Motacilla_alba',] # TRUE 
Tr['Motacilla_flava_flavissima',] == Tr['Motacilla_flava',] # MIGRATION FALSE 
Tr <- Tr[!(rownames(Tr) %in% c('Motacilla_alba_yarrellii','Motacilla_flava_flavissima')), ]
rownames(Tr)[rownames(Tr) == "Acanthis_flammea_cabaret"] <- "Acanthis_flammea"

# --> LOAD STUDY DESIGN 
Design <- read.csv(file.path(input,"data/1_preprocessing/design/studyDesign.csv"),row.names=5)

# remove sites without data 
Design <- Design[row.names(Design) %in% sites_actual,]

# sort
Design <- Design[sort(row.names(Design)), ]

# convert to factors
Design$site <- as.factor(Design$site)
Design$atlas <- as.factor(Design$atlas)
Design$year[Design$atlas == '1'] <- 1971
Design$year[Design$atlas == '2'] <- 1992
Design$year[Design$atlas == '3'] <- 2014


# phylogeny
phy <- read.tree(file.path(input,'data/1_preprocessing/Taxonomy/tree_fromPD.tre'))
# check pd
pd_matrix <- cophenetic.phylo(phy)
pd_matrix <- pd_matrix[sort(rownames(pd_matrix)), sort(colnames(pd_matrix))]

# check if species-lists are identical
table(rownames(pd_matrix) == names(Y))
# stunning 

# get xycoords
xycoords <- Design[,colnames(Design) %in% c('lat','lon','site')]
xycoords <- distinct(xycoords)
rownames(xycoords) <- xycoords$site
xycoords <- xycoords[,colnames(xycoords) %in% c('lat','lon')]

# PREPARING MODEL BUILD ---------------------------------------------------
# Define model types: 
mtVec = c(1:4) # 1-nngp, 2-gpp, 3-enviro, 4-full
mtSuffix = c('nng','gpp','env','ful')

date <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
for(mt in mtVec){
  ### SAME ACROSS ALL MODELS
  # Define model formulas for environmental and trait data
  XFormula <- as.formula(paste("~", paste(colnames(X), collapse = "+"), sep = " "))
  TrFormula <- as.formula(paste("~", paste(colnames(Tr), collapse = "+"), sep = " "))
  
  Design$year <- as.numeric(as.character(Design$year))
  time <- data.frame(year = sort(unique(Design$year)))
  rownames(time) <- sort(unique(Design$year))
  struc_time <- HmscRandomLevel(sData = time)
  Design$year <- as.factor(Design$year)
  
  if(mt==1){
    xycoords <- as.matrix(xycoords)
    struc_space <- HmscRandomLevel(sData = xycoords, sMethod = "NNGP", nNeighbours = 10)
  }else if(mt==2){
    xycoords <- as.matrix(xycoords)
    xyKnots = constructKnots(xycoords,knotDist = 0.5, minKnotDist = 0.3)
    plot(xycoords[,2],xycoords[,1],pch=18, asp=1)
    points(xyKnots[,2],xyKnots[,1],col='red',pch=18)
    struc_space <- HmscRandomLevel(sData = xycoords, sMethod = "GPP",
                                   sKnot = xyKnots)
  }else if(mt==3){
    struc_space <- HmscRandomLevel(sData = xycoords, sMethod = "Full")
    ## EDIT FORMULA FOR LESS ENVIRONMENTAL PARAMS 
    Xsub <- X %>% 
      select(tmean_year,tmean_winter,tmean_breeding,
             prec_year,prec_winter,prec_breeding,
             hh,unique)
    XFormula <- as.formula(paste("~", paste(colnames(Xsub), collapse = "+"), sep = " "))
    
  }else if(mt==4){
    struc_space <- HmscRandomLevel(sData = xycoords, sMethod = "Full")
  }
  # define m 
  m <-Hmsc(Y = Y, XData = X, XFormula = XFormula, TrData = Tr,
           TrFormula = TrFormula, phyloTree = phy,
           studyDesign = Design[,c(1,5)], 
           ranLevels = list('site'=struc_space,'year'=struc_time),
           distr='probit')
  
  ### IN RSTUDIO START SAMPLING 
  if(flagFitR){
    print('Rstudio stuff executed')
  }
  ### IN HPC ENVIORNMENT SET UP INIT
  if(flagInit){
  # initiate mcmc sampling 
  init_obj <- sampleMcmc(m, samples=nSamples, thin=thin,
                         transient=transient, nChains=nChains,
                         verbose = verbose,
                         engine="HPC")
  
  dir_name <- paste0(date,'_',mtSuffix[mt])
  dir.create(file.path(input,'tmp_rds',dir_name))
  
  init_file_path = file.path(input,'tmp_rds',dir_name, "init_file.rds")
  m_file_path = file.path(input,'tmp_rds',dir_name, "m_object.rds")
  param_file_path = file.path(input,'tmp_rds',dir_name, "params.rds")
  lines <- paste(names(params), unlist(params), sep = "=")
  writeLines(lines, file.path(input,'tmp_rds',dir_name, "params.txt"))
  
  
  # save as json 
  saveRDS(to_json(init_obj), file=init_file_path)
  saveRDS(m,file=m_file_path)
  saveRDS(params,file=param_file_path)
  
  
  # operates in python, so formulate the required call 
  post_file_path = file.path(input,'tmp_rds',dir_name, "post_file.rds")
  python_cmd_args = paste("-m hmsc.run_gibbs_sampler",
                          "--input", shQuote(init_file_path),
                          "--output", shQuote(post_file_path),
                          "--samples", nSamples,
                          "--transient", transient,
                          "--thin", thin,
                          "--verbose", verbose)
  cat(paste(shQuote(python), python_cmd_args), "\n")
  print('Init files created')
  }
}








