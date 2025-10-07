rm(list = ls())

# Define MCMC settings
env_vars <- c('tmean_year')
spatial_method <- c('GPP','Full')

nChains <- 4
thin <- 10
nSamples <- 250
transient <- 50000
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
  library(dplyr,lib="~/Rlibs")
  library(withr,lib="~/Rlibs")
  input <- '~/home/projects/hmsc-danishbirds'
  python <- '~/opt/software/anaconda3/2021.05/bin/python'
  flagInit = 1
  flagFitR = 0
}

# check if accurately installed
summary(TD)

# LOADING DATA -----------------------------------------------------

### ENVIRONMENT
X <- read.csv(file.path(input,'data/1_preprocessing/X_environmental/X_Environmental.csv'),row.names=1)
X <- X[sort(row.names(X)),]

# keep only sites with data  
X <- na.omit(X)

# keep only atlas 3
atlas3 <- X[rownames(X)[grep("_3$", rownames(X))],,drop=F]

# grab tmean_year

sites_actual <- row.names(atlas3)

### OCCURRENCES 
Y <- read.csv(file.path(input,'data/1_preprocessing/Y_occurrences/Y_occurrences.csv'),row.names=1)

# remove sites without data 
Y <- Y[row.names(Y) %in% sites_actual,]

# grab 13 warblers
genera <- c('Phylloscopus','Curruca','Sylvia','Acrocephalus','Hippolais')
keep <- sapply(strsplit(colnames(Y),'_'),head,1) %in% genera
Y_warblers <- Y[,keep]

# check if they have enough records 
colSums(Y_warblers, na.rm =T)
# barred warbler is absent so we remove it 
Y_warblers <- Y_warblers[colnames(Y_warblers) != "Curruca_nisoria"]
# fixed
colSums(Y_warblers, na.rm =T)


# --> LOAD STUDY DESIGN 
Design <- read.csv(file.path(input,"data/1_preprocessing/design/studyDesign.csv"),row.names=5)

# remove sites without data 
Design <- Design[row.names(Design) %in% sites_actual,]

# sort
Design <- Design[sort(row.names(Design)), ]
Design3 <- Design[rownames(Design)[grep("_3$", rownames(Design))],,drop=F]
Design3 <- Design3[,c('site','lat','lon')]

head(Design3)

# convert to factors
Design3$site <- as.factor(Design3$site)

xycoords <- data.frame(lon = Design3$lon, lat = Design3$lat)
rownames(xycoords) <- Design3$site

v <- vect(xycoords, geom = c("lon","lat"), crs = "EPSG:4326")
v_proj <- project(v, "EPSG:23032")

proj_xycoords <- crds(v_proj)
rownames(proj_xycoords) <- rownames(xycoords)
head(proj_xycoords)


# LOAD KNOTS  -------------------------------------------------------------
knots <- vect(file.path(input,"data/1_preprocessing/knots_all/knots.shp"))
knots_proj <- project(knots,"EPSG:23032")
plot(knots_proj,col='green')
plot(v_proj,add=T,pch=19,cex=0.5)

# Get centroids
centroids <- centroids(knots_proj)
# get coords
knots_xycoords <- crds(centroids)

#plot
plot(knots_xycoords,col='red',pch=19,cex=1)
plot(v_proj,cex=0.5,add=T)


# PREPARING MODEL BUILD ---------------------------------------------------
# Define model types: 

date <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
for(spatial_p in spatial_method){
for(env in env_vars){
  ### SAME ACROSS ALL MODELS
  # Define model formulas for environmental and trait data
  X <- atlas3[,env,drop=F]
  XFormula <- as.formula(paste("~", paste(colnames(X), collapse = "+"), sep = " "))
  }
  
  if(spatial_p == 'GPP'){
  struc_space <- HmscRandomLevel(sData = proj_xycoords, sMethod = "GPP",
                                 sKnot = knots_xycoords)
  }else if(spatial_p == 'Full'){
  struc_space <- HmscRandomLevel(sData = proj_xycoords, sMethod = "Full",
                                 sKnot = knots_xycoords)
  }
  # set distances to min and max distances between sites 
  freq <- c(0.5,rep(0.005,100))
  samples <- c(0,seq(from = 4999, to = 477312, length.out = 100))
  small <- cbind(samples,freq)
  
  struc_space_small <- setPriors(struc_space,alphapw=small)

m <-Hmsc(Y = Y_warblers, 
         XData = X, 
         XFormula = XFormula,
         studyDesign = Design3[,c('site'),drop=F], 
         ranLevels = list('site'=struc_space_small),
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
  
  dir_name <- paste0(date,'_singleev_',env,'_spatialmethod_','spatialp')
  print(dir_name)
  dir.create(file.path(input,'tmp_rds',dir_name))
  
  init_file_path = file.path(input,'tmp_rds',dir_name, "init_file.rds")
  m_file_path = file.path(input,'tmp_rds',dir_name, "m_object.rds")
  param_file_path = file.path(input,'tmp_rds',dir_name, "params.rds")
  lines <- paste(names(params), unlist(params), sep = "=")
  writeLines(lines, file.path(input,'tmp_rds',dir_name, "params.txt"))
  if(spatial_p == 'GPP'){
  pdf(file.path(input,'tmp_rds',dir_name,'knots.pdf'))
  plot(knots_xycoords,col='red',pch=19,cex=1)
  plot(v_proj,cex=0.5,add=T)
  dev.off()
  }
  
  
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

