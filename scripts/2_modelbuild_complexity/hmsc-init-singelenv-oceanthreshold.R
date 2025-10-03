rm(list = ls())

# Define MCMC settings
env_vars <- c('tmean_year')

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
  library(sp)
  library(terra)
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
  library(sp,lib="~/Rlibs")
  library(terra,lib="~/Rlibs")

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
genera <- c('Phylloscopus','Curruca','Sylvia','Acrocephalus')
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

# PREPARING MODEL BUILD ---------------------------------------------------
# Define model types: 

date <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")

for(knotdist in knot_dists)
for(env in env_vars){
  ### SAME ACROSS ALL MODELS
  # Define model formulas for environmental and trait data
  if(env == 'NULL'){
      XFormula <- ~1
        X <- atlas3[,'tmean_winter',drop=F]

  }else{
    X <- atlas3[,env,drop=F]
  XFormula <- as.formula(paste("~", paste(colnames(X), collapse = "+"), sep = " "))
  }

  xyKnots = constructKnots(proj_xycoords,knotDist = knotdist,minKnotDist =  knotdist)
  nKnots <- nrow(xyKnots)
  params[["nKnots"]] <- nKnots
  params[["knotdist"]] <- knotdist
  
  struc_space <- HmscRandomLevel(sData = proj_xycoords, sMethod = "GPP",
                                 sKnot = xyKnots)
                                 
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
  
  dir_name <- paste0(date,'_singleev_',env,'_knots_',nKnots)
  print(dir_name)
  dir.create(file.path(input,'tmp_rds',dir_name))
  
  init_file_path = file.path(input,'tmp_rds',dir_name, "init_file.rds")
  m_file_path = file.path(input,'tmp_rds',dir_name, "m_object.rds")
  param_file_path = file.path(input,'tmp_rds',dir_name, "params.rds")
  lines <- paste(names(params), unlist(params), sep = "=")
  writeLines(lines, file.path(input,'tmp_rds',dir_name, "params.txt"))
  pdf(file.path(input,'tmp_rds',dir_name,'knots.pdf'))
    plot(proj_xycoords[,1],proj_xycoords[,2],pch=18, asp=1,
       main=paste0('Number of knots: ',nKnots))
  points(xyKnots[,1],xyKnots[,2],col='red',pch=18)
  print(nKnots)
  dev.off()
  
  
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