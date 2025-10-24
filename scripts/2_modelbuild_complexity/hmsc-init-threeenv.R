rm(list = ls())

# Define MCMC settings
env_vars <- c('tmean_year','prec_year','hh')
chars <- c('all')
atlases <- list(c('1','2'),
                c('1','3'),
                c('2','3'),
                c('1','2','3')
)

nChains <- 4
thin <- 10
nSamples <- 250
transient <- 100000
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
  library(sf)
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
  library(sf,lib="~/Rlibs")
  
  input <- '~/home/projects/hmsc-danishbirds'
  python <- '/maps/projects/cmec/people/bhr597/projects/hmsc-danishbirds/hmsc-venv'
  flagInit = 1
  flagFitR = 0
}

# check if accurately installed
summary(TD)

# LOADING DATA -----------------------------------------------------

### ENVIRONMENT
X <- read.csv(file.path(input,'data/1_preprocessing/X_environmental/X_Environmental.csv'),row.names=1)
X <- X[sort(row.names(X)),]

# get ocean thresholds
grids_thresholds <- st_read(file.path(input,'data/1_preprocessing/atlas-grids/grids-ocean-thresholds/grids_ocean_thresholds.shp'))
thresholds <- grids_thresholds$kvdrtkd[grids_thresholds$pct_lnd>=25]

X <- X[sub("_[123]$", "", rownames(X)) %in% thresholds,]
# any NAs?
table(is.na(X))
X <- X[,env_vars,drop=F]
table(is.na(X)) # still som NAs
X <- na.omit(X)
table(is.na(X)) # fixed

# grab tmean_year
sites_actual <- row.names(X)

### OCCURRENCES 
Y <- read.csv(file.path(input,'data/1_preprocessing/Y_occurrences/Y_occurrences.csv'),row.names=1)

# remove sites without data 
Y <- Y[row.names(Y) %in% sites_actual,]

# grab 13 warblers
#genera <- c('Phylloscopus','Curruca','Sylvia','Acrocephalus','Hippolais','Locustella')
#keep <- sapply(strsplit(colnames(Y),'_'),head,1) %in% genera
#Y_warblers <- Y[,keep]

# barred warbler is absent so we remove it 
#Y_warblers <- Y_warblers[colnames(Y_warblers) != "Curruca_nisoria"]
#Y_warblers <- Y_warblers[colnames(Y_warblers) != "Locustella_fluviatilis"]
#Y_warblers <- Y_warblers[colnames(Y_warblers) != "Phylloscopus_trochiloides"]

# phylloscopus also very absent in some of the thresholds now 
#Y_warblers <- Y_warblers[colnames(Y_warblers) != "Phylloscopus_trochiloides"]
for(number in c('1','2','3')){
  Y_sub <- Y[rownames(Y)[grep(paste0("_",number,"$"), rownames(Y))],,drop=F]
  if(any(colSums(Y_sub, na.rm =T)<5)){
    print(paste0('In atlas ',number,' these species: '))
    print(names(which(colSums(Y_sub, na.rm =T)<5)))
    tofilter <- names(which(colSums(Y_sub, na.rm =T)<5))
    print('have less than 5 occurrences')
    Y <- Y[, !(colnames(Y) %in% tofilter)]
    print('and are now filtered ')
  }
}

for(number in c('1','2','3')){
  Y_sub <- Y[rownames(Y)[grep(paste0("_",number,"$"), rownames(Y))],,drop=F]
  if(any(colSums(Y_sub, na.rm =T)<5)){
    print(paste0('In atlas ',number,' these species: '))
    print(names(which(colSums(Y_sub, na.rm =T)<5)))
    tofilter <- names(which(colSums(Y_sub, na.rm =T)<5))
    print('have less than 5 occurrences')
    stop('stop')
  }
}


# fixed
# --> LOAD TRAITS 
Tr <- read.csv(file.path(input,"data/1_preprocessing/Tr_aits/traits-guild_migration.csv"),row.names = 2)[,c(2,3)]
# sort by Y
Tr <- Tr[colnames(Y), , drop = FALSE]

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

xycoords <- data.frame(lon = Design$lon, lat = Design$lat)
#rownames(xycoords) <- Design3$site

v <- vect(xycoords, geom = c("lon","lat"), crs = "EPSG:4326")
v_proj <- project(v, "EPSG:23032")

proj_xycoords <- crds(v_proj)
Design$lon <- proj_xycoords[,1]
Design$lat <- proj_xycoords[,2]
# plot 
plot(Design$lat~Design$lon)

# load phylo
# phylogeny
phy <- read.tree(file.path(input,'data/1_preprocessing/Taxonomy/tree_fromPD.tre'))
phy <- keep.tip(phy,colnames(Y))

pd_matrix <- cophenetic.phylo(phy)
pd_matrix <- pd_matrix[sort(rownames(pd_matrix)), sort(colnames(pd_matrix))]
pd_matrix

# check if species-lists are identical
setdiff(rownames(pd_matrix),names(Y))
# stunning 

# PREPARING MODEL BUILD ---------------------------------------------------
# Define model types: 
date <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")

atlasnr <- c('1','2')
# loop over different atlases 
for(atlasnr in atlases){
  ### SAME ACROSS ALL MODELS
  # Define model formulas for environmental and trait data
  XFormula <- as.formula(paste("~", paste(colnames(X), collapse = "+"), sep = " "))
  TrFormula <- as.formula(paste("~", paste(colnames(Tr), collapse = "+"), sep = " "))
  
  # get random effects for space
  proj_xycoords_unique <- distinct(data.frame(X = Design$lon,
                                              Y = Design$lat))
  rownames(proj_xycoords_unique) <- unique(Design$site) 
  struc_space <- HmscRandomLevel(sData = proj_xycoords_unique, sMethod = "Full")
  
  # keep only atlas 1,2,3 
  pattern <- paste0("_(", paste(atlasnr, collapse = "|"), ")$")
  
  Y_sub <- Y[rownames(Y)[grep(pattern, rownames(Y))],,drop=F]
  X_sub <- X[rownames(X)[grep(pattern, rownames(X))],,drop=F]
  Design_sub <- Design[rownames(Design)[grep(pattern, rownames(Design))],,drop=F]
  Design_sub$atlas <- droplevels(Design_sub$atlas)
  #Design_sub <- Design_sub[,c('site','year'),drop=F]
  
  # and time
  years_unique <- distinct(data.frame(Year = Design_sub$year))
  rownames(years_unique) <- unique(Design_sub$atlas) 
  struc_time <- HmscRandomLevel(sData = years_unique, sMethod = "Full")
  
  m <-Hmsc(Y = Y_sub, 
           XData = X_sub,
           XFormula = XFormula,
           TrData = Tr,
           TrFormula = TrFormula,
           phyloTree = phy,
           studyDesign = Design_sub[,c('site','atlas'),drop=F], 
           ranLevels = list('site'=struc_space,
                            'atlas'=struc_time),
           distr='probit')
  
  
  print(head(Y_sub[,1:5]))
  print(tail(Y_sub[,1:5]))
  
  summary(m)
  m$rLNames
  m$ranLevels
  m$ranLevels$site
  m$ranLevels$atlas
  m$studyDesign
  ### IN RSTUDIO START SAMPLING 
  if(flagFitR){
    print('Rstudio stuff executed')
    init_obj <- sampleMcmc(m, samples=nSamples, thin=thin,
                           transient=transient, nChains=nChains,
                           verbose = verbose,
                           engine="HPC")
  }
  ### IN HPC ENVIORNMENT SET UP INIT
  if(flagInit){
    # initiate mcmc sampling 
    init_obj <- sampleMcmc(m, samples=nSamples, thin=thin,
                           transient=transient, nChains=nChains,
                           verbose = verbose,
                           engine="HPC")
    
    dir_name <- paste0(date,'_threeenv_allspecies_atlas_',  paste(atlasnr, collapse = ""))
    dir.create(file.path(input,'tmp_rds',dir_name))
    
    init_file_path = file.path(input,'tmp_rds',dir_name, "init_file.rds")
    m_file_path = file.path(input,'tmp_rds',dir_name, "m_object.rds")
    param_file_path = file.path(input,'tmp_rds',dir_name, "params.rds")
    lines <- paste(names(params), format(unlist(params), scientific = FALSE, trim=T), sep = "=")
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
                            "--transient", format(transient,scientific=F),
                            "--thin", thin,
                            "--verbose", verbose)
    cat(paste(shQuote(python), python_cmd_args), "\n")
    print('Init files created')
    
    
    # --- Auto-backup of the running script ----------------------------------------
    
    # Try to detect the current script path
    args <- commandArgs(trailingOnly = FALSE)
    script_path <- NULL
    
    if ("--file=" %in% substr(args, 1, 7)) {
      # Works for Rscript or batch jobs
      script_path <- normalizePath(sub("--file=", "", args[grep("--file=", args)]))
    } else if (!is.null(sys.frames()) && !is.na(sys.frames()[[1]]$ofile)) {
      # Fallback for interactive runs in RStudio
      script_path <- normalizePath(sys.frames()[[1]]$ofile)
    }
    
    if (!is.null(script_path) && file.exists(script_path)) {
      backup_file <- file.path(
        file.path(input,'tmp_rds',dir_name),
        paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S_"), basename(script_path))
      )
      
      file.copy(script_path, backup_file, overwrite = TRUE)
      message("✅ Script copied to: ", backup_file)
    } else {
      warning("⚠️ Could not determine script path — are you running interactively?")
    }
  }
  
}

