rm(list = ls())

# Define MCMC settings
env_vars <- c('tmean_year','prec_year','hh')
characteristics <- c('migration','guild','taxonomy')

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

# keep only sites with data  
X <- na.omit(X)

# keep only atlas 3
atlas3 <- X[rownames(X)[grep("_3$", rownames(X))],,drop=F]

# get ocean thresholds
grids_thresholds <- st_read(file.path(input,'data/1_preprocessing/atlas-grids/grids-ocean-thresholds/grids_ocean_thresholds.shp'))
thresholds <- grids_thresholds$kvdrtkd[grids_thresholds$pct_lnd>=25]

atlas3 <- atlas3[sub("_[123]$", "", rownames(atlas3)) %in% thresholds,]

# grab tmean_year
sites_actual <- row.names(atlas3)

### OCCURRENCES 
Y <- read.csv(file.path(input,'data/1_preprocessing/Y_occurrences/Y_occurrences.csv'),row.names=1)

# remove sites without data 
Y <- Y[row.names(Y) %in% sites_actual,]

# grab 13 warblers
genera <- c('Phylloscopus','Curruca','Sylvia','Acrocephalus','Hippolais','Locustella')
keep <- sapply(strsplit(colnames(Y),'_'),head,1) %in% genera
Y_warblers <- Y[,keep]

# check if they have enough records 
colSums(Y_warblers, na.rm =T)
# barred warbler is absent so we remove it 
Y_warblers <- Y_warblers[colnames(Y_warblers) != "Curruca_nisoria"]
# phylloscopus also very absent in some of the thresholds now 
#Y_warblers <- Y_warblers[colnames(Y_warblers) != "Phylloscopus_trochiloides"]
if(any(colSums(Y_warblers, na.rm =T)<5)){
  stop('At least one species has less than 5 occurrences')
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
i <- 'guild'
for(i in characteristics){
  ### SAME ACROSS ALL MODELS
  # Define model formulas for environmental and trait data
  
  X <- atlas3[,env_vars,drop=F]
  XFormula <- as.formula(paste("~", paste(colnames(X), collapse = "+"), sep = " "))
  
  struc_space <- HmscRandomLevel(sData = proj_xycoords, sMethod = "Full")
  
  # set distances to min and max distances between sites 
  # freq <- c(0.5,rep(0.005,100))
  # samples <- c(0,seq(from = 4999, to = 477312, length.out = 100))
  # small <- cbind(samples,freq)
  # 
  # struc_space_small <- setPriors(struc_space,alphapw=small)
  mig_guild <- list('migration' = 'Migration_a3_DOF',
                   'guild' = 'foraging_guild_consensus')
  if(i %in% c('migration','guild')){
    
    Tr_solo <- Tr[,mig_guild[[i]],drop=F]
    Tr_solo <- Tr_solo[rownames(Tr_solo) %in% colnames(Y_warblers),,drop=F]
    
    Tr_solo
    
    TrFormula <- as.formula(paste("~", paste(colnames(Tr_solo), collapse = "+"), sep = " "))
    
    # init the model 
    m <-Hmsc(Y = Y_warblers, 
             XData = X,
             XFormula = XFormula,
             TrData = Tr_solo,
             TrFormula = TrFormula,
             studyDesign = Design3[,c('site'),drop=F], 
             ranLevels = list('site'=struc_space),
             distr='probit')
  }else if(i %in% c('taxonomy')){
    # phylogeny
    phy <- read.tree(file.path(input,'data/1_preprocessing/Taxonomy/tree_fromPD.tre'))
    phy_warblers <- keep.tip(phy,colnames(Y_warblers))
    plot(phy_warblers)

    pd_matrix <- cophenetic.phylo(phy_warblers)
    pd_matrix <- pd_matrix[sort(rownames(pd_matrix)), sort(colnames(pd_matrix))]
    pd_matrix

    # check if species-lists are identical
    setdiff(rownames(pd_matrix),names(Y))
    # stunning 
    
    # init the phylo model 
    m <-Hmsc(Y = Y_warblers, 
             XData = X,
             XFormula = XFormula,
             phyloTree = phy_warblers,
             studyDesign = Design3[,c('site'),drop=F], 
             ranLevels = list('site'=struc_space),
             distr='probit')
  }
    
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

dir_name <- paste0(date,'_threeenv_',i)
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
