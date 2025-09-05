mod <- '2025-08-29_15-52-09_samples_1000_thin_100'
easy_mode <- 0

# GETTING STARTED ---------------------------------------------------------
if (interactive() && Sys.getenv("RSTUDIO") == "1") {
  message("Running in RStudio")
  rstudio = easy_mode
  library(Hmsc)
  library(jsonify)
  library(knitr)
  library(corrplot)
  library(colorspace)
  library(vioplot)
  library(dplyr)
  input <- file.path('./tmp_rds/mods-single',mod)
  source_path <- file.path('./scripts/3_modeldiagnostics/plotting-scripts')
  
} else {
  message("Running from terminal or non-interactive environment")
  rstudio = 0
  library(RColorBrewer,lib="~/Rlibs")
  library(farver,lib="~/Rlibs")
  library(scales,lib="~/Rlibs")
  library(jsonify,lib="~/Rlibs")
  library(ape,lib="~/Rlibs")
  library(dplyr,lib="~/Rlibs")
  library(Hmsc,lib="~/Rlibs")
  library(colorspace,lib="~/Rlibs")
  library(vioplot,lib="~/Rlibs")
  input <- file.path('~/home/projects/hmsc-danishbirds/tmp-rds',mod)
  source_path <- file.path('~/home/projects/hmsc-danishbirds/scripts/3_modeldiagnostics/plotting-scripts')
  
}

# LOADING DATA --------------------------------------------------------
# load unfitted object
m <- readRDS(file.path(input,'m_object.rds'))
summary(m)
# load params 
params <- readRDS(file.path(input,'params.rds'))
nChains <- params$nChains
nSamples <- params$nSamples
thin <- params$thin
transient <- params$transient

# loading the chains 
chainList = vector("list", nChains)
for(cInd in 1:nChains){
  chain_file_path = file.path(input, sprintf("post_chain%.2d_file.rds", cInd-1))
  print(chain_file_path)
  if(file.exists(chain_file_path)) {
    chainList[[cInd]] = from_json(readRDS(file = chain_file_path)[[1]])[[1]]
  }
}

# use only 4 chains if running on local computer
if(rstudio == 1){
  filteredList <- list(chainList[[1]],chainList[[2]])
}else if(rstudio == 0){
  filteredList <- chainList
}

fitSepTF = importPosteriorFromHPC(m, filteredList, nSamples, thin, transient)

# PSRF ----------------------------------------------------
# Convert model output to coda format for MCMC diagnostics
# make smaller for ease in R

if(rstudio == 1){
  start <- nSamples - 1
  mpost <- convertToCodaObject(fitSepTF,start=start)
}else if(rstudio == 0){
  mpost <- convertToCodaObject(fitSepTF,start=1)
}

# PLOTTING PARAMETERS 
# official names 
names <- c('beta','gamma','rho','V',
           'alpha_space','alpha_time',
           'omega_raw','omega_resid')

# plot names 
full_names <- c("Beta - Species x Environment",
                "Gamma - Traits x Environment",
                "Rho - Phylogeny",
                "V - Variation",
                "Alpha - Spatial",
                "Alpha - Time",
                "Omega1 - Raw associations",
                "Omega2 - Residual associations")

# plot colors 
cols <- c('green4','purple4','yellow4','orange4','ivory','ivory','blue4','blue4')

# objects in the list 
params <- list(mpost$Beta,mpost$Gamma,mpost$Rho,mpost$V,
               mpost$Alpha[[1]],mpost$Alpha[[2]],
               mpost$Omega[[1]],mpost$Omega[[2]])

diags <- list(psrf = list(), ess = list())
chunk_size <- 10

for(j in seq_along(params)){
  print(j)
  mat <- params[[j]]
  ncols <- ncol(mat[[1]]) # get ncol 
  
  # process in chunks if ncols > 1
  if(length(ncols)){
    idx_chunks <- split(1:ncols, ceiling(seq_along(1:ncols)/chunk_size))
    
    psrf_list <- list()
    ess_list  <- list()
    
    # parallelize
    results <- parallel::mclapply(idx_chunks, function(cols) {
      list(
        psrf = gelman.diag(mat[,cols], multivariate=FALSE)$psrf[,1],
        ess  = effectiveSize(mat[,cols])
      )
    }, mc.cores = 4)  # adjust cores to your HPC
    
    # Combine
    psrf_list <- lapply(results, `[[`, "psrf")
    ess_list  <- lapply(results, `[[`, "ess")
    
    diags$psrf[[j]] <- unlist(psrf_list)
    diags$ess[[j]]  <- unlist(ess_list)
  }else{
    diags$psrf[[j]] <- gelman.diag(mat, multivariate=FALSE)$psrf[,1]
    diags$ess[[j]] <- effectiveSize(mat)
  }
  
}

source(file.path(source_path,'psrf-ess-plots.R'))
source(file.path(source_path,'psrf-ess-singles.R'))

