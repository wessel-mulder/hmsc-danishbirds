args <- commandArgs(trailingOnly = TRUE)

mod = args[1]
psrfess_flag <- 0 
fit_flag <- 1
VP_flag <- 1
pred_flag <- 1


# GETTING STARTED ---------------------------------------------------------
if (interactive() && Sys.getenv("RSTUDIO") == "1") {
  message("Running in RStudio")
  library(Hmsc)
  library(jsonify)
  library(knitr)
  library(corrplot)
  library(colorspace)
  library(vioplot)
  library(dplyr)
  mod <- '2025-09-08_17-32-13_samples_1000_thin_100'
  input <- file.path('./tmp_rds/mods-single',mod)
  source_path <- file.path('./scripts/3_modeldiagnostics/plotting-scripts')
  
} else {
  message("Running from terminal or non-interactive environment")
  library(RColorBrewer,lib="~/Rlibs")
  library(farver,lib="~/Rlibs")
  library(scales,lib="~/Rlibs")
  library(jsonify,lib="~/Rlibs")
  library(ape,lib="~/Rlibs")
  library(dplyr,lib="~/Rlibs")
  library(Hmsc,lib="~/Rlibs")
  library(colorspace,lib="~/Rlibs")
  library(vioplot,lib="~/Rlibs")
  input <- file.path('~/home/projects/hmsc-danishbirds/tmp_rds',mod)
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



chainList = vector("list", nChains)
for(cInd in 1:nChains){
  chain_file_path = file.path(input, sprintf("post_chain%.2d_file.rds", cInd-1))
  print(chain_file_path)
  if(file.exists(chain_file_path)) {
    chainList[[cInd]] = from_json(readRDS(file = chain_file_path)[[1]])[[1]]
  }
}

filteredList <- chainList

fitSepTF = importPosteriorFromHPC(m, filteredList, nSamples, thin, transient)

print('model succesfully loaded')
# PSRF / ESS  ----------------------------------------------------
# Convert model output to coda format for MCMC diagnostics
# make smaller for ease in R
if(psrfess_flag == 1){
mpost <- convertToCodaObject(fitSepTF)

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
    }, mc.cores = 4)  # adjust cores to liking
    
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

# init pdf 
if(!dir.exists(file.path(input,'model-outputs'))) {dir.create(file.path(input,'model-outputs'))}
saveRDS(diags,file=file.path(input,'model-outputs','psrf-ess.rds'))

print('psrf and ess succesfully saved')
}else{
  print('psrf and ess skipped')
}
# AUC-TJUR ----------------------------------------------------------------
preds  <- computePredictedValuesParallel(fitSepTF,expected=T)
print('preds succesfully made')
MF <- evaluateModelFit(hM=fitSepTF, predY=preds)
saveRDS(MF,file=file.path(input,'model-outputs','model-fit.rds'))
print('model fit succesfully saved')

# VP ----------------------------------------------------------------
VP = computeVariancePartitioning(fitSepTF)

# split by groups 
#names <- VP$groupnames
VP_split = computeVariancePartitioning(fitSepTF,
                                       group = c(1,1,1,
                                                 2,2,2,
                                                 3,3,
                                                 rep(4,8)),
                                       groupnames = c('temperature',
                                                      'precipitation',
                                                      'landscape',
                                                      'land-use classes'))
# split by seasons
VP_season = computeVariancePartitioning(fitSepTF,
                                        group = c(1,2,3,
                                                  1,2,3,
                                                  4,4,
                                                  rep(5,8)),
                                        groupnames = c('year',
                                                       'winter',
                                                       'breeding',
                                                       'landscape',
                                                       'land-use classes'))

saveRDS(VP,file=file.path(input,'model-outputs','VP-full.rds'))
saveRDS(VP_split,file=file.path(input,'model-outputs','VP-split.rds'))
saveRDS(VP_season,file=file.path(input,'model-outputs','VP-season.rds'))

print('variance partitions succesfully saved')


# TEST A PREDICTION -------------------------------------------------------
covariates <- c('tmean_year','prec_year')
mclapply(covariates, function(covariate) {

  print('starting gradient')
  Gradient <- constructGradient(m, focalVariable = covariate, ngrid = 30)
  print('starting predictions')

  predY <- predict(fitSepTF, Gradient = Gradient, expected = TRUE)
  
  # Save each covariate separately
  saveRDS(list(predY = predY, Gradient = Gradient),
          file = file.path(input,'model-outputs', paste0("pred_", covariate, ".rds")))
  
}, mc.cores = 2)

print('predictions succesfully saved')



