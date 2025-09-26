args <- commandArgs(trailingOnly = TRUE)

mod = args[1]
n_cores = args[2]

# other flags
psrfess_flag <- args[3]
fit_flag <- args[4]
VP_flag <- args[5]
pred_flag <- args[6]
sp_pred_flag <- args[7]
post_estimates_flag <- args[8]


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
  library(abind)
  mod <- '2025-09-15_14-10-17_spaceonly'
  input <- file.path('./tmp_rds/mods-space-nonspace',mod)
  source_path <- file.path('./scripts/3_modeldiagnostics/plotting-scripts')
  # other flags
  psrfess_flag <- 1
  fit_flag <- 1
  VP_flag <- 1
  pred_flag <- 1
  sp_pred_flag <- 1
  post_estimates_flag <- 1
  n_cores <-4
  
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
  library(abind,lib="~/Rlibs")
  
  input <- file.path('~/home/projects/hmsc-danishbirds/tmp_rds',mod)
  source_path <- file.path('~/home/projects/hmsc-danishbirds/scripts/3_modeldiagnostics/plotting-scripts')
  
}
# make dir for outputs 
if(!dir.exists(file.path(input,'model-outputs'))) {dir.create(file.path(input,'model-outputs'))}
if(!dir.exists(file.path(input,'results'))) {dir.create(file.path(input,'results'))}


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
mpost <- convertToCodaObject(fitSepTF)


print('model succesfully loaded')
# PSRF / ESS  ----------------------------------------------------
# Convert model output to coda format for MCMC diagnostics
# make smaller for ease in R
psrf_ess_output<-file.path(input,'model-outputs','psrf-ess.rds')
if (psrfess_flag == 0 ){
  print('psrf and ess skipped')
}else{
mpost <- convertToCodaObject(fitSepTF)
summary(mpost)

# objects in the list 
params <- list(mpost$Beta,mpost$Gamma,mpost$Rho,mpost$V,mpost$Sigma,
               mpost$Eta[[1]],mpost$Eta[[2]],
               mpost$Alpha[[1]],mpost$Alpha[[2]],
               mpost$Omega[[1]],mpost$Omega[[2]],
               mpost$Lambda[[1]],mpost$Lambda[[2]],
               mpost$Psi[[1]],mpost$Psi[[2]],
               mpost$Delta[[1]],mpost$Delta[[2]])

diags <- list(psrf = list(), ess = list())
chunk_size <- 10

for(j in seq_along(params)){
  print(j)
  mat <- params[[j]]
  ncols <- ncol(mat[[1]]) # get ncol 
  print(ncols)

  # process in chunks if ncols > 1
  if(length(ncols)){
    idx_chunks <- split(1:ncols, ceiling(seq_along(1:ncols)/chunk_size))
    
    psrf_list <- list()
    ess_list  <- list()
    
    # parallelize
    results <- parallel::mclapply(idx_chunks, function(cols) {
      list(
        psrf = gelman.diag(mat[,cols], multivariate=FALSE,autoburnin=F)$psrf[,1],
        ess  = effectiveSize(mat[,cols])
      )
    }, mc.cores = n_cores)  # adjust cores to liking
    
    # Combine
    psrf_list <- lapply(results, `[[`, "psrf")
    ess_list  <- lapply(results, `[[`, "ess")
    
    diags$psrf[[j]] <- unlist(psrf_list)
    diags$ess[[j]]  <- unlist(ess_list)
  }else{
    diags$psrf[[j]] <- gelman.diag(mat, multivariate=FALSE,autoburnin=F)$psrf[,1]
    diags$ess[[j]] <- effectiveSize(mat)
  }
  
}

# init pdf 
saveRDS(diags,file=psrf_ess_output)
print('psrf and ess succesfully saved')
}

source(file.path(source_path,'psrf-ess-plots.R'))
source(file.path(source_path,'psrf-ess-singles.R'))

# AUC-TJUR ----------------------------------------------------------------
fit_output <- file.path(input,'model-outputs','model-fit.rds')

### HEAVY 
if (fit_flag == 1 || !file.exists(fit_output)) {
preds  <- pcomputePredictedValues(fitSepTF,expected=T)
MF <- evaluateModelFit(hM=fitSepTF, predY=preds)
saveRDS(MF,file=fit_output)

print('model fit succesfully saved')
}else{
  print('model fit skipped')
}
# ### CHUNKY
# chunk_size <- 10
# 
# # get postlist
# postlist <- poolMcmcChains(fitSepTF$postList)
# idx_chunks[1]
# summary(postlist[idx_chunks[[2]]])
# postlist[[2]]
# 
# length_pl <- length(postlist)
# 
# # get chunks
# idx_chunks <- split(1:length_pl, ceiling(seq_along(1:length_pl)/chunk_size))
# 
# # parallelize
# idx_chunks <- idx_chunks[1:2]
# 
# preds_chunk <- parallel::mclapply(idx_chunks, function(chunk) {
#   sub_postlist <- postlist[chunk]
#   predict(fitSepTF,sub_postlist)
# }, mc.cores = 4)  # adjust cores to liking
# 
# abind(preds_chunk,along=3)
# 
# MF <- evaluateModelFit(hM=fitSepTF, predY=preds)
# ns = dim(fitSepTF$Y)[2]
# RMSE = rep(NA,ns)
# for (i in 1:ns){
#   RMSE[i] = sqrt(mean((fitSepTF$Y[,i]-preds[,i])^2, na.rm=TRUE))
# }
# 
# median2 = function(x){return (median(x,na.rm=TRUE))}
# mean2 = function(x){return (mean(x,na.rm=TRUE))}
# 
# mPredY = matrix(NA, nrow=fitSepTF$ny, ncol=fitSepTF$ns)
# sel = fitSepTF$distr[,1]==3
# if (sum(sel)>0){
#   mPredY[,sel] = as.matrix(apply(abind(preds[,sel,,drop=FALSE], along=3),
#                                  c(1,2), median2))
# }
# sel = !fitSepTF$distr[,1]==3
# if (sum(sel)>0){
#   mPredY[,sel] = as.matrix(apply(abind(preds[,sel,,drop=FALSE], along=3),
#                                  c(1,2), mean2))
# }


# VP ----------------------------------------------------------------
VP_full_output <- file.path(input,'model-outputs','VP-full.rds')
VP_split_output <- file.path(input,'model-outputs','VP-split.rds')
VP_season_output <- file.path(input,'model-outputs','VP-season.rds')

VP_output <- c(VP_full_output, VP_split_output, VP_season_output)

if (fit_flag == 1 || all(!file.exists(VP_output))) {
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

saveRDS(VP,file=VP_full_output)
saveRDS(VP_split,file=VP_split_output)
saveRDS(VP_season,file=VP_season_output)

print('variance partitions succesfully saved')
}else{
  print('variance partitions skipped')
}

# TEST A PREDICTION -------------------------------------------------------
if(pred_flag==1){
covariates <- c('tmean_year','prec_year')
covariate <- covariates[1]
#parallel::mclapply(covariates, function(covariate) {

  print('starting gradient')
  Gradient <- constructGradient(m, focalVariable = covariate, ngrid = 5)
  print('starting predictions')

  predY <- predict(fitSepTF, Gradient = Gradient, expected = TRUE,
  nParallel = n_cores)
  
  # Save each covariate separately
  saveRDS(list(predY = predY, Gradient = Gradient),
          file = file.path(input,'model-outputs', paste0("pred_", covariate, ".rds")))
  
#}, mc.cores = 2)

print('predictions succesfully saved')
}else{
  print('predictions skipped')
}

# GET POSTERIOR ESTIMATES 
print('starting posteriors')
if(post_estimates_flag==1){
  for(parameter in c('Beta','Gamma','Omega','OmegaCor')){
  posterior = getPostEstimate(fitSepTF, parName = parameter)
  saveRDS(posterior,file=file.path(input,'model-outputs',paste0('posterior-',parameter,'.rds')))
  }
}



