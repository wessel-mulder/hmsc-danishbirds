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
taxonomy_flag <- args[9]
spatial_flag <- args[10]
temporal_flag <- args[11]


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
  mod <- '2025-09-29_20-30-44_singleev_tmean_year_knots_212'
  input <- file.path('./tmp_rds/mods-complexity',mod)
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

#params <- list(mpost$Beta,mpost$Gamma,mpost$Rho,mpost$V,mpost$Sigma,
 #              mpost$Eta[[1]],mpost$Eta[[2]],
 #              mpost$Alpha[[1]],mpost$Alpha[[2]],
 #              mpost$Omega[[1]],mpost$Omega[[2]],
 #              mpost$Lambda[[1]],mpost$Lambda[[2]],
 #              mpost$Psi[[1]],mpost$Psi[[2]],
 #              mpost$Delta[[1]],mpost$Delta[[2]])
params <- list(
  beta = mpost$Beta,
  gamma = mpost$Gamma,
  V = mpost$V,
  sigma = mpost$Sigma
)
 if(taxonomy_flag == 1){
  params$Rho <- mpost$Rho
 }
  if(spatial_flag == 1){
    params$eta <- mpost$Eta[[1]]
    params$alpha = mpost$Alpha[[1]]
    params$omega = mpost$Omega[[1]]
    params$lambda = mpost$Lambda[[1]]
    params$psi = mpost$Psi[[1]]
    params$delta = mpost$Delta[[1]]
  }
 if(temporal_flag == 1){
    params$eta2 <- mpost$Eta[[2]]
    params$alpha2 = mpost$Alpha[[2]]
    params$omega2 = mpost$Omega[[2]]
    params$lambda2 = mpost$Lambda[[2]]
    params$psi2 = mpost$Psi[[2]]
    params$delta2 = mpost$Delta[[2]]
 }


diags <- list(psrf = list(), ess = list())
chunk_size <- 10

for (name in names(params)) {
  cat("Processing:", name, "\n")
  mat <- params[[name]]
  
  # Handle objects with and without columns
  ncols <- tryCatch(ncol(mat[[1]]), error = function(e) NULL)
  
  if (length(ncols) && !is.null(ncols)) {
    idx_chunks <- split(1:ncols, ceiling(seq_along(1:ncols)/chunk_size))
    
    results <- parallel::mclapply(idx_chunks, function(cols) {
      list(
        psrf = gelman.diag(mat[, cols], multivariate = FALSE, autoburnin = FALSE)$psrf[, 1],
        ess  = effectiveSize(mat[, cols])
      )
    }, mc.cores = n_cores)
    
    diags$psrf[[name]] <- unlist(lapply(results, `[[`, "psrf"))
    diags$ess[[name]]  <- unlist(lapply(results, `[[`, "ess"))
    
  } else {
    diags$psrf[[name]] <- gelman.diag(mat, multivariate = FALSE, autoburnin = FALSE)$psrf[, 1]
    diags$ess[[name]]  <- effectiveSize(mat)
  }
}

# init pdf 
saveRDS(diags,file=psrf_ess_output)
print('psrf and ess succesfully saved')
}


# AUC-TJUR ----------------------------------------------------------------
fit_output <- file.path(input,'model-outputs','model-fit.rds')
print('starting model-fit')
### HEAVY 
if (fit_flag == 1) {
preds_expected  <- pcomputePredictedValues(fitSepTF,expected=T)
#preds  <- pcomputePredictedValues(fitSepTF,expected=F)

mean_expected <- apply(preds_expected,c(1,2),mean)
#mean <- apply(preds,c(1,2),mean)

saveRDS(mean_expected,file=file.path(input,'model-outputs','pred-vals.rds'))

MF <- evaluateModelFit(hM=fitSepTF, predY=preds_expected)
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

if (VP_flag == 1) {
VP = computeVariancePartitioning(fitSepTF)
saveRDS(VP,file=file.path(input,'model-outputs','VP.rds'))
print('variance partitions succesfully saved')
}else{
  print('variance partitions skipped')
}

# POSTERIOR ESTIMATES  ----------------------------------------------------

if(post_estimates_flag==1){
  print('starting posteriors')
  if(spatial_flag == 1){
    param_list <- c("Beta", "Gamma", "V", "rho", "sigma","wRRR","Eta", "Lambda", "Psi", "Delta","Alpha","Omega","OmegaCor")
  }
  for(parameter in param_list){
  posterior = getPostEstimate(fitSepTF, parName = parameter)
  if(!dir.exists(file.path(input,'model-outputs','posterior-estimates'))) {dir.create(file.path(input,'model-outputs','posterior-estimates'))}
  saveRDS(posterior,file=file.path(input,'model-outputs','posterior-estimates',paste0('posterior-',parameter,'.rds')))
  }
  print('posteriors finished')
}else{
  print('posteriors skipped')
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



# SPATIAL PREDICTION ------------------------------------------------------

if(sp_pred_flag==1){
  grid <- fitSepTF$XData
  covariates <- colnames(grid)
  
  # get design and xycoords of sites  
  design <- fitSepTF$studyDesign
  xycoords <- tryCatch(
    {
      # First attempt
      data.frame(fitSepTF$ranLevels$site$s@coords[drop = FALSE])
    },
    error = function(e) {
      # Fallback if the above errors
      data.frame(fitSepTF$ranLevels$site$s)
    }
  )
  
  head(xycoords)
  # rename X and Y consistentely 
  col_with_max <- colnames(xycoords)[which.max(lapply(xycoords, max))]
  names(xycoords)[names(xycoords) == col_with_max] <- "Y"  # Rename it to "Y"
  col_with_min <- colnames(xycoords)[which.min(lapply(xycoords, min))]
  names(xycoords)[names(xycoords) == col_with_min] <- "X"  # Rename it to "Y"
  head(xycoords)
  # organize columns
  xycoords <- xycoords[,c('X','Y')]
  head(xycoords)
  # merge by rownmaes 
  merge <- merge(design,xycoords, by.x = 'site', by.y = 'row.names', all.x =T)
  rownames(merge) <- rownames(design)
  highlight <- merge$highlight <- grepl("BF", merge$site)
  
  # Plot all points in black
  #plot(merge[,c('X','Y')])
  # Add highlighted points in red on top
  #points(merge$X[highlight], merge$Y[highlight], col = "red")
  
  merge_env <- merge(merge,grid,by='row.names',all.x=T)
  # Plot all points in black
  #plot(merge_env[,c('X','Y')])
  # Add highlighted points in red on top
  # overlay highlighted sites, colored by temperature
  #points(merge_env$X[highlight], merge_env$Y[highlight],
  #       pch = 19, col = heat.colors(100)[cut(merge_env$tmean_year[highlight], 100)])

  #construct gradient 
  chunk_size <- 10
  nr_prds <- nrow(merge_env)
  
  idx_chunks <- split(1:nr_prds, ceiling(seq_along(1:nr_prds)/chunk_size))

    # parallelize
  print(covariates)
  results <- parallel::mclapply(idx_chunks, function(preds_subset) {
    print(preds_subset)
    data_subset <- merge_env[preds_subset,] # subset data to chunks 
    
    # grab coordinates for 'new' sites predictions 
    xy_subset <- data_subset[,c('X','Y')]
    env_subset <- data_subset[,covariates,drop=F]
    design_subset <- data_subset[,'Row.names']
    
    gradient_spat <- prepareGradient(fitSepTF,
                                     XDataNew = env_subset,
                                     sDataNew = list(site = xy_subset))
    predY = predict(fitSepTF,Gradient = gradient_spat,
                    predictEtaMean=F,predictEtaMeanField=F,
                    expected=F)
    EpredY = Reduce("+",predY)/length(predY)
    rownames(EpredY) <- design_subset
    EpredY

  },mc.cores=n_cores)
  
  combined <- do.call(rbind,results)
  saveRDS(combined,file=file.path(input,'model-outputs','spatialpreds.rds'))
  print('predictions finished')

}else{
  print('predictions skipped')
}







