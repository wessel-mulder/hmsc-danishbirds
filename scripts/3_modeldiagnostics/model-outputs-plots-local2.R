
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
  mod <- '2025-09-10_11-14-06_samples_250_thin_1000'
  input <- file.path('./tmp_rds/mods-single',mod)
  source_path <- file.path('./scripts/3_modeldiagnostics/plotting-scripts')
  # other flags
  psrfess_flag <- 1
  fit_flag <- 1
  VP_flag <- 1
  pred_flag <- 1
  sp_pred_flag <- 1
  post_estimates_flag <- 1
  n_cores <- 6
  
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
if (psrfess_flag == 0 || file.exists(psrf_ess_output)){
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

diags <- readRDS(psrf_ess_output)
lapply(diags, function(x) lapply(x, summary))

source(file.path(source_path,'psrf-ess-plots.R'))

#source(file.path(source_path,'psrf-ess-singles.R'))

# Beta 
mpost <- convertToCodaObject(fitSepTF)
# get indices of worst chains 
pdf(file = file.path(input,'results','chains.pdf'),
    width = 6,
    height = 8)
for(i in 1:17){
print(i)
params <- list(mpost$Beta,mpost$Gamma,mpost$Rho,mpost$V,mpost$Sigma,
               mpost$Eta[[1]],mpost$Eta[[2]],
               mpost$Alpha[[1]],mpost$Alpha[[2]],
               mpost$Omega[[1]],mpost$Omega[[2]],
               mpost$Lambda[[1]],mpost$Lambda[[2]],
               mpost$Psi[[1]],mpost$Psi[[2]],
               mpost$Delta[[1]],mpost$Delta[[2]])
if(length(diags$psrf[[i]])>3){
top3_idx <- order(diags$psrf[[i]], decreasing = TRUE)[1:3]
plot(params[[i]][,top3_idx])
}else{
  plot(params[[i]])
}
}
dev.off()

fitSepTF

pdf(file = file.path(input,'results','chains-alpha.pdf'),
    width = 6,
    height = 8)
traceplot(mpost$Alpha[[1]])
traceplot(mpost$Alpha[[2]])
traceplot(mpost$Alpha[[3]])


dev.off()

top3_idx <- order(diags$psrf[[1]], decreasing = F)[1:3]
plot(mpost$Beta[,top3_idx])

# CHAINS ------------------------------------------------------------------

summary(mpost$Alpha[[1]])



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

MF <- readRDS(fit_output)
mean(MF$TjurR2)
source(file.path(source_path,'auc-tjur-plots.R'))
# VP ------------------------------------------------
----------------
VP_full_output <- file.path(input,'model-outputs','VP-full.rds')
VP_split_output <- file.path(input,'model-outputs','VP-split.rds')
VP_season_output <- file.path(input,'model-outputs','VP-season.rds')

VP_output <- c(VP_full_output, VP_split_output, VP_season_output)

if (fit_flag == 1 || all(!file.exists(VP_output))) {
  VP = computeVariancePartitioning(fitSepTF,start=150)
  
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

readRDS(VP)

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

fitSepTF
# TRACE PLOTS -------------------------------------------------------------
mpost$Eta[[1]]
postEta <- getPostEstimate(fitSepTF,parName='Eta')
xy <- fitSepTF$rL[[1]]$s

mpost$Delta[[1]]
t <- mpost$Delta[[1]]
summary(mpost$Delta[[1]])
summary(mpost$Alpha[[1]])

fitSepTF$XData[[1]]


traceplot(mpost$Eta[[1]][,1:2])

library(ggplot2)
for(i in 1:ncol(postEta$mean)){
  eta <- postEta$mean[,i]
  eta2 <- cbind(eta,xy)
  p<-ggplot(eta2,aes(x=lon,y=lat,col=eta))+
    geom_point(cex=3)+
    labs(title=i)+
    scale_colour_gradient2(
      low = "blue",      # color for negative values
      mid = "ivory",     # color for zero
      high = "red",      # color for positive values
      midpoint = 0
    )
  print(p)
}


### CLIMATE DATA 
landuse <- 'tmean_year'
# check spatial distribution
tmean <- data.frame(temp=fitSepTF$X[,landuse,drop=F])
colnames(tmean) <- 'env'
# keep only rows where rownames end in "_3"
tmean_a3 <- tmean[grepl(paste0("_",3,"$"), rownames(tmean)), , drop = FALSE]
# strip the "_3" from the rownames
rownames(tmean_a3) <- sub(paste0("_",3,"$"), "", rownames(tmean_a3))
#mergs
tmean_space <- merge(xy,tmean_a3,by='row.names')
# plot 
p<-ggplot(tmean_space,
          aes(x=lon,y=lat,color=env))+
  geom_point(cex=3)+
  scale_color_gradientn(colors=c('goldenrod','coral','firebrick'))+
  labs(title='Yearly temperatures - Atlas III')
p

landuse <- 'prec_year'
# check spatial distribution
tmean <- data.frame(temp=fitSepTF$X[,landuse,drop=F])
colnames(tmean) <- 'env'
# keep only rows where rownames end in "_3"
tmean_a3 <- tmean[grepl(paste0("_",3,"$"), rownames(tmean)), , drop = FALSE]
# strip the "_3" from the rownames
rownames(tmean_a3) <- sub(paste0("_",3,"$"), "", rownames(tmean_a3))
#mergs
tmean_space <- merge(xy,tmean_a3,by='row.names')
# plot 
p<-ggplot(tmean_space,
          aes(x=lon,y=lat,color=env))+
  geom_point(cex=3)+
  scale_color_gradientn(colors=c('lightblue','blue','darkblue'))+
  labs(title='Yearly precipitation (mm) - Atlas III')
p

landuse <- 'LULC_0'
# check spatial distribution
tmean <- data.frame(temp=fitSepTF$X[,landuse,drop=F])
colnames(tmean) <- 'env'
# keep only rows where rownames end in "_3"
tmean_a3 <- tmean[grepl(paste0("_",3,"$"), rownames(tmean)), , drop = FALSE]
# strip the "_3" from the rownames
rownames(tmean_a3) <- sub(paste0("_",3,"$"), "", rownames(tmean_a3))
#mergs
tmean_space <- merge(xy,tmean_a3,by='row.names')
# plot 
p<-ggplot(tmean_space,
          aes(x=lon,y=lat,color=env))+
  geom_point(cex=3)+
  scale_color_gradientn(colors=c('white','darkblue'))+
  labs(title='Ocean proportional coverage - Atlas III')
p

