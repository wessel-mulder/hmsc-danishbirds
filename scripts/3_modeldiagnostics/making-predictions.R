mod <- '2025-08-22_10-09-28_samples_250_thin_100'
# GETTING STARTED ---------------------------------------------------------
if (interactive() && Sys.getenv("RSTUDIO") == "1") {
  message("Running in RStudio")
  library(Hmsc)
  library(cli)
  input <- file.path('./tmp_rds/mods-tuning-hypeparam',mod)
   
} else {
  message("Running from terminal or non-interactive environment")
  library(Hmsc,lib="~/Rlibs")
  input <- file.path('~/home/projects/hmsc-danishbirds/tmp-rds',mod)
  
}

# init output dir 
output_dir <- file.path(input,'predictions')
if(!dir.exists(output_dir)){dir.create(output_dir)}

# log file
log_file <- file.path(output_dir, "predictions_log.txt")
# Open a connection to the log file
log_con <- file(log_file, open = "wt")

# Redirect stdout
sink(log_con, type = "output", split = TRUE)

# SET UP LOGFILES  --------------------------------------------------------



# LOADING THE MODEL -------------------------------------------------------
# load unfitted object
m <- readRDS(file.path(input,'m_object.rds'))
# load params 
params <- readRDS(file.path(input,'params.rds'))
nChains <- params$nChains
nSamples <- params$nSamples
thin <- params$thin
transient <- params$transient
# loading the chains 
# chainList = vector("list", nChains)
# for(cInd in 1:nChains){
#   chain_file_path = file.path(input, sprintf("post_chain%.2d_file.rds", cInd-1))
#   if(file.exists(chain_file_path)) {
#     chainList[[cInd]] = from_json(readRDS(file = chain_file_path)[[1]])[[1]]
#   }
# }
# check for any mistakes in chains and remove if necessary 
#filteredList <- chainList[!sapply(chainList, is.null)]
cat(paste0(length(filteredList),' chains loaded succesfully \n'))
# import posterior 
#fitSepTF = importPosteriorFromHPC(m, filteredList, nSamples, thin, transient)
cat('Posterior succesfully imported \n')

# PREDICTIONS ACROSS GRADIENTS --------------------------------------------
# construct gradient 
covariates <- colnames(fitSepTF$XData)
if(length(covariates)>0){
  # load output file 
  outfile <- file.path(output_dir,'preds.Rdata')
  
  cat("Making predictions \n")
  if(file.exists(file.path(outfile))){
    cat("Predictions already calculated \n")
    load(outfile)
  } else {
    preds <- vector('list',length(covariates))
    mclapply(covariates, function(covariate) {
      ptm <- proc.time()  # start timing
      
      Gradient <- constructGradient(m, focalVariable = covariate, ngrid = 30)
      predY <- predict(fitSepTF, Gradient = Gradient, expected = TRUE)
      
      # Save each covariate separately
      saveRDS(list(predY = predY, Gradient = Gradient),
              file = file.path(output_dir, paste0("Pred_", covariate, ".rds")))
      
      #Compute elapsed time
      computational.time <- proc.time() - ptm
      cat(sprintf("[%s] Time taken: %.2f s | Current time: %s\n",
                  covariate, computational.time[3], format(Sys.time(), "%H:%M:%S")))
      
      
      return(covariate)  # optional: just to keep track in the list
    }, mc.cores = nCores)

  }
}

# Gradient = constructGradient(fitSepTF,focalVariable = "tmean_winter",ngrid=30)
# pred = predict(fitSepTF,Gradient = Gradient,expected = T)
# 
# length(unique(fitSepTF$TrData$foraging_guild_consensus))
# Gradient$XDataNew
# 
# predY = predict(fitSepTF, XData=Gradient$XDataNew, studyDesign=Gradient$studyDesignNew,
#                 ranLevels=Gradient$rLNew, expected=TRUE)
# plotGradient(fitSepTF, Gradient, pred=predY, measure="S", showData = TRUE)
# 
# preds <- computePredictedValues(fitSepTF,partition.sp = 1,thin=20)

# Close sinks when done
sink(type = "output")
close(log_con)
