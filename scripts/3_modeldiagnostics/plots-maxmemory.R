args <- commandArgs(trailingOnly = TRUE)

mod = args[1]
easy_mode <- 1

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
  mod <- '2025-09-08_17-32-13_samples_1000_thin_100'
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

# PSRF / ESS  ----------------------------------------------------
# Convert model output to coda format for MCMC diagnostics
# make smaller for ease in R

if(rstudio == 1){
  start <- 1
}else if(rstudio == 0){
  start <- 1
} 

mpost <- convertToCodaObject(fitSepTF,start=start)

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

print(diags$psrf[[3]])
print(diags$ess[[3]])

source(file.path(source_path,'psrf-ess-plots.R'))
source(file.path(source_path,'psrf-ess-singles.R'))

# AUC / TJUR --------------------------------------------------------------
# get preds 
preds  <- computePredictedValues(fitSepTF, start = 1)
MF <- evaluateModelFit(hM=fitSepTF, predY=preds)

source(file.path(source_path,'auc-tjur-plots.R'))


# VP ----------------------------------------------------------------------
VP = computeVariancePartitioning(fitSepTF,start = 1)

# split by groups 
#names <- VP$groupnames
VP_split = computeVariancePartitioning(fitSepTF,start = 1,
                                            group = c(1,1,1,
                                                      2,2,2,
                                                      3,3,
                                                      rep(4,8)),
                                            groupnames = c('temperature',
                                                           'precipitation',
                                                           'landscape',
                                                           'land-use classes'))
# split by seasons
VP_season = computeVariancePartitioning(fitSepTF,start = 1,
                                             group = c(1,2,3,
                                                       1,2,3,
                                                       4,4,
                                                       rep(5,8)),
                                             groupnames = c('year',
                                                            'winter',
                                                            'breeding',
                                                            'landscape',
                                                            'land-use classes'))

source(file.path(source_path,'VP-plots.R'))

# VP BY GUILD / STRATEGY  ------------------------------------------------------------
source(file.path(source_path,'VP-guild-plots.R'))
source(file.path(source_path,'VP-migration-plots.R'))

# VP SORTED BY TJUR  ---------------------------------------------------------------
# get number of species occurrences
occs <- data.frame(occs = colSums(m$Y,na.rm=T),
                   species = colnames(m$Y))

# get species and trait information to identify poorly modelled 
# species and traits 
sp <- fitSepTF$spNames
tr <- fitSepTF$TrData

sp_fits <- data.frame(
  auc = MF$AUC,
  tjur = MF$TjurR2
)
rownames(sp_fits) <- sp

head(sp_fits)
head(tr)

# sanity check that these are the same 
setdiff(rownames(sp_fits),rownames(tr))

# merge traits + occurrences 
semi <- merge(sp_fits,tr,by='row.names')
names(semi)[names(semi)=='Row.names'] <- 'species'
join <- merge(semi,occs,by='species')
head(join)

# get the 10 best & the 10 worst 
sort_asc_tjur <- join[order(join$tjur,decreasing = T),][1:10,]
sort_desc_tjur <- join[order(join$tjur,decreasing = F),][1:10,]

# get plots for these species 









# VP SORTED BY CLASSES  ---------------------------------------------------------------
# get number of species occurrences
VP_split
which(max(VP_split[1,]))
rs <- rownames(VP_split$vals)
for(i in rs){
focal_species <- names(sort(VP_split$vals[i,],decreasing=T)[1:10])
indices <-  which(colnames(VP_split$vals)%in%focal_species)
VP_guild <- VP_split
VP_guild$vals <- VP_split$vals[,indices]
main <- paste0('Variance Partitioning - ',i,' - Grouped by category')
Hmsc::plotVariancePartitioning(fitSepTF,VP_guild,
                               cols = c('firebrick3',
                                        'dodgerblue3',
                                        'goldenrod2',
                                        'springgreen4',
                                        'cornsilk2',
                                        'cornsilk3'),
                               main = main,
                               las = 2,
                               border = NA,
                               space=0,
                               axisnames=T,
                               ann=T,
                               args.legend = list(x = 'topright',
                                                  inset=c(-0.25,0)))
}

# VP SORTED BY MOST IMPORTANT TO A SPECIES  -------------------------------
pdf(file=file.path(input,'results','VP_mostimportantvars.pdf'),
    width = 14,
    height = 6)

par(mar = c(15,5,5,15))

# TEMPERATURE SPECIES 
info <- apply(VP_split$vals,2,which.max)
table(info)

focal_species <- names(info[info==1])
indices <-  which(colnames(VP_split$vals)%in%focal_species)
VP_guild <- VP_split
VP_guild$vals <- VP_split$vals[,indices]
main <- paste0('Variance Partitioning - ','Temperature most important',' - Grouped by category')
Hmsc::plotVariancePartitioning(fitSepTF,VP_guild,
                               cols = c('firebrick3',
                                        'dodgerblue3',
                                        'goldenrod2',
                                        'springgreen4',
                                        'cornsilk2',
                                        'cornsilk3'),
                               main = main,
                               las = 2,
                               border = NA,
                               space=0,
                               axisnames=T,
                               ann=T,
                               args.legend = list(x = 'topright',
                                                  inset=c(-0.25,0)))

# WINTER SPECIES 
info <- apply(VP_season$vals,2,which.max)
table(info)

focal_species <- names(info[info==2])
indices <-  which(colnames(VP_season$vals)%in%focal_species)
VP_guild <- VP_season
VP_guild$vals <- VP_season$vals[,indices]
main <- paste0('Variance Partitioning - ','Winter most important',' - Grouped by season')
Hmsc::plotVariancePartitioning(fitSepTF,VP_guild,
                               cols = c('coral',
                                        'lightblue',
                                        'lightgreen',
                                        'goldenrod2',
                                        'springgreen4',
                                        'cornsilk2',
                                        'cornsilk3'),
                               main = main,
                               las = 2,
                               border = NA,
                               space=0,
                               axisnames=T,
                               ann=T,
                               args.legend = list(x = 'topright',
                                                  inset=c(-0.25,0)))

# YEARLY TEMPERATURE MOST IMPORTANT
info <- apply(VP$vals,2,which.max)
table(info)

focal_species <- names(info[info==1])
indices <-  which(colnames(VP$vals)%in%focal_species)
VP_guild <- VP
VP_guild$vals <- VP$vals[,indices]
main <- paste0('Variance Partitioning - ','Yearly temperature most important')
Hmsc::plotVariancePartitioning(fitSepTF,VP_guild,
                               cols = c('firebrick3',
                                        'firebrick2',
                                        'firebrick1',
                                        'dodgerblue3',
                                        'dodgerblue2',
                                        'dodgerblue1',
                                        'goldenrod2',
                                        'goldenrod1',
                                        colorRampPalette(c("springgreen4", "springgreen1"))(8),
                                        'cornsilk2',
                                        'cornsilk3'),
                               main = main,
                               las = 2,
                               border = NA,
                               space=0,
                               axisnames=T,
                               ann=T,
                               args.legend = list(x = 'topright',
                                                  inset=c(-0.25,0)))

dev.off()


