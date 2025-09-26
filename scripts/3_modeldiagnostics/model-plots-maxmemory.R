args <- commandArgs(trailingOnly = TRUE)

mod = args[1]

# other flags
psrfess_flag <- args[2]
fit_flag <- args[3]
VP_flag <- args[4]
pred_flag <- args[5]
sp_pred_flag <- args[6]

# GETTING STARTED ---------------------------------------------------------
if (interactive() && Sys.getenv("RSTUDIO") == "1") {
  message("Running in RStudio")
  #rstudio = easy_mode
  library(Hmsc)
  library(jsonify)
  library(knitr)
  library(corrplot)
  library(colorspace)
  library(vioplot)
  library(dplyr)
  #mod <- '2025-09-08_17-32-13_samples_1000_thin_100'
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

# load model j
m <- readRDS(file.path(input,'m_object.rds'))
if(!dir.exists(file.path(input,'results'))) {dir.create(file.path(input,'results'))}

# LOADING DATA --------------------------------------------------------
print('starting psrf-ess plots')
if(psrfess_flag == 1){
diags <- readRDS(file.path(input,'model-outputs','psrf-ess.rds'))
source(file.path(source_path,'psrf-ess-plots.R'))
source(file.path(source_path,'psrf-ess-singles.R'))
}
# AUC / TJUR --------------------------------------------------------------
# get preds 
print('starting fit-tjur plots')
if(fit_flag == 1){

MF <- readRDS(file.path(input,'model-outputs','model-fit.rds'))
source(file.path(source_path,'auc-tjur-plots.R'))
}

# VP ----------------------------------------------------------------------
print('starting VP plots')
if(VP_flag==1){
VP <- readRDS(file.path(input,'model-outputs','VP-full.rds'))
VP_split <- readRDS(file.path(input,'model-outputs','VP-split.rds'))
VP_season <- readRDS(file.path(input,'model-outputs','VP-season.rds'))
source(file.path(source_path,'VP-plots.R'))

# VP BY GUILD / STRATEGY  ------------------------------------------------------------
print('starting VP guild/migration plots')

source(file.path(source_path,'VP-guild-plots.R'))
source(file.path(source_path,'VP-migration-plots.R'))

# VP SORTED BY CLASSES  ---------------------------------------------------------------
# get number of species occurrences
# sorted by classes 
print('starting VP other plots')

rs <- rownames(VP_split$vals)
for(i in rs){
focal_species <- names(sort(VP_split$vals[i,],decreasing=T)[1:10])
indices <-  which(colnames(VP_split$vals)%in%focal_species)
VP_guild <- VP_split
VP_guild$vals <- VP_split$vals[,indices]
main <- paste0('Variance Partitioning - ',i,' - Grouped by category')
Hmsc::plotVariancePartitioning(m,VP_guild,
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
# most important vars 

print('starting VP omosti mportant VAR plots')

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
Hmsc::plotVariancePartitioning(m,VP_guild,
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

print('temp succesfull')
# WINTER SPECIES 
info <- apply(VP_season$vals,2,which.max)
table(info)

focal_species <- names(info[info==2])
indices <-  which(colnames(VP_season$vals)%in%focal_species)
VP_guild <- VP_season
VP_guild$vals <- VP_season$vals[,indices]
main <- paste0('Variance Partitioning - ','Winter most important',' - Grouped by season')
Hmsc::plotVariancePartitioning(m,VP_guild,
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

print('winter succesfull')

# YEARLY TEMPERATURE MOST IMPORTANT
info <- apply(VP$vals,2,which.max)
table(info)

focal_species <- names(info[info==1])
indices <-  which(colnames(VP$vals)%in%focal_species)
VP_guild <- VP
VP_guild$vals <- VP$vals[,indices]
main <- paste0('Variance Partitioning - ','Yearly temperature most important')
Hmsc::plotVariancePartitioning(m,VP_guild,
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

print('yearly succesfull')



dev.off()


}


# PARAMETER ESTIMATES -----------------------------------------------------
if(pred_flag==1){
postBeta <- readRDS(file.path(input,'model-outputs','posterior-Beta.rds'))
plotBeta(m,post = postBeta, 
         param = "Sign", supportLevel = 0.95)
postGamma <- readRDS(file.path(input,'model-outputs','posterior-Gamma.rds'))
plotGamma(m,post = postGamma, 
         param = "Sign", supportLevel = 0.95)

getPostEstimate(fitSepTF)
}

