rm(list = ls())
args <- commandArgs(trailingOnly = TRUE)

between <- 'mods-complexity-v1'
dirs <- list.dirs(file.path('./tmp_rds',between),recursive=F)
inaloop <- F
for(dir in seq_along(dirs)){
  if(!dir %in% c(2,8)){
  inaloop <- T

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
  library(fields)
  library(RColorBrewer)
  library(reshape2)
  library(ggplot2)
  

  mod <- '2025-09-26_16-15-15_singleev_tmean_year'
  between <- 'mods-complexity-v1'
  input <- file.path('./tmp_rds',between,mod)
  if(inaloop){
  input <- dirs[dir]
  print(input)
  }
  source_path <- file.path('./scripts/3_modeldiagnostics/plotting-scripts')
  
} else {
  message("Running from terminal or non-interactive environment")
  rstudio = 0
  library(RColorBrewer,lib="~/Rlibs")
  library(fields,lib="~/Rlibs")
  library(farver,lib="~/Rlibs")
  library(scales,lib="~/Rlibs")
  library(jsonify,lib="~/Rlibs")
  library(ape,lib="~/Rlibs")
  library(dplyr,lib="~/Rlibs")
  library(Hmsc,lib="~/Rlibs")
  library(colorspace,lib="~/Rlibs")
  library(vioplot,lib="~/Rlibs")
  library(corrplot,lib='~/Rlibs')
  library(reshape2,lib='~/Rlibs')
  library(ggplot2,lib='~/Rlibs')
  
  
  input <- file.path('~/home/projects/hmsc-danishbirds/tmp_rds',mod)
  source_path <- file.path('~/home/projects/hmsc-danishbirds/scripts/3_modeldiagnostics/plotting-scripts')
  
}

# load model j
m <- readRDS(file.path(input,'m_object.rds'))
if(!dir.exists(file.path(input,'results'))) {dir.create(file.path(input,'results'))}

# LOADING MODEL  ----------------------------------------------------------
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
mpost <-convertToCodaObject(fitSepTF)
print('model succesfully loaded')

# LOADING DATA --------------------------------------------------------
print('starting psrf-ess plots')

diags <- readRDS(file.path(input,'model-outputs','psrf-ess.rds'))

if(between == 'mods-complexity-v1'){
  source(file.path(source_path,'psrf-ess-plots-complexity-v1.R'))
}else{
  source(file.path(source_path,'psrf-ess-plots.R'))
  source(file.path(source_path,'psrf-ess-singles.R'))
}
# AUC / TJUR --------------------------------------------------------------
# get preds 
print('starting fit-tjur plots')
MF <- readRDS(file.path(input,'model-outputs','model-fit.rds'))
if(between == 'mods-complexity-v1'){
  source(file.path(source_path,'auc-tjur-plots-complexity-v1.R'))
}else{
  source(file.path(source_path,'auc-tjur-plots.R'))
}


# VP ----------------------------------------------------------------------
print('starting VP plots')

VP <- readRDS(file.path(input,'model-outputs','VP.rds'))
if(between == 'mods-complexity-v1'){
  source(file.path(source_path,'VP-plots-complexity-v1.R'))
}else{
#VP_split <- readRDS(file.path(input,'model-outputs','VP-split.rds'))
#VP_season <- readRDS(file.path(input,'model-outputs','VP-season.rds'))
source(file.path(source_path,'VP-plots.R'))
}


# SPATIAL PREDICTIONS  ----------------------------------------------------
print('starting spatial preds ')
preds <- readRDS(file.path(input,'model-outputs','pred-vals.rds'))
print(head(preds))
if(between == 'mods-complexity-v1'){
  source(file.path(source_path,'spatial-preds-complexity-v1.R'))
  source(file.path(source_path,'spatial-preds-species-complexity-v1.R'))
  
}else{

}

# POSTERIOR ESTIMATES  ----------------------------------------------------
print('starting post estimates')
postBeta <- readRDS(file.path(input,'model-outputs','posterior-Beta.rds'))
if(between == 'mods-complexity-v1'){
  print('start omega')
  source(file.path(source_path,'posterior-omega-corrplot-complexity-v1.R'))
  source(file.path(source_path,'posterior-beta-complexity-v1.R'))
  
}else{

}

# CHAINS ------------------------------------------------------------------


# VP BY GUILD / STRATEGY  ------------------------------------------------------------
print('starting VP guild/migration plots')
if(between == 'mods-complexity-v1'){
}else{
  source(file.path(source_path,'VP-guild-plots.R'))
  source(file.path(source_path,'VP-migration-plots.R'))
}


# VP SORTED BY CLASSES  ---------------------------------------------------------------
if(between == 'mods-complexity-v1'){
}else{
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


dev.new()
}

# PARAMETER ESTIMATES -----------------------------------------------------
# plot.new()
# postBeta <- readRDS(file.path(input,'model-outputs','posterior-Beta.rds'))
# plotBeta(m,post = postBeta,
#          param = "Sign", supportLevel = 0.95)
# postGamma <- readRDS(file.path(input,'model-outputs','posterior-Gamma.rds'))
# plotGamma(m,post = postGamma,
#          param = "Mean", supportLevel = 0.95)
# 
# getPostEstimate(fitSepTF)
# 
# postAlpha <- getPostEstimate(fitSepTF,parName='Alpha',r=2)
# mpost<-convertToCodaObject(fitSepTF)
# 
# # spatial structure
# summary(mpost$Alpha[[1]])[2]$quantiles[1,]
# # temporal structure
# summary(mpost$Alpha[[2]])[2]$quantiles[1,]
# 
# summary(mpost$Alpha[[2]])
# 
# fitSepTF
# 
# postEta <- getPostEstimate(fitSepTF,parName='Eta')
# summary(postEta$support)
# 
# 
# i <- 1
# for(i in 1:ncol(postEta$mean)){
#   eta <- postEta$mean[,i]
#   eta2 <- cbind(eta,xy)
#   p<-ggplot(eta2,aes(x=lon,y=lat,col=eta))+
#     geom_point(cex=3)+
#     labs(title=i)+
#     scale_colour_gradient2(
#       low = "blue",      # color for negative values
#       mid = "ivory",     # color for zero
#       high = "red",      # color for positive values
#       midpoint = 0
#     )
#   print(p)
# }
# 
# xy <- fitSepTF$rL[[1]]$s
# 
# summary(mpost$Alpha[[1]])[2]
# 
# 
# 
# eta2_df <- cbind(eta2,xy)
# ggplot(eta2_df,aes(x=lon,y=lat,col=eta2))+
#   geom_point(cex=3)+
#   scale_colour_gradient2(
#     low = "blue",      # color for negative values
#     mid = "ivory",     # color for zero
#     high = "red",      # color for positive values
#     midpoint = 0
#   )
# 
# eta3_df <- cbind(eta3,xy)
# ggplot(eta3_df,aes(x=lon,y=lat,col=eta3))+
#   geom_point(cex=3)+
#   scale_colour_gradient2(
#     low = "blue",      # color for negative values
#     mid = "ivory",     # color for zero
#     high = "red",      # color for positive values
#     midpoint = 0
#   )
# 
# eta4_df <- cbind(eta4,xy)
# ggplot(eta4_df,aes(x=lon,y=lat,col=eta4))+
#   geom_point(cex=3)+
#   scale_colour_gradient2(
#     low = "blue",      # color for negative values
#     mid = "ivory",     # color for zero
#     high = "red",      # color for positive values
#     midpoint = 0
#   )
# 
# eta5_df <- cbind(eta5,xy)
# ggplot(eta5_df,aes(x=lon,y=lat,col=eta5))+
#   geom_point(cex=3)+
#   scale_colour_gradient2(
#     low = "blue",      # color for negative values
#     mid = "ivory",     # color for zero
#     high = "red",      # color for positive values
#     midpoint = 0
#   )
# 
# eta6_df <- cbind(eta6,xy)
# ggplot(eta6_df,aes(x=lon,y=lat,col=eta6))+
#   geom_point(cex=3)+
#   scale_colour_gradient2(
#     low = "blue",      # color for negative values
#     mid = "ivory",     # color for zero
#     high = "red",      # color for positive values
#     midpoint = 0
#   )


# TRACE PLOTS -------------------------------------------------------------
if(between == 'mods-complexity-v1'){
}else{
par(mfrow=c(1,1))
library(rethinking)

mpost$Alpha[[1]]

mpost<-convertToCodaObject(fitSepTF)
summary(mpost)


plot(mpost$Alpha[[1]])
plot(mpost$Alpha[[2]])

plot(mpost$Alpha[[1]],auto.layout=F)

plot(mpost$Alpha[[1]][,1:4])

traceplot(mpost$Alpha[[1]])
plot(mpost$Beta[,1:10],auto.layout=F)
dev.new()
}
  }
}

