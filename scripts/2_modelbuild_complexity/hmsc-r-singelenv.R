rm(list = ls())
n_cores=4

# GETTING STARTED ---------------------------------------------------------
if (interactive() && Sys.getenv("RSTUDIO") == "1") {
  message("Running in RStudio")
  library(Hmsc)
  library(jsonify)
  library(knitr)
  library(corrplot)
  library(ape)
  library(dplyr)
  library(MASS)
  library(vioplot)
  start <- '.'
  python <- file.path(getwd(), 'hmsc-hpc-main',"hmsc-venv", "bin", "python3.11")

} 

# make dirs 
mod <- 'v1_tmean'
input <- file.path('.','tmp_rds','mods-complexity','v1_tmean')
if(!dir.exists(input)) {dir.create(input)}
if(!dir.exists(file.path(input,'model-outputs'))){dir.create(file.path(input,'model-outputs'))}
if(!dir.exists(file.path(input,'results'))){dir.create(file.path(input,'results'))}

# check if accurately installed
summary(TD)


# LOADING BASIC DATA - TEMPERATURE MODEL   ------------------------------------------------------------
### ENVIRONMENT
X <- read.csv(file.path(start,'data/1_preprocessing/X_environmental/X_Environmental.csv'),row.names=1)
X <- X[sort(row.names(X)),]

# keep only sites with data  
X <- na.omit(X)

# keep only atlas 3
atlas3 <- X[rownames(X)[grep("_3$", rownames(X))],,drop=F]

# grab tmean_year
atlas3_temp <- atlas3[,'tmean_year',drop=F]

sites_actual <- row.names(atlas3_temp)

### OCCURRENCES 
Y <- read.csv(file.path(start,'data/1_preprocessing/Y_occurrences/Y_occurrences.csv'),row.names=1)

# remove sites without data 
Y <- Y[row.names(Y) %in% sites_actual,]

# grab 13 warblers
genera <- c('Phylloscopus','Curruca','Sylvia','Acrocephalus')
keep <- sapply(strsplit(colnames(Y),'_'),head,1) %in% genera
Y_warblers <- Y[,keep]

# check if they have enough records 
colSums(Y_warblers, na.rm =T)
# barred warbler is absent so we remove it 
Y_warblers <- Y_warblers[colnames(Y_warblers) != "Curruca_nisoria"]
# fixed
colSums(Y_warblers, na.rm =T)


# --> LOAD STUDY DESIGN 
Design <- read.csv(file.path(start,"data/1_preprocessing/design/studyDesign.csv"),row.names=5)

# remove sites without data 
Design <- Design[row.names(Design) %in% sites_actual,]

# sort
Design <- Design[sort(row.names(Design)), ]
Design3 <- Design[rownames(Design)[grep("_3$", rownames(Design))],,drop=F]
Design3 <- Design3[,c('site','lat','lon')]

head(Design3)

# convert to factors
Design3$site <- as.factor(Design3$site)

# get xycoords
xycoords <- Design3 
rownames(xycoords) <- xycoords$site
xycoords <- xycoords[,colnames(xycoords) %in% c('lat','lon')]

# SANITY CHECKS -----------------------------------------------------------
head(Y_warblers)
### XY COORDS OF BIRDS 
# check with atlas3 to make sure everything is correct
leaf <- data.frame(pa=Y_warblers[,'Phylloscopus_trochiloides',drop=F])
# strip the "_3" from the rownames
rownames(leaf) <- sub("_3$", "", rownames(leaf))
#mergs
pa_leaf <- merge(xycoords,leaf,by='row.names')
# plot 
library(ggplot2)
ggplot(pa_leaf,
       aes(x=lon,y=lat,color=Phylloscopus_trochiloides))+
  geom_point()

### XY DISTRIBUTION OF ENVIRONMENT
tmean <- data.frame(temp=atlas3_temp[,'tmean_year',drop=F])
# strip the "_3" from the rownames
rownames(tmean) <- sub(paste0("_",3,"$"), "", rownames(tmean))
#mergs
tmean_space <- merge(xycoords,tmean,by='row.names')
# plot 
ggplot(tmean_space,
          aes(x=lon,y=lat,color=tmean_year))+
  geom_point()+
  scale_color_gradientn(,colors=c('blue','red'))

# TEMPERATURE MODEL  ------------------------------------------------------

Formula_temp <- as.formula(paste("~", paste(colnames(atlas3_temp), collapse = "+"), sep = " "))
struc_space <- HmscRandomLevel(sData = xycoords,sMethod = 'NNGP',nNeighbours=10)
#struc_space <- HmscRandomLevel(sData = xycoords,sMethod = 'Full',longlat = T)

struc_space$alphapw
# checking out priors
# standard frequencies 
freq <- c(rep(0.01,100))
#under 1 degree
samples <- c(seq(from = 0.07, to = 1, length.out = 100))
small <- cbind(samples,freq)

struc_space_small <- setPriors(struc_space,alphapw=small)
struc_space_small$alphapw

# define m 
m.temp <-Hmsc(Y = Y_warblers, 
         XData = atlas3_temp, 
         XFormula = Formula_temp,
         studyDesign = Design3[,c('site'),drop=F], 
         ranLevels = list('site'=struc_space_small),
         distr='probit')

m.temp$ranLevels[[1]]$alphapw

nChains = 1
test.run = T
if(test.run){
  thin = 1
  samples = 100
  transient = 100
  verbose =1
}else{
  thin = 10
  samples = 100
  transient = 1
  verbose = 10
}

time <- system.time({
m.temp.sampled = sampleMcmc(m.temp, thin = thin, samples = samples, transient = transient,
                       nChains = nChains, nParallel = nChains, verbose = verbose)
})
time
mpost <- convertToCodaObject(m.temp.sampled)
plot(mpost$Alpha[[1]])

summary(mpost$Alpha[[1]])
summary(mpost$Eta[[1]])


# storing mod 
mod <- 'v1_tmean'
input <- file.path('.','tmp_rds','mods-complexity','v1_tmean')
if(!dir.exists(input)) {dir.create(input)}
saveRDS(m.temp.sampled,file.path(input,'mod.rds'))
saveRDS(m.temp,file.path(input,'m.rds'))

# LOOKING AT PSRF & ESS  ---------------------------------------------------------
m.temp.sampled <- readRDS(file.path(input,'mod.rds'))
mpost <- convertToCodaObject(m.temp.sampled)

summary(mpost)

# objects in the list 
params <- list(mpost$Beta,mpost$Gamma,mpost$V,mpost$Sigma,
               mpost$Eta[[1]],
               mpost$Alpha[[1]],
               mpost$Omega[[1]],
               mpost$Lambda[[1]],
               mpost$Psi[[1]],
               mpost$Delta[[1]])

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
saveRDS(diags,file.path(input,'model-outputs','psrf-ess.rds'))


# AND PLOTTING THEM  ------------------------------------------------------
# PLOTTING PARAMETERS 
diags <- readRDS(file.path(file.path(input,'model-outputs','psrf-ess.rds')))

params <- list(mpost$Beta,mpost$Gamma,mpost$V,mpost$Sigma,
               mpost$Eta[[1]],
               mpost$Alpha[[1]],
               mpost$Omega[[1]],
               mpost$Lambda[[1]],
               mpost$Psi[[1]],
               mpost$Delta[[1]])



# plot names 
full_names <- c("Species x Environment",
                "Gamma intercept",
                "Residual covariance in species x environment associations",
                'Residual variance in species occurrences',
                'Site loadings',
                'Decay factor of latent, spatial process',
                'Species associations',
                "Species loadings",
                'Local shrinkage - species loadings',
                'Global shrinkage - species loadings')

cols <- c('green4','purple4','orange4','firebrick',
          'lightgreen','ivory','blue4',
          'cornsilk2','cyan1','cyan3')

all_psrf <- list()
all_ess  <- list()

for(i in seq_along(full_names)){
  print(i)
  print(as.numeric(diags$psrf[[i]]))
  all_psrf[[ full_names[i] ]] <- as.numeric(diags$psrf[[i]])
  all_ess[[ full_names[i] ]]  <- as.numeric(diags$ess[[i]])
}

# process a bit to make plotting visible

all_psrf <- lapply(all_psrf, function(x) {
  summary(x)
  x[!is.finite(x)] <- 10
  x[is.na(x)] <- 10
  if(!length(x)){
    x[2] <- x[1] + 0.001
  } else if (length(unique(x)) == 1) {
    # All values are identical
    x[1] <- x[1] + 0.001
    x[2] <- x[2] - 0.001
  }
  x
})

all_ess <- lapply(all_ess, function(x) {
  summary(x)
  x[!is.finite(x)] <- 1
  x[is.na(x)] <- 1
  if(!length(x)){
    x[2] <- x[1] + -1
  } else if (length(unique(x)) == 1) {
    # All values are identical
    x[1] <- x[1] + 1
    x[2] <- x[2] - 2
  }
  x
})

# init pdf 
pdf(file=file.path(input, "results", "PSRF_ESS_combined.pdf"), width=10, height=6)

par(mfrow=c(1,1))

# PSRF plot
print('psrf plot full')

vioplot(all_psrf, col=cols, main="PSRF across parameters",
        ylim=c(0.99, 11),
        xaxt = 'n',ann =  T)
abline(h=1.1, lty=1, col="red") # warning threshold
#axis(1, at=seq_along(all_psrf), labels=names(all_psrf), las=2, cex.axis=0.7)
legend('topleft',
       fill = cols,
       full_names)
print('psrf plot threshold')
vioplot(all_psrf, col=cols, main="PSRF across parameters",
        ylim=c(0.99, 1.1),
        xaxt = 'n',ann =  T)
abline(h=1.1, lty=1, col="red") # warning threshold
#axis(1, at=seq_along(all_psrf), labels=names(all_psrf), las=2, cex.axis=0.7)


# ESS plot
vioplot(all_ess, col=cols, main="ESS across parameters",
        xaxt = 'n',ann=T)
abline(h=100, lty=1, col="red") # warning threshold
#axis(1, at=seq_along(all_ess), labels=names(all_ess), las=3, cex.axis=0.7)
# ESS plot
vioplot(all_ess, col=cols, main="ESS across parameters",
        ylim=c(0, 1000),
        xaxt = 'n',ann=T)
abline(h=100, lty=1, col="red") # warning threshold
#axis(1, at=seq_along(all_ess), labels=names(all_ess), las=3, cex.axis=0.7)

dev.off()


# LETS LOOK AT TRACE PLOTS ------------------------------------------------
# get indices of worst chains 
pdf(file = file.path(input,'results','chains.pdf'),
    width = 6,
    height = 8)
for(i in 1:10){
  print(i)
  params <- list(mpost$Beta,mpost$Gamma,mpost$V,mpost$Sigma,
                 mpost$Eta[[1]],
                 mpost$Alpha[[1]],
                 mpost$Omega[[1]],
                 mpost$Lambda[[1]],
                 mpost$Psi[[1]],
                 mpost$Delta[[1]])
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

mpost <- convertToCodaObject(m.temp.sampled)
plot(mpost$Alpha[[1]])
# MODEL FIT ---------------------------------------------------------------
preds <- computePredictedValues(m.temp.sampled)
MF <- evaluateModelFit(m.temp.sampled,preds)
saveRDS(MF,file=file.path(input,'model-outputs','model-fit.rds'))

# plot results 
MF <- readRDS(file.path(input,'model-outputs','model-fit.rds'))

# get number of species occurrences
occs <- data.frame(occs = colSums(m.spatial$Y,na.rm=T),
                   species = colnames(m.spatial$Y))

sp_fits <- data.frame(
  auc = MF$AUC,
  tjur = MF$TjurR2,
  rmse = MF$RMSE
)
rownames(sp_fits) <- row.names(occs)
head(occs)
head(sp_fits)

join <- merge(sp_fits,occs,by='row.names')

par(mfrow=c(1,1),
    mar = c(10,4,4,2))

join$species_short <- sub("^([A-Za-z])[A-Za-z]+_([a-z]+)$", "\\1. \\2", join$species)


plot.new()
# make pdf 
pdf(file=file.path(input,'results','AUC_TjurR2.pdf'),
    width = 8,
    height = 6)
par(mfrow=c(1,1))

plot(auc~occs,
     data = join,
     main = 'AUC ~ nr of occurrences',
     xlab = 'Occurrences',
     ylab = 'AUC',
     xlim=c(0,2500))
join$pos <- 2
join$pos[join$occs<1000] <- 4
join$pos[join$species_short%in%c('S. borin','C. communis')] <- 1
join$pos[join$species_short%in%c('P. collybita','C. curruca','P. sibilatrix','A. schoenobaenus','A. scirpaceus')] <- 3
join$pos[join$species_short%in%c('P. collybita','P. trochilus')] <- 4
text(join$occs, join$auc, labels=join$species_short, cex= 1,pos = join$pos)
text(0,0.99,
     paste0('Mean AUC = ',round(mean(join$auc),2)),
     pos=4,
     cex=1.5)

plot(tjur~occs,
     data = join,
     main = 'TjurR2 ~ nr of occurrences',
     xlab = 'Occurrences',
     ylab = 'TjurR2',
     ylim=c(0,0.5),
     xlim=c(0,2500))
join$postjur <- 2
join$postjur[join$occs<1000] <- 4
join$postjur[join$species_short%in%c(,'S. atricapilla')] <- 1
join$postjur[join$species_short%in%c('P. sibilatrix','S. borin','C. curruca','A. schoenobaenus','A. scirpaceus','A. palustris','P. collybita','C. communis')] <- 3

text(join$occs, join$tjur, labels=join$species_short, cex= 1,pos = join$postjur)
text(0,0.5,
     paste0('Mean TjurR2 = ',round(mean(join$tjur),2)),
     pos=4,
     cex=1.5)
dev.off()


# VARIANCE PARTITIONING ---------------------------------------------------
par(mfrow=c(1,1))
m.temp.sampled <- readRDS(file.path(input,'mod.rds'))
VP <- computeVariancePartitioning(m.temp.sampled)

# edit function
plotVariancePartitioning2 =
  function (hM, VP, cols=NULL, main = "Variance Partitioning", ...)
  {
    ng = dim(VP$vals)[1]
    if(is.null(cols)){
      cols = heat.colors(ng, alpha = 1)
    }
    leg = VP$groupnames
    for (r in seq_len(hM$nr)) {
      leg = c(leg, paste("Random: ", hM$rLNames[r], sep = ""))
    }
    means = round(100 * rowMeans(VP$vals), 1)
    for (i in 1:ng) {
      leg[i] = paste(leg[i], " (mean = ", toString(means[i]),
                     ")", sep = "")
    }
    
    barplot(VP$vals, main = main, xlab= " ", ylab = "Total variance", las = 1,
            legend = leg, col = cols,...)
    #   mtext("Species", 1,line = 1)
  }
# and plot
{
pdf(file=file.path(input,'results','VP.pdf'),
    width = 10,
    height = 5)
par(mar = c(10,5,5,2))

VP_2 <- VP
VP_2$groupnames <- c('Yearly temp.')
tjurs <- MF$TjurR2
VP_2$vals <- VP$vals * tjurs
colnames(VP_2$vals) <- sub("^([A-Za-z])[A-Za-z]+_([a-z]+)$", "\\1. \\2", colnames(VP_2$vals))

plotVariancePartitioning2(m.temp.sampled,
                          cols = c('firebrick3','cornsilk2'),
                         VP_2,
                         main = 'Total variance explained',
                         las = 2,
                         border = NA,
                         space=0.01,
                         ann=T,
                         ylim=c(0,1),
                         args.legend = list(x = 'topright'))
                                            #inset=c(-0.25,0)))
dev.off()
}
# PREDICTS ----------------------------------------------------------------
covariate <- c('tmean_year')
#parallel::mclapply(covariates, function(covariate) {

print('starting gradient')
Gradient <- constructGradient(m.temp, focalVariable = covariate, ngrid = 10)
print('starting predictions')

predY <- predict(m.temp.sampled, Gradient = Gradient, expected = TRUE,
                 nParallel = n_cores)

# Save each covariate separately
saveRDS(list(predY = predY, Gradient = Gradient),
        file = file.path(input,'model-outputs', paste0("pred_", covariate, ".rds")))

plotGradient(m.temp,Gradient,predY,measure='S')


