input <- 'tmp_rds/mods-single/2025-08-29_15-52-09_samples_1000_thin_100/'
subset_sppairs <- T
reset <- 1
if(reset){
beta <- NULL
gamma <- NULL
V <- NULL
rho <- NULL
alpha <- NULL
omega_1 <- NULL
omega_2 <- NULL
}


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

filteredList <- chainList[!sapply(chainList, is.null)]


fitSepTF = importPosteriorFromHPC(m, filteredList, nSamples, thin, transient)


# PSRF ----------------------------------------------------
# Convert model output to coda format for MCMC diagnostics
mpost <- convertToCodaObject(fitSepTF,start=1)

# Get number of species (used for dimensions if needed)
ns <- ncol(fitSepTF$Y)

summary(mpost)

### OMEGA PREPROCESSING
omegas = mpost$Omega
maxOmega = 1000
omega_pre_1=omegas[[1]]
omega_pre_2=omegas[[2]]

z = ncol(omega_pre_1[[1]])  # assume same structure for all chains
if (z > maxOmega) {
  sel = sample(1:z, size = maxOmega)
  omega_k_1 = lapply(omega_pre_1, function(chain) chain[, sel, drop=FALSE])
  omega_k_2= lapply(omega_pre_2, function(chain) chain[, sel, drop=FALSE])
} else{
  omega_k_1 <- omega_pre_1
  omega_k_2 <- omega_pre_2
}
gelman.diag(mpost$Omega[[1]])

# PLOTTING PARAMETERS 
params <- list(mpost$Beta,mpost$Gamma,mpost$Rho,mpost$V,
            mpost$Alpha[[1]],mpost$Alpha[[2]],
            omega_k_1,omega_k_2)

diags <- list(psrf = list(), ess = list())
for(i in seq_along(params)){
  # BETA - species x environment  
  print(i)
  # get psrfs 
  psrf <- gelman.diag(params[[i]], multivariate = F)$psrf
  ge <- psrf[,1]
  diags$psrf[[i]] <- ge
  
  # get ess 
  ess <- effectiveSize(params[[i]])
  diags$ess[[i]] <- ess
}

names <- c('beta','gamma','rho','V',
           'alpha_space','alpha_time',
           'omega_raw','omega_resid')
full_names <- c("Beta - Species x Environment",
                "Gamma - Traits x Environment",
                "Rho - Phylogeny",
                "V - Variation",
                "Alpha - Spatial",
                "Alpha - Time",
                "Omega1 - Raw associations",
                "Omega2 - Residual associations")

cols <- c('green4','purple4','yellow4','orange4','ivory','ivory','blue4','blue4')

for(i in seq_along(params)){
  ### PSRFS 
  if(i == 1){
    # init output documents 
    if(!dir.exists(file.path(input,'results'))) {dir.create(file.path(input,'results'))}
    text.file = file.path(input,'results','PSRF_ESS.txt')
    pdf(file=file.path(input,'results','PSRF_ESS.pdf'),
        width = 8,
        height = 6)
    
    cat("MCMC Convergence statistics\n\n",file=text.file,sep="")
    cat(c("\n",input,"\n\n"),file=text.file,sep="",append=TRUE)
    
  }
  
  # concatenate output 
  cat(paste0("\n",'PSRF: ',full_names[i],"\n\n"),file=text.file,sep="",append=TRUE)
  cat(capture.output(summary(diags$psrf[[i]])), file = text.file, sep = "\n", append = TRUE)  
  cat(paste0("\n",'ESS: ',full_names[i],"\n\n"),file=text.file,sep="",append=TRUE)
  cat(capture.output(summary(diags$ess[[i]])), file = text.file, sep = "\n", append = TRUE)  
  
  # plot PSRF
  par(mfrow=c(1,2))
  diags$psrf[[i]][diags$psrf[[i]] > 10] <- 10
  vioplot(diags$psrf[[i]],col=cols[i],ylim=c(1,max(diags$psrf[[i]],na.rm = T)),
          main=paste0('PSRF: ',full_names[i]))
  abline(h=1.1,col='red')
  vioplot(diags$psrf[[i]],col=cols[i],ylim=c(1,1.1),main=paste0('Mean = ',round(mean(diags$psrf[[i]],na.rm = T),2)))
  abline(h=1.1,col='red')
  par(mfrow=c(1,1))
  hist(diags$ess[[i]],breaks = seq(0,max(diags$ess[[i]])+100,by=100),
       xlab = 'Number of effective samples',
       main=paste0('ESS: ',full_names[i]))
  abline(v=100,col='red',lty=1)
  abline(v=200,col='red',lty=2)
  abline(v=1000,col='red',lty=3)
  legend(
    "topright",
    legend = c(
      paste0("Proportion of effective samples >100 = ",
             round(length(diags$ess[[i]][diags$ess[[i]] > 100]) / length(diags$ess[[i]]), 2)),
      paste0("Proportion of effective samples >200 = ",
             round(length(diags$ess[[i]][diags$ess[[i]] > 200]) / length(diags$ess[[i]]), 2)),
      paste0("Proportion of effective samples >1000 = ",
             round(length(diags$ess[[i]][diags$ess[[i]] > 1000]) / length(diags$ess[[i]]), 2))
    ),
    col = "red",
    lty = c(1, 2, 3)
  )

  if(i == length(params)){dev.off()}
  if(i == length(params)){
    # --- 1. Collect all psrf and ess into one structure ---
    all_psrf <- list()
    all_ess  <- list()
    
    for(i in seq_along(params)){
      all_psrf[[ full_names[i] ]] <- as.numeric(diags$psrf[[i]])
      all_ess[[ full_names[i] ]]  <- as.numeric(diags$ess[[i]])
    }
    
    # --- 2. Open PDF for output ---
    pdf(file=file.path(input, "results", "PSRF_ESS_combined.pdf"), width=10, height=6)
    
    # --- 3. Violinplot for PSRF ---
    par(mfrow=c(1,1))  # left = PSRF, right = ESS
    
    # PSRF plot
    vioplot(all_psrf, col=cols, main="PSRF across parameters",
            ylim=c(0.99, max(unlist(all_psrf), na.rm=TRUE)),
            xaxt = 'n',ann =  T)
    abline(h=1.1, lty=1, col="red") # warning threshold
    #axis(1, at=seq_along(all_psrf), labels=names(all_psrf), las=2, cex.axis=0.7)
    legend('topleft',
           fill = cols,
           full_names)
    
    vioplot(all_psrf, col=cols, main="PSRF across parameters",
            ylim=c(0.99, 1.1),
            xaxt = 'n',ann =  T)
    abline(h=1.1, lty=1, col="red") # warning threshold
    #axis(1, at=seq_along(all_psrf), labels=names(all_psrf), las=2, cex.axis=0.7)

    
    # ESS plot
    vioplot(all_ess, col=cols, main="ESS across parameters",
            ylim=c(0, max(unlist(all_ess), na.rm=TRUE)),
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
    
  }
  
}
diags$psrf[[5]]

# EVALUATE MODEL FIT  -----------------------------------------------------
# AUC AND TJUR 
# for(i in c(100,500,1000,5000,50)){
# preds <- computePredictedValues(fitSepTF,thin=i)
# MF <- evaluateModelFit(hM=fitSepTF, predY=preds)
# par(mfrow=c(1,2))
# mean <- mean(MF$AUC)
# hist(MF$AUC,main=paste0(i,' = ',mean))
# mean <- mean(MF$TjurR2)
# hist(MF$TjurR2,main=paste0(i,' = ',mean))
# }
rm(chainList,filteredList)

dev.off()
# get auc and tjur from model fits 
preds_10 <- computePredictedValues(fitSepTF,thin=10)
MF_10 <- evaluateModelFit(hM=fitSepTF, predY=preds_10)
hist(MF_10$AUC,
     main=round(mean(MF_10$AUC),2))
hist(MF_10$TjurR2,
     main=round(mean(MF_10$TjurR2),2))

# get number of species occurrences
occs <- data.frame(occs = colSums(m$Y,na.rm=T),
                   species = colnames(m$Y))


# get species and trait information to identify poorly modelled 
# species and traits 
sp <- fitSepTF$spNames
tr <- fitSepTF$TrData

sp_fits <- data.frame(
  auc = MF_10$AUC,
  tjur = MF_10$TjurR2
)
rownames(sp_fits) <- sp

semi <- merge(sp_fits,tr,by='row.names')
names(semi)[names(semi)=='Row.names'] <- 'species'
join <- merge(semi,occs,by='species')

par(mfrow=c(1,1),
    mar = c(10,4,4,2))

sort_desc_tjur <- join[order(join$tjur,decreasing = F),]


### ORDER BY NR OF SPECIES IN GUILD 
join_foraging <- join %>%
  group_by(foraging_guild_consensus) %>%
  summarise(auc_mu = mean(auc),
            auc_sd = sd(auc),
            tjur_sd = sd(tjur),
            tjur_mu = mean(tjur),
            n = dplyr::n(),
            .groups = "drop") %>%
  ungroup()

join_foraging <- join_foraging[order(join_foraging$n, decreasing = TRUE), ]
tail(join_foraging)

join$foraging_guild_consensus <- factor(join$foraging_guild_consensus,
                                        levels = join_foraging$foraging_guild_consensus)


plot.new()
# make pdf 
{
pdf(file=file.path(input,'results','AUC_TjurR2.pdf'),
    width = 8,
    height = 6)
par(mfrow=c(1,1))
mean <- mean(MF_10$AUC)
hist(MF_10$AUC,main=paste0('Mean AUC across species = ',floor(mean*100)/100),
     ylab = 'Number of species',
     xlab = 'AUC')
mean <- mean(MF_10$TjurR2)
hist(MF_10$TjurR2,main=paste0('Mean TjurR2 across species = ',floor(mean*100)/100),
     ylab = 'Number of species',
     xlab = 'TjurR2')


plot(auc~occs,
       data = join,
     main = 'AUC ~ nr of occurrences',
     xlab = 'Occurrences',
     ylab = 'AUC',
     ylim=c(0.75,1),
     xlim=c(0,6500))

plot(tjur~occs,
     data = join,
     main = 'TjurR2 ~ nr of occurrences',
     xlab = 'Occurrences',
     ylab = 'TjurR2',
     ylim = c(0,0.65),
     xlim=c(0,6500))

boxplot(auc~foraging_guild_consensus,data=join,
        las = 3,                # rotate labels if long
        cex.names = 0.5,       # smaller axis tick labels
        ylab = "Mean AUC",
        xlab = NULL)
boxplot(tjur~foraging_guild_consensus,data=join,
        las = 3,                # rotate labels if long
        cex.names = 0.5,       # smaller axis tick labels
        ylab = "Mean AUC",
        xlab = NULL)


bp <- barplot(
  height = join_foraging$auc_mu,
  names.arg = join_foraging$foraging_guild_consensus,
  ylim = c(0,1.1),
  las = 3,                # rotate labels if long
  cex.names = 0.75,       # smaller axis tick labels
  ylab = "Mean AUC"
)
# Add SD whiskers
bp+arrows(
  x0 = bp, y0 = join_foraging$auc_mu - join_foraging$auc_sd,
  x1 = bp, y1 = join_foraging$auc_mu + join_foraging$auc_sd,
  angle = 90, code = 3, length = 0.05
)


bp <- barplot(
  height = join_foraging$tjur_mu,
  names.arg = join_foraging$foraging_guild_consensus,
  ylim=c(0,0.6),
  las = 3,                # rotate labels if long
  cex.names = 0.75,       # smaller axis tick labels
  ylab = "Mean TjurR2",
)
bp+arrows(
  x0 = bp, y0 = join_foraging$tjur_mu - join_foraging$tjur_sd,
  x1 = bp, y1 = join_foraging$tjur_mu + join_foraging$tjur_sd,
  angle = 90, code = 3, length = 0.05
)



# ### MIGRATORY
# join_migratory <- join %>%
#   group_by(Migration_AVONET) %>%
#   summarise(auc_mu = mean(auc),
#             tjur_mu = mean(tjur),
#             n = dplyr::n(),
#             .groups = "drop") %>%
#   ungroup()
# 
# join_migratory <- join_migratory[order(join_migratory$n, decreasing = TRUE), ]
# tail(join_migratory)
# barplot(
#   height = join_migratory$auc_mu,
#   names.arg = join_migratory$Migration_AVONET,
#   las = 1,                # rotate labels if long
#   ylab = "Mean AUC"
# )
# barplot(
#   height = join_migratory$tjur_mu,
#   #names.arg = join_migratory$Migration_AVONET,
#   #las = 2,                # rotate labels if long
#   ylab = "Mean TjurR2",
# )

dev.off()
}

# PLOTTING ACROSS GRADIENTS -----------------------------------------------
Gradient = constructGradient(fitSepTF,focalVariable = "tmean_winter",
                             ngrid=30,
                             non.focalVariables = list(tmean_breeding=list(1),
                                                       prec_winter=list(1),
                                                       prec_breeding=list(1),
                                                       hh=list(1),
                                                       unique=list(1)))
Gradient$rLNew
post <- poolMcmcChains(filteredList)
pred = predict(fitSepTF,post=post,
               Gradient = Gradient,expected = T,predictEtaMean = T)

for(i in 1:3){
  plotGradient(fitSepTF,Gradient,pred=pred,measure = 'T',
               index=i)
}

length(unique(fitSepTF$TrData$foraging_guild_consensus))
Gradient$XDataNew

predY = predict(fitSepTF, XData=Gradient$XDataNew, studyDesign=Gradient$studyDesignNew,
                ranLevels=Gradient$rLNew, expected=TRUE)
plotGradient(fitSepTF, Gradient, pred=predY, measure="S", showData = TRUE)

preds <- computePredictedValues(fitSepTF,partition.sp = 1,thin=20)



# VARIANCE PARTITIONG FULL ENV -----------------------------------------------------
# by groups 


rm(chainList,filteredList,preds)
VP_1950 = computeVariancePartitioning(fitSepTF,start = 900)
# split by groups 
names <- VP_1950$groupnames
VP_split_1950 = computeVariancePartitioning(fitSepTF,start = 900,
                                            group = c(1,1,1,
                                                      2,2,2,
                                                      3,3,
                                                      rep(4,8)),
                                            groupnames = c('temperature',
                                                           'precipitation',
                                                           'landscape',
                                                           'land-use classes'))
# split by seasons
VP_season_1950 = computeVariancePartitioning(fitSepTF,start = 900,
                                             group = c(1,2,3,
                                                       1,2,3,
                                                       4,4,
                                                       rep(5,8)),
                                             groupnames = c('year',
                                                            'winter',
                                                            'breeding',
                                                            'landscape',
                                                            'land-use classes'))

{
pdf(file=file.path(input,'results','VP_full.pdf'),
    width = 8,
    height = 6)

default_margins <- par("mar")
default_margins
#> [1] 5.1 4.1 4.1 2.1

# Make the bottom margin larger
new_margins <- default_margins + c(0, 0, 0, 0)
par(mar = new_margins)

greens <- colorRampPalette(c("springgreen4", "springgreen1"))(8)
greens

plotVariancePartitioning(fitSepTF,
                         VP_1950,
                         cols = c('firebrick3',
                                  'firebrick2',
                                  'firebrick1',
                                  'dodgerblue3',
                                  'dodgerblue2',
                                  'dodgerblue1',
                                  'goldenrod2',
                                  'goldenrod1',
                                  greens,
                                  'cornsilk2',
                                  'cornsilk3'),
                         las = 2,
                         border = NA,
                         space=0,
                         axisnames=F,
                         legend.text=F,
                         ann=T)
                          
plotVariancePartitioning(fitSepTF,VP_split_1950,
                         cols = c('firebrick3',
                                  'dodgerblue3',
                                  'goldenrod2',
                                  'springgreen4',
                                  'cornsilk2',
                                  'cornsilk3'),
                         las = 2,
                         border = NA,
                         space=0,
                         axisnames=F,
                         legend.text=F,
                         ann=T)

plotVariancePartitioning(fitSepTF,VP_season_1950,
                         cols = c('coral',
                                  'lightblue',
                                  'lightgreen',
                                  'goldenrod2',
                                  'springgreen4',
                                  'cornsilk2',
                                  'cornsilk3'),
                         las = 2,
                         border = NA,
                         space=0,
                         axisnames=F,
                         legend.text=F,
                         ann=T)

dev.off()
}

t <- VP_split$vals

sp <- colnames(t)
t[]


# VARIANCE PARTITIONG SEASONS + LANDSCAPE -----------------------------------------------------
# by groups 


rm(chainList,filteredList,preds)
VP_1950 = computeVariancePartitioning(fitSepTF,start=1950)
# split by groups 
names <- VP_1950$groupnames
VP_split_1950 = computeVariancePartitioning(fitSepTF,start = 1950,
                                            group = c(1,1,
                                                      2,2,
                                                      3,3),
                                            groupnames = c('temperature',
                                                           'precipitation',
                                                           'landscape'))
# split by seasons
VP_season_1950 = computeVariancePartitioning(fitSepTF,start = 1950,
                                             group = c(1,2,
                                                       1,2,
                                                       3,3),
                                             groupnames = c('winter',
                                                            'breeding',
                                                            'landscape'))

{
  pdf(file=file.path(input,'results','VP_full.pdf'),
      width = 8,
      height = 6)
  
  default_margins <- par("mar")
  default_margins
  #> [1] 5.1 4.1 4.1 2.1
  
  # Make the bottom margin larger
  new_margins <- default_margins + c(0, 0, 0, 0)
  par(mar = new_margins)
  
  plotVariancePartitioning(fitSepTF,
                           VP_1950,
                           cols = c('firebrick3',
                                    'firebrick1',
                                    'dodgerblue3',
                                    'dodgerblue1',
                                    'goldenrod2',
                                    'goldenrod1',
                                    'cornsilk2',
                                    'cornsilk3'),
                           las = 2,
                           border = NA,
                           space=0,
                           axisnames=F,
                           #legend.text=F,
                           ann=T)
  
  plotVariancePartitioning(fitSepTF,VP_split_1950,
                           cols = c('firebrick3',
                                    'dodgerblue3',
                                    'goldenrod2',
                                    'cornsilk2',
                                    'cornsilk3'),
                           las = 2,
                           border = NA,
                           space=0,
                           axisnames=F,
                           #legend.text=F,
                           ann=T)
  
  plotVariancePartitioning(fitSepTF,VP_season_1950,
                           cols = c('lightblue',
                                    'lightgreen',
                                    'goldenrod2',
                                    'cornsilk2',
                                    'cornsilk3'),
                           las = 2,
                           border = NA,
                           space=0,
                           axisnames=F,
                           #legend.text=F,
                           ann=T)
  
  dev.off()
}

t <- VP_split$vals

sp <- colnames(t)
t[]


# VARIANCE PARTITIONING PER GUILD  ----------------------------------------
colnames(VP_1900$vals)

# get species and trait information to identify poorly modelled 
# species and traits 
sp <- fitSepTF$spNames
tr <- fitSepTF$TrData
head(tr)

rownames(tr)[tr$foraging_guild_consensus=='Low flycatching feeders']
unique(tr$foraging_guild_consensus)

# GUILDS
for(focal_guild in unique(tr$Migration_AVONET)){
focal_species <- rownames(tr)[tr$Migration_AVONET==focal_guild]
focal_matrices <- VP_split_1950$vals[,focal_species]
greens <- colorRampPalette(c("springgreen4", "springgreen1"))(8)
greens

barplot(focal_matrices,
                         col = c('firebrick3',
                                   'dodgerblue3',
                                   'goldenrod2',
                                   'springgreen4',
                                   'cornsilk2',
                                   'cornsilk3'),
                         las = 2,
                         border = NA,
                         space=0,
                         axisnames=T,
                         legend.text=T,
                         ann=T,
        main = focal_guild)
}

# MIGRATION
for(focal_guild in unique(tr$Migration_AVONET)){
  focal_species <- rownames(tr)[tr$Migration_AVONET==focal_guild]
  focal_matrices <- VP_season_1950$vals[,focal_species]
  greens <- colorRampPalette(c("springgreen4", "springgreen1"))(8)
  greens
  
  barplot(focal_matrices,
          col = c('coral',
                   'lightblue',
                   'lightgreen',
                   'goldenrod2',
                   'springgreen4',
                   'cornsilk2',
                   'cornsilk3'),
          las = 2,
          border = NA,
          space=0,
          axisnames=T,
          legend.text=T,
          ann=T,
          main = focal_guild)
}




join <- merge(sp_fits,tr,by='row.names')
join

head(VP_1900$vals)

# LOOKING AT PARAMETERS ---------------------------------------------------

# how much traits explain their niches
# mostly the different proportional land-use classes 
VP_1950$R2T$Beta

# how much of the traits propagate into explaining the distributions 
# 13% of occurrence is explained by the traits 
VP_1950$R2T$Y

# parameter estimates for environments 
postBeta = getPostEstimate(fitSepTF, parName = "Beta")
postBeta$mean[,1]
postBeta$support[,1]
postBeta$support[,1]

pdf(file=file.path(input, "results", "beta.pdf"), 
    width=10, height=6)

plotBeta(fitSepTF,post = postBeta, 
         param = "Sign", supportLevel = 0.9)
dev.off()

postGamma = getPostEstimate(fitSepTF,parName = 'Gamma')
plotGamma(fitSepTF, post=postGamma, param="Support", supportLevel = 0.95)
#plotGamma(fitSepTF, post=postGamma, param="Mean", supportLevel = 0.95)

postGamma_select <- lapply(postGamma,function(x) x[2:8,])
plotGamma(fitSepTF, post=postGamma, param="Sign", supportLevel = 0.95)

#postOmega = getPostEstimate(fitSepTF,parName = 'Omega')
#plotGamma(fitSepTF, post=postGamma, param="Sign", supportLevel = 0.95)
#plotGamma(fitSepTF, post=postGamma, param="Mean", supportLevel = 0.95)

# residuals
OmegaCor = computeAssociations(fitSepTF)
supportLevel = 0.95
# effect of site 
toPlot = ((OmegaCor[[1]]$support>supportLevel)
          + (OmegaCor[[1]]$support<(1-supportLevel))>0)*OmegaCor[[1]]$mean
corrplot(toPlot, method = "color",
         col=colorRampPalette(c("blue","white","red"))(200),
         tl.cex=.6, tl.col="black",
         title=paste("random effect level:", m$rLNames[1]), mar=c(0,0,1,0))

# effect of atlas 
toPlot = ((OmegaCor[[1]]$support>supportLevel)
          + (OmegaCor[[1]]$support<(1-supportLevel))>0)*OmegaCor[[2]]$mean
corrplot(toPlot, method = "color",
         col=colorRampPalette(c("blue","white","red"))(200),
         tl.cex=.6, tl.col="black",
         title=paste("random effect level:", m$rLNames[2]), mar=c(0,0,1,0))

# 
summary(mpost$Rho)$quant # strong phylogenetic signal 


# PLOTTING ACROSS GRADIENTS -----------------------------------------------
Gradient = constructGradient(fitSepTF,focalVariable = "tmean_year")
predY = predict(fitSepTF,Gradient = Gradient, expected = TRUE,
                post = mpost,
                predictEtaMean = T)




# PLOTTING PARAMETERS -----------------------------------------------------
#' plotPostEstimate Function
#'
#' Generates visualizations for posterior estimates of model parameters from Hmsc models.
#' Supports plotting of "Beta", "Gamma", and "Omega" parameters with options to visualize
#' either the sign or the mean of the estimates. The function can return either a ggplot object
#' or a data frame.
#'
#' @param m An Hmsc model object containing fitted model results.
#' @param parName A character string specifying the parameter to plot. Options are "Beta", "Gamma", or "Omega".
#' @param plotType A character string specifying the type of plot. Options are "Sign" for plotting the sign of the estimates or "Mean" for plotting the mean values.
#' @param returnType A character string specifying the return type. Options are "ggplot" for returning a ggplot object or "dataframe" for returning a data frame.
#' @param supportLevel A numeric value between 0 and 1 specifying the support level to filter the estimates (default is 0.95).
#' @param newEnvNames A character vector of new names for environmental variables (used only if `parName` is "Beta").
#' @param newTraitNames A character vector of new names for traits (used only if `parName` is "Gamma").
#'
#' @details
#' \itemize{
#'   \item{\strong{Beta}:} Plots the environmental covariates vs species. Intercept row is removed if present.
#'   \item{\strong{Gamma}:} Plots species vs traits. Intercept row is removed if present.
#'   \item{\strong{Omega}:} Plots species vs species.
#' }
#'
#' For "Sign" plotType, the function plots the sign of the estimates, and for "Mean", it plots the magnitude with zero values replaced by NA.
#'
#' @return Returns a ggplot object if `returnType` is "ggplot", or a data frame if `returnType` is "dataframe".
#'
#' @examples
#' # Plotting Beta Estimates with Default Names
#' plotPostEstimate(m = TCP_Fit, parName = "Beta", plotType = "Sign", returnType = "ggplot")
#'
#' # Plotting Gamma Estimates with New Trait Names
#' newTraitNames <- c("Trait1", "Trait2", "Trait3")
#' plotPostEstimate(m = TCP_Fit, parName = "Gamma", plotType = "Mean", returnType = "ggplot", newTraitNames = newTraitNames)
#'
#' # Returning Data Frame for Omega Estimates
#' plotPostEstimate(m = TCP_Fit, parName = "Omega", plotType = "Sign", returnType = "dataframe")
#'
#' @export

parName <- 'Beta'
plotType <- 'Sign'
m <- fitSepTF
plotPostEstimate <- function(m,
                             parName = c("Beta", "Gamma", "Omega"),
                             plotType = c("Sign", "Mean"),
                             returnPlot = TRUE,
                             supportLevel = 0.95,
                             newEnvNames = NULL,
                             newTraitNames = NULL,
                             spVector = NULL,
                             main = NULL) {
  # Needs
  require(Hmsc)
  require(tidyverse)
  require(reshape2) 
  
  # Handle Omega parameter case with computeAssociations
  if (parName == "Omega") {
    # Only plots the first ranndom level, can add a parameter here if wanted
    OmegaCor <- computeAssociations(m)
    meanEst <- OmegaCor[[1]]$mean  # Get mean correlation matrix
    supportEst <- OmegaCor[[1]]$support  # Get support matrix
  } else {
    # For Beta, Gamma, use getPostEstimate
    postEstimate <- getPostEstimate(m, parName = parName)
    # subset if spVector is supplied 
    if(length(spVector)){
      postEstimate <- lapply(postEstimate,function(x) {
        x[,spVector]
      }
      )
    }
    meanEst <- postEstimate$mean
    supportEst <- postEstimate$support
  }
  
  # Identify the structure based on the parameter being plotted
  if (parName == "Beta") {
    rowNames <- m$covNames  # Covariates (use new names if provided)
    if(length(spVector)){
      colNames <- spVector
    }else{
      colNames <- m$spNames
    }
    rowLabel <- "Covariates"
    colLabel <- "Species"
    
    # Format the data as sign or mean
    if (plotType == "Sign") {
      toPlot = sign(meanEst)
      toPlot = toPlot * ((supportEst > supportLevel) + (supportEst < (1 - supportLevel)) > 0)
    } else if (plotType == "Mean") {
      toPlot = meanEst
      toPlot = toPlot * ((supportEst > supportLevel) + (supportEst < (1 - supportLevel)) > 0)
      toPlot[toPlot == 0] <- NA  # Replace zeros with NA
    } else {
      stop("Invalid plotType. Choose either 'Sign' or 'Mean'.")
    }
    
    # Remove intercept for Beta (if it's present)
    toPlot <- toPlot[-1, ]  # Remove intercept row
    rowNames <- rowNames[-1]  # Adjust row names
    
    # Remove intercept and set row and column names just before plotting
    if (!is.null(newEnvNames)) {
      rowNames <- newEnvNames  # Use provided new names
    }
    
  } else if (parName == "Gamma") {
    rowNames <- m$covNames   # Species
    colNames <- m$trNames  # Traits (use new names if provided)
    rowLabel <- "Covariates"
    colLabel <- "Traits"
    
    # Format the data as sign or mean
    if (plotType == "Sign") {
      toPlot = sign(meanEst)
      toPlot = toPlot * ((supportEst > supportLevel) + (supportEst < (1 - supportLevel)) > 0)
    } else if (plotType == "Mean") {
      toPlot = meanEst
      toPlot = toPlot * ((supportEst > supportLevel) + (supportEst < (1 - supportLevel)) > 0)
      toPlot[toPlot == 0] <- NA  # Replace zeros with NA
    } else {
      stop("Invalid plotType. Choose either 'Sign' or 'Mean'.")
    }
    
    # Remove intercept for Gamma (if it's present)
    toPlot <- toPlot[-1, ]  # Remove intercept row
    rowNames <- rowNames[-1]  # Adjust row names
    
    # Set names if needed
    if (!is.null(newTraitNames)) {
      colNames <- newTraitNames  # Use provided new names
    }
    
    # Set names if needed
    if (!is.null(newEnvNames)) {
      colNames <- newTraitNames  # Use provided new names
    }
    
  } else if (parName == "Omega") {
    rowNames <- m$spNames   # Species (rows and columns)
    colNames <- m$spNames
    rowLabel <- "Species (Rows)"
    colLabel <- "Species (Columns)"
    
    # Format the data as sign or mean
    if (plotType == "Sign") {
      toPlot = sign(meanEst)
      toPlot = toPlot * ((supportEst > supportLevel) + (supportEst < (1 - supportLevel)) > 0)
    } else if (plotType == "Mean") {
      toPlot = meanEst
      toPlot = toPlot * ((supportEst > supportLevel) + (supportEst < (1 - supportLevel)) > 0)
      toPlot[toPlot == 0] <- NA  # Replace zeros with NA
    } else {
      stop("Invalid plotType. Choose either 'Sign' or 'Mean'.")
    }
    
    # Remove intercept and set row and column names just before plotting
    if (!is.null(newEnvNames)) {
      rowNames <- newEnvNames  # Use provided new names
      colNames <- newEnvNames  # Use provided new names
    }
    
  } else {
    stop("Invalid parName. Choose 'Beta', 'Gamma', or 'Omega'.")
  }
  
  # Format the data as matrix and add column and row names
  plotMat = matrix(toPlot, nrow = length(rowNames), ncol = length(colNames))
  colnames(plotMat) <- colNames
  rownames(plotMat) <- rowNames
  
  # Reorder matrix based on parameter type
  if (parName != "Omega") {
    plotMat <- plotMat[, rev(order(colnames(plotMat)))]  # Reorder columns
  } else {
    plotMat <- plotMat[rev(order(rownames(plotMat))), rev(order(colnames(plotMat)))]  # Reorder rows and columns
  }
  
  # Return data frame if requested
  if (isFALSE(returnPlot)) {
    return(plotMat)
  }
  
  # Prepare data for ggplot
  plotMatmelt <- as.data.frame(melt(as.matrix(plotMat)))
  
  # Define color scheme
  colors = colorRampPalette(c('#4269D0FF', '#fdfdf8', '#fdfdf8', '#EFB118FF'))
  colorLevels <- 3
  cols_to_use = colors(colorLevels)
  
  # Plot based on Sign or Mean
  if (plotType == "Sign") {
    plot <- ggplot(plotMatmelt, aes(x = Var1, y = Var2, fill = factor(value))) +
      labs(x = rowLabel, y = colLabel, fill = "Sign",title = main) +
      geom_tile(color = 'gray60') +
      scale_fill_manual(breaks = levels(factor(plotMatmelt$value)),
                        values = cols_to_use,
                        labels = c('\u2013', '0', '+')) +
      theme(plot.background = element_rect(fill = "white"),
            legend.background = element_rect(fill = "white"),
            panel.border = element_rect(fill = NA, color = NA),
            legend.margin = margin(l = 1, unit = 'cm'),
            legend.title = element_text(hjust = 0.1),
            legend.key.width = unit(1, 'cm'),
            legend.key.height = unit(4, 'cm'),
            legend.text = element_text(size = 12),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      scale_x_discrete(expand = c(0, 0)) +

      scale_y_discrete(expand = c(0, 0))
  } else if (plotType == "Mean") {
    plot <- ggplot(plotMatmelt, aes(x = Var1, y = Var2, fill = value)) +
      labs(x = rowLabel, y = colLabel, fill = "Mean",title = main) +
      geom_tile(color = 'gray60') +
      scale_fill_steps2(n.breaks = 7,
                        high = cols_to_use[3],
                        mid = 'white',
                        low = cols_to_use[1],
                        nice.breaks = TRUE,
                        na.value = "white") +
      theme(plot.background = element_rect(fill = "white"),
            legend.background = element_rect(fill = "white"),
            panel.border = element_rect(fill = NA, color = NA),
            legend.margin = margin(l = 1, unit = 'cm'),
            legend.title = element_text(hjust = 0.08),
            legend.key.width = unit(1, 'cm'),
            legend.key.height = unit(2.25, 'cm'),
            legend.text = element_text(size = 12),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0))
  }
  
  return(plot)
}

rownames(tr)[tr$foraging_guild_consensus==focal_guild]
unique(tr$foraging_guild_consensus)

# GUILDS
for(focal_guild in unique(tr$foraging_guild_consensus)){
  focal_species <- rownames(tr)[tr$foraging_guild_consensus==focal_guild]
  focal_matrices <- lapply(postBeta,function(x) {
    x[,focal_species]
  }
)
}

plotPostEstimate(fitSepTF, parName = "Omega", plotType = "Sign")

plotBeta(fitSepTF,postBeta,
         param='Sign',supportLevel=0.95)

plotGamma(fitSepTF,postGamma,
         param='Sign',supportLevel=0.95)

tr <- fitSepTF$TrData

### MIGRATION SPECIFIC 
focal_migrate <- 'sedentary'
focal_species <- rownames(tr)[tr$Migration_AVONET==focal_migrate]
plotPostEstimate(fitSepTF, parName = "Beta", plotType = "Sign",
                 spVector = focal_species)
plotPostEstimate(fitSepTF, parName = "Gamma", plotType = "Sign",
                 supportLevel = 0.95)

tree <- fitSepTF$phyloTree
plot(tree,
     cex = 1)


