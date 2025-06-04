<<<<<<< HEAD
input <- 'tmp_rds/2025-05-22_16-39-28'

# GETTING STARTED ---------------------------------------------------------
#library(RColorBrewer,lib="~/Rlibs")
#library(farver,lib="~/Rlibs")
#library(scales,lib="~/Rlibs")
#library(jsonify,lib="~/Rlibs")
#library(ape,lib="~/Rlibs")
#library(dplyr,lib="~/Rlibs")
#library(Hmsc,lib="~/Rlibs")
library(Hmsc)
library(jsonify)
library(knitr)
library(corrplot)
=======
input <- '/home/bhr597/home/projects/hmsc-danishbirds/tmp_rds/2025-05-22_16-39-28'

# GETTING STARTED ---------------------------------------------------------
library(RColorBrewer,lib="~/Rlibs")
library(farver,lib="~/Rlibs")
library(scales,lib="~/Rlibs")
library(jsonify,lib="~/Rlibs")
library(ape,lib="~/Rlibs")
library(dplyr,lib="~/Rlibs")
library(Hmsc,lib="~/Rlibs")

>>>>>>> 2d629c8d347170c0dea1f64a425d2a60dc98a8a5
# GETTING STARTED --------------------------------------------------------
# load unfitted object
m <- readRDS(file.path(input,'m_object.rds'))
params <- readRDS(file.path(input,'params.rds'))

<<<<<<< HEAD
nChains <- params$nChains
=======

>>>>>>> 2d629c8d347170c0dea1f64a425d2a60dc98a8a5
# loading the chains 
chainList = vector("list", nChains)
for(cInd in 1:nChains){
    chain_file_path = file.path(input, sprintf("post_chain%.2d_file.rds", cInd-1))
    chainList[[cInd]] = from_json(readRDS(file = chain_file_path)[[1]])[[1]]
}

<<<<<<< HEAD
# import fitted model 
fitSepTF = importPosteriorFromHPC(m, chainList, 100, 10, 1000)


# CHECKING CONVERGENCE ----------------------------------------------------
# Convert model output to coda format for MCMC diagnostics
mpost <- convertToCodaObject(fitSepTF)

# Get number of species (used for dimensions if needed)
ns <- ncol(fitSepTF$Y)

# Effective sample size and Gelman-Rubin diagnostics for all key parameters
# Beta: species responses to environment
es.beta <- effectiveSize(mpost$Beta)
ge.beta <- gelman.diag(mpost$Beta, multivariate = F)$psrf

# Gamma: trait effects
es.gamma <- effectiveSize(mpost$Gamma)
ge.gamma <- gelman.diag(mpost$Gamma, multivariate = F)$psrf

# Rho: phylogenetic signal
es.rho <- effectiveSize(mpost$Rho)
ge.rho <- gelman.diag(mpost$Rho, multivariate = F)$psrf

# V: latent variables
es.V <- effectiveSize(mpost$V)
ge.V <- gelman.diag(mpost$V, multivariate = F)$psrf

# Alpha: spatial structure (assumes only 1 spatial level)
es.alpha <- effectiveSize(mpost$Alpha[[1]])
ge.alpha <- gelman.diag(mpost$Alpha[[1]], multivariate = F)$psrf

# Omega: residual species associations
# (subsets to the first 50 species to limit memory load)
sppairs = matrix(sample(x = 1:ns^2, size = 50))
tmp = mpost$Omega[[1]]
for (chain in 1:length(tmp)){
  tmp[[chain]] = tmp[[chain]][,sppairs]
}

es.omega <- effectiveSize(mpost$temp)
ge.omega <- gelman.diag(mpost$temp, multivariate = F)$psrf

# Combine diagnostics into one list for easy access or export
mixing <- list(
  es.beta = es.beta, ge.beta = ge.beta,
  es.gamma = es.gamma, ge.gamma = ge.gamma,
  es.rho = es.rho, ge.rho = ge.rho,
  es.V = es.V, ge.V = ge.V,
  es.alpha = es.alpha, ge.alpha = ge.alpha,
  es.omega = es.omega, ge.omega = ge.omega
)

# ---- Optional: Trace plot setup ----
# Visual inspection (just a few Beta params as example)
par(mfrow = c(2, 2))
traceplot(mpost$Rho)


# PLOTTING ACROSS GRADIENTS -----------------------------------------------
Gradient = constructGradient(fitSepTF,focalVariable = "tmean_year")
Gradient$XDataNew

predY = predict(fitSepTF, XData=Gradient$XDataNew, studyDesign=Gradient$studyDesignNew,
                ranLevels=Gradient$rLNew, expected=TRUE)
plotGradient(fitSepTF, Gradient, pred=predY, measure="S", showData = TRUE)



# LOOKING AT PARAMETERS ---------------------------------------------------
VP_cluster <- computeVariancePartitioning(fitSepTF,
                                  group = c(rep(1,6),rep(2,10)),
                                  groupnames = c('climate','landuse'))


VP <- computeVariancePartitioning(fitSepTF)

plotVariancePartitioning(fitSepTF, VP, args.legend=list(x="bottomright"))

# how much traits explain their variation
kable(VP$R2T$Beta)

# how much of the traits propagate into explaining the distributions 
VP$R2T$Y

# parameter estimates for environments 
postBeta = getPostEstimate(fitSepTF, parName = "Beta")
plotBeta(m, post = postBeta, param = "Support",
         plotTree = TRUE, supportLevel = 0.95, split=.4, spNamesNumbers = c(F,F))

postGamma = getPostEstimate(fitSepTF,parName = 'Gamma')
plotGamma(fitSepTF, post=postGamma, param="Support", supportLevel = 0.5)

# get predictions 
preds <- computePredictedValues(fitSepTF)
MF <- evaluateModelFit(hM=fitSepTF, predY=preds)
hist(MF$TjurR2, xlim = c(0,1), main=paste0("Mean = ", round(mean(MF$TjurR2),2)))

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
summary(mpost$Rho) # strong phylogenetic signal 


# SUBSETTING BY GUILDS ----------------------------------------------------
Y <- m$Y
Tr <- m$Tr

migrate <- Tr[,2]
guild <- Tr[,c(3:ncol(Tr))]

library(tidyverse)
guild_df <- guild %>%
  as.data.frame() %>%
  rownames_to_column(var = "species") %>%
  pivot_longer(-species, names_to = "guild", values_to = "membership") %>%
  filter(membership == 1) %>%
  mutate(guild = str_remove(guild, "^foraging_guild_consensus")) %>%
  select(species, guild)

guild_df$guild <- as.factor(guild_df$guild)


=======
fitSepTF = importPosteriorFromHPC(m, chainList, 100, 10, 1000)


>>>>>>> 2d629c8d347170c0dea1f64a425d2a60dc98a8a5
