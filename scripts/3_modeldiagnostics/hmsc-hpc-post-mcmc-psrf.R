input <- 'tmp_rds/2025-07-08_12-02-50_gpp/'
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

### OMEGA PREPROCESSING
omegas = mpost$Omega
maxOmega = 100
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

# PLOTTING PARAMETERS 
params <- list(mpost$Beta,mpost$Gamma,mpost$Rho,mpost$V,
            mpost$Alpha[[1]],mpost$Alpha[[2]],
            omega_k_1,omega_k_2)
names <- c('beta','gamma','rho','V',
           'alpha_space','alpha_time',
           'omega_space','omega_time')
full_names <- c("PSRF: Beta - Species x Environment",
                "PSRF: Gamma - Traits x Environment",
                "PSRF: Rho - Phylogeny",
                "PSRF: V - autocorrelation?",
                "PSRF: Alpha - Spatial",
                "PSRF: Alpha - Time",
                "PSRF: Omega1 - Spec. Assoc. Space",
                "PSRF: Omega2 - Spec. Assoc. Time")

cols <- c('green4','purple4','yellow4','orange4','ivory','ivory','blue4','blue4')

for(i in seq_along(params)){
  if(i == 1){
    # init output documents 
    if(!dir.exists(file.path(input,'results'))) {dir.create(file.path(input,'results'))}
    text.file = file.path(input,'results','MCMC_convergence.txt')
    pdf(file=file.path(input,'results','MCMC_convergence.pdf'),
        width = 8,
        height = 6)
    
    cat("MCMC Convergence statistics\n\n",file=text.file,sep="")
    cat(c("\n",input,"\n\n"),file=text.file,sep="",append=TRUE)
    

  }
  
  # BETA - species x environment  
  beta <- gelman.diag(params[[i]], multivariate = F)$psrf
  ge.beta <- beta[,1]
  # concatenate output 
  cat(paste0("\n",full_names[i],"\n\n"),file=text.file,sep="",append=TRUE)
  cat(capture.output(summary(ge.beta)), file = text.file, sep = "\n", append = TRUE)  # plot 
  par(mfrow=c(1,2))
  vioplot(ge.beta,col=cols[i],ylim=c(0,max(ge.beta)),main=full_names[i])
  vioplot(ge.beta,col=cols[i],ylim=c(0.9,1.1),main="")
  
  if(i == length(params)){dev.off()}
}


# EFFECTIVE SAMPLE SIZES  -------------------------------------------------
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

# get auc and tjur from model fits 
preds <- computePredictedValues(fitSepTF,thin=50)
MF <- evaluateModelFit(hM=fitSepTF, predY=preds)

# make pdf 
{
pdf(file=file.path(input,'results','AUC_TjurR2.pdf'),
    width = 8,
    height = 6)
par(mfrow=c(1,2))
mean <- mean(MF$AUC)
hist(MF$AUC,main=paste0('AUC = ',floor(mean*100)/100))
mean <- mean(MF$TjurR2)
hist(MF$TjurR2,main=paste0('TjurR2 = ',floor(mean*100)/100))

# get species and trait information to identify poorly modelled 
# species and traits 
sp <- fitSepTF$spNames
tr <- fitSepTF$TrData

sp_fits <- data.frame(
  auc = MF$AUC,
  tjur = MF$TjurR2
)
rownames(sp_fits) <- sp

join <- merge(sp_fits,tr,by='row.names')
join

par(mfrow=c(2,1))

### FORAGING GUILD 
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
barplot(
  height = join_foraging$auc_mu,
  names.arg = join_foraging$foraging_guild_consensus,
  las = 3,                # rotate labels if long
  cex.names = 0.75,       # smaller axis tick labels
  ylab = "Mean AUC"
)
barplot(
  height = join_foraging$tjur_mu,
  #names.arg = join_foraging$foraging_guild_consensus,
  #las = 3,                # rotate labels if long
  #cex.names = 0.75,       # smaller axis tick labels
  ylab = "Mean TjurR2",
)

### MIGRATORY
join_migratory <- join %>%
  group_by(Migration_AVONET) %>%
  summarise(auc_mu = mean(auc),
            tjur_mu = mean(tjur),
            n = dplyr::n(),
            .groups = "drop") %>%
  ungroup()

join_migratory <- join_migratory[order(join_migratory$n, decreasing = TRUE), ]
tail(join_migratory)
barplot(
  height = join_migratory$auc_mu,
  names.arg = join_migratory$Migration_AVONET,
  las = 1,                # rotate labels if long
  ylab = "Mean AUC"
)
barplot(
  height = join_migratory$tjur_mu,
  #names.arg = join_migratory$Migration_AVONET,
  #las = 2,                # rotate labels if long
  ylab = "Mean TjurR2",
)

dev.off()
}




# PLOTTING ACROSS GRADIENTS -----------------------------------------------
Gradient = constructGradient(fitSepTF,focalVariable = "tmean_year")
pred = predict(fitSepTF,Gradient = Gradient)

Gradient$XDataNew

predY = predict(fitSepTF, XData=Gradient$XDataNew, studyDesign=Gradient$studyDesignNew,
                ranLevels=Gradient$rLNew, expected=TRUE)
plotGradient(fitSepTF, Gradient, pred=predY, measure="S", showData = TRUE)

preds <- computePredictedValues(fitSepTF,partition.sp = 1,thin=20)



# VARIANCE PARTITIONG -----------------------------------------------------
# by groups 
VP = computeVariancePartitioning(fitSepTF,start = 1900)
#saveRDS(VP,file.path(input,'results','VP'))
#VP <- readRDS(file.path(input,'results','VP'))
plotVariancePartitioning(fitSepTF,VP)

VP$group

# split by groups 
VP_split = computeVariancePartitioning(fitSepTF,start = 1900,
                                 group = c(1,1,1,
                                           2,2,2,
                                           3,3,
                                           4,4,4,4,4,4,4,4),
                                 groupnames = c('temperature',
                                                'precipitation',
                                                'landscape',
                                                'land-use classes'))
plotVariancePartitioning(fitSepTF,VP_split,
                         cols = c('red2',
                                  'blue2',
                                  'yellow2',
                                  'green2',
                                  'orange3',
                                  'orange4'))

t <- VP_split$vals

sp <- colnames(t)
t[]

# LOOKING AT PARAMETERS ---------------------------------------------------

# how much traits explain their niches
# mostly the different proportional land-use classes 
VP$R2T$Beta

# how much of the traits propagate into explaining the distributions 
# 13% of occurrent is explained by the traits 
VP$R2T$Y

# parameter estimates for environments 
postBeta = getPostEstimate(fitSepTF, parName = "Beta")
postBeta$mean[,1]
postBeta$support[,1]
postBeta$support[,1]




plotBeta(m, post = postBeta, param = "Sign",
         plotTree = F, supportLevel = 0.95, split=.4, spNamesNumbers = c(F,F))
plotBeta(m, post = postBeta, param = "Mean",
         plotTree = F, supportLevel = 0.95, split=.4, spNamesNumbers = c(F,F),
         colorLevels = c(100))

postGamma = getPostEstimate(fitSepTF,parName = 'Gamma')
plotGamma(fitSepTF, post=postGamma, param="Sign", supportLevel = 0.95)
plotGamma(fitSepTF, post=postGamma, param="Mean", supportLevel = 0.95)

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



fitSepTF = importPosteriorFromHPC(m, chainList, 100, 10, 1000)



input <- '/home/bhr597/home/projects/hmsc-danishbirds/tmp_rds/2025-05-22_16-39-28'




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
plotPostEstimate <- function(m,
                             parName = c("Beta", "Gamma", "Omega"),
                             plotType = c("Sign", "Mean"),
                             returnPlot = TRUE,
                             supportLevel = 0.95,
                             newEnvNames = NULL,
                             newTraitNames = NULL) {
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
    meanEst <- postEstimate$mean
    supportEst <- postEstimate$support
  }
  
  # Identify the structure based on the parameter being plotted
  if (parName == "Beta") {
    rowNames <- m$covNames  # Covariates (use new names if provided)
    colNames <- m$spNames   # Species
    rowLabel <- "Environmental Niche"
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
      labs(x = rowLabel, y = colLabel, fill = "Sign") +
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
      labs(x = rowLabel, y = colLabel, fill = "Mean") +
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


plotPostEstimate(fitSepTF, parName = "Beta", plotType = "Sign")
plotPostEstimate(fitSepTF, parName = "Beta", plotType = "Mean")
plotPostEstimate(fitSepTF, parName = "Gamma", plotType = "Sign")
plotPostEstimate(fitSepTF, parName = "Gamma", plotType = "Mean")
plotPostEstimate(fitSepTF, parName = "Omega", plotType = "Sign")
plotPostEstimate(fitSepTF, parName = "Omega", plotType = "Mean")

