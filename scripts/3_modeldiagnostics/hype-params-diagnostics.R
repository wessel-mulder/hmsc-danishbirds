# DEFINING INPUT ----------------------------------------------------------
input <- 'tmp_rds/mods-tuning-hypeparam'
dirs <- list.dirs(input)[-1] #drop toplevel

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
# make empty vectors to store results 

for(x in seq_along(dirs)){
  if(x == 1){
    # make empty df to store results 
    df <- data.frame()
  }
  print(paste0('Starting with ',dir,'...'))

  # get 
  dir <- dirs[x]
  # load unfitted object
  m <- readRDS(file.path(dir,'m_object.rds'))
  summary(m)
  # load params 
  params <- readRDS(file.path(dir,'params.rds'))
  
  # loading the chains 
  chainList = vector("list", params$nChains)
  for(cInd in 1:params$nChains){
      chain_file_path = file.path(dir, sprintf("post_chain%.2d_file.rds", cInd-1))
      #print(chain_file_path)
      if(file.exists(chain_file_path)) {
        chainList[[cInd]] = from_json(readRDS(file = chain_file_path)[[1]])[[1]]
        }
      }
  
  filteredList <- chainList[!sapply(chainList, is.null)]
  
  
  fitSepTF = importPosteriorFromHPC(m, filteredList,  params$nSamples, params$thin, params$transient)
  
  print('Posterior import, converting to Coda Object...')
  
  # Convert model output to coda format for MCMC diagnostics
  mpost <- convertToCodaObject(fitSepTF,start=1)
  
  # Get number of species (used for dimensions if needed)
  ns <- ncol(fitSepTF$Y)
  
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
  
  # PLOTTING PARAMETERS 
  params2 <- list(mpost$Beta,mpost$Gamma,mpost$Rho,mpost$V,
              mpost$Alpha[[1]],mpost$Alpha[[2]],
              omega_k_1,omega_k_2)
  names <- c('beta','gamma','rho','V',
             'alpha_space','alpha_time',
             'omega_raw','omega_resid')
  
  print('Working on parameters....')
  diags <- list(psrf = list(), ess = list())
  for(i in seq_along(params2)){
    # BETA - species x environment  
    # get psrfs 
    psrf <- gelman.diag(params2[[i]], multivariate = F)$psrf
    ge <- psrf[,1]
    df <- rbind(df,c(
      params$nSamples,
      params$thin,
      'psrf',
      names[i],
      mean(ge,na.rm=T),
      median(ge,na.rm=T)
    ))
  
    # get ess 
    ess <- effectiveSize(params2[[i]])
    df <- rbind(df,c(
      params$nSamples,
      params$thin,
      'ess',
      names[i],
      mean(ess,na.rm=T),
      median(ess,na.rm=T)
    ))
  }
}

colnames(df) <- c('samples','thin','metric','parameter','val_mean','val_median')
df$samples <- as.numeric(df$samples)
df$thin <- as.numeric(df$thin)
df$metric <- as.factor(df$metric)
df$parameter <- as.factor(df$parameter)
df$val_mean <- as.numeric(df$val_mean)
df$val_median <- as.numeric(df$val_median)

write.csv(df,'tmp_rds/mods-tuning-hypeparam/diagnostics.csv')

# PLOTTING ----------------------------------------------------------------
df <- read.csv('tmp_rds/mods-tuning-hypeparam/diagnostics.csv',
               row.names = 1)
head(df)
df$thin <- factor(df$thin)
#df$samples <- factor(df$samples)

df$metric <- factor(df$metric,
                        levels = c('psrf','ess'))
df$parameter <- factor(df$parameter,
                    levels = c('beta','gamma',
                               'rho','V',
                               'alpha_space','alpha_time',
                               'omega_raw','omega_resid'))
df$norm <- df$val_mean
df$norm[df$metric=='ess'] <- df$norm[df$metric=='ess'] / 
  as.numeric(df$samples[df$metric=='ess'])


pdf(file='tmp_rds/mods-tuning-hypeparam/psrf_ess_thin_samp.pdf',
    width=10,height=8)

### ALL OF THEM 
ggplot() +
  # PSRF on primary y-axis
  geom_line(data = df, aes(x = samples, y = val_mean,colour = thin)) +
  geom_point(data = df, aes(x = samples, y = val_mean,colour=thin)) +
  facet_grid(metric~parameter,scales='free_y',
             space='free_x')+
  theme_minimal()



### ONES THAT LOOK GOOD 
sub <- subset(df,parameter %in% c('beta','gamma','rho',
                                  'V','alpha_time'))
ggplot() +
  # PSRF on primary y-axis
  geom_line(data = sub, aes(x = samples, y = val_mean,colour = thin)) +
  geom_point(data = sub, aes(x = samples, y = val_mean,colour=thin)) +
  facet_grid(metric~parameter,scales='free_y',
             space='free_x')+
  theme_minimal()

### OMEGAS
sub <- subset(df,parameter %in% c('omega_raw','omega_resid'))
ggplot() +
  # PSRF on primary y-axis
  geom_line(data = sub, aes(x = samples, y = val_mean,colour = thin)) +
  geom_point(data = sub, aes(x = samples, y = val_mean,colour=thin)) +
  facet_grid(metric~parameter,scales='free_y',
             space='free_x')+
  theme_minimal()

### ALPHA SPACE
sub <- subset(df,parameter %in% c('alpha_space'))
ggplot() +
  # PSRF on primary y-axis
  geom_line(data = sub, aes(x = samples, y = val_mean,colour = thin)) +
  geom_point(data = sub, aes(x = samples, y = val_mean,colour=thin)) +
  facet_wrap(~metric,nrow=2,scales='free_y')+
  theme_minimal()

dev.off()



head(df)



plots <- list()
count <- 0
for(i in c('beta','gamma','rho','V')){
  for(j in c('psrf','ess')){
    count <- count + 1
    sub <- subset(df,parameter==i)
    sub <- subset(sub,metric==j)
    plot <- ggplot()+
      geom_tile(data=sub,
                aes(x = samples,y=thin,fill=val_mean)) +
      if(j == 'psrf'){
      scale_fill_gradient2(low='forestgreen',
                          high = 'firebrick',
                          midpoint = 1.1,
                          limits=c(1,1.2)) 
      }else{
          scale_fill_gradient2(low='firebrick',
                               high = 'forestgreen',
                               midpoint = 400,
                               limits=c(0,2000))
        }
    
    plots[[count]] <- plot
  }
}

plots[[1]]
plots[[2]]
plots[[3]]
plots[[4]]





# OLD ---------------------------------------------------------------------



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
  vioplot(diags$psrf[[i]],col=cols[i],ylim=c(1,max(diags$psrf[[i]])),
          main=paste0('PSRF: ',full_names[i]))
  abline(h=mean(diags$psrf[[i]]),col='red')
  vioplot(diags$psrf[[i]],col=cols[i],ylim=c(1,1.1),main=paste0('Mean = ',round(mean(diags$psrf[[i]]),2)))
  abline(h=mean(diags$psrf[[i]]),col='red')
  par(mfrow=c(1,1))
  hist(diags$ess[[i]],breaks = seq(0,max(diags$ess[[i]])+100,by=100),
       xlab = 'Number of effective samples',
       main=paste0('ESS: ',full_names[i]))
  legend('topright',
         c(paste0('Proportion of effective samples <100 = ',
                round(length(diags$ess[[i]][diags$ess[[i]] < 100]) / length(diags$ess[[i]]),2)),
           paste0('Proportion of effective samples <200 = ',
                  round(length(diags$ess[[i]][diags$ess[[i]] < 200]) / length(diags$ess[[i]]),2))),
         bty = 'n')

  
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

semi <- merge(sp_fits,tr,by='row.names')
names(semi)[names(semi)=='Row.names'] <- 'species'
join <- merge(semi,occs,by='species')

par(mfrow=c(1,1),
    mar = c(10,4,4,2))

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
mean <- mean(MF$AUC)
hist(MF$AUC,main=paste0('Mean AUC across species = ',floor(mean*100)/100),
     ylab = 'Number of species',
     xlab = 'AUC')
mean <- mean(MF$TjurR2)
hist(MF$TjurR2,main=paste0('Mean TjurR2 across species = ',floor(mean*100)/100),
     ylab = 'Number of species',
     xlab = 'TjurR2')


plot(auc~occs,
       data = join,
     main = 'AUC ~ nr of occurrences',
     xlab = 'Occurrences',
     ylab = 'AUC')
plot(tjur~occs,
     data = join,
     main = 'TjurR2 ~ nr of occurrences',
     xlab = 'Occurrences',
     ylab = 'TjurR2')


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

