# PLOTTING PARAMETERS 
params <- list(mpost$Beta,mpost$Gamma,mpost$Rho,mpost$V,mpost$Sigma,
               mpost$Eta[[1]],mpost$Eta[[2]],
               mpost$Alpha[[1]],mpost$Alpha[[2]],
               mpost$Omega[[1]],mpost$Omega[[2]],
               mpost$Lambda[[1]],mpost$Lambda[[2]],
               mpost$Psi[[1]],mpost$Psi[[2]],
               mpost$Delta[[1]],mpost$Delta[[2]])


# plot names 
full_names <- c("Species x Environment",
                "Traits x Environment",
                "Phylogeny",
                "Residual covariance (species niches)",
                'Residual variance (occurrences)',
                'Site loadings 1',
                'Site loadings 2',
                'Scale of latent factors 1',
                'Scale of latent factors 2',
                "Species associations 1",
                "Species associations 2",
                "Species loadings 1",
                "Species loadings 2",
                'Local shrinkage - species loadings 1',
                'Local shrinkage - species loadings 2',
                'Global shrinkage - species loadings 1',
                'Global shrinkage - species loadings 2')
print(full_names)
                
                
                 
                

# plot colors 
cols <- c('green4','purple4','yellow4','orange4','firebrick',
          'lightgreen','lightgreen','ivory','ivory','blue4','blue4',
          'cornsilk2','cornsilk2','cyan1','cyan1','cyan3','cyan3')

all_psrf <- list()
all_ess  <- list()

for(i in seq_along(full_names)){
  print(i)
  print(as.numeric(diags$psrf[[i]]))
  all_psrf[[ full_names[i] ]] <- as.numeric(diags$psrf[[i]])
  all_ess[[ full_names[i] ]]  <- as.numeric(diags$ess[[i]])
}


# init pdf 
pdf(file=file.path(input, "results", "PSRF_ESS_combined.pdf"), width=10, height=6)

par(mfrow=c(1,1))

# PSRF plot
print('psrf plot full')

all_psrf <- lapply(all_psrf, function(x) {
  summary(x)
  x[!is.finite(x)] <- 10
  x[is.na(x)] <- 10
  if(!length(x)){
    x[2] <- x[1] + 0.001
  }
  x
})

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

