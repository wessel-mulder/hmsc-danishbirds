# PLOTTING PARAMETERS 
params <- list(mpost$Beta,mpost$Gamma,mpost$V,mpost$Sigma,
               mpost$Eta[[1]],
               mpost$Alpha[[1]],
               mpost$Omega[[1]],
               mpost$Lambda[[1]],
               mpost$Psi[[1]],
               mpost$Delta[[1]])


# plot names 
full_names <- c("Species x Environment",
                "Traits x Environment",
                "Residual covariance (species niches)",
                "Residual variance (species occurrences)",
                'Site loadings 1',
                'Scale of latent factors 1',
                "Species associations 1",
                "Species loadings 1",
                'Local shrinkage - species loadings 1',
                'Global shrinkage - species loadings 1')

# plot colors 
cols <- c('green4','purple4','orange4','firebrick',
          'cornsilk3','ivory','blue4',
          'lightblue','cyan3','cyan4')

all_psrf <- list()
all_ess  <- list()

for(i in seq_along(full_names)){
  all_psrf[[ full_names[i] ]] <- as.numeric(diags$psrf[[i]])
  all_ess[[ full_names[i] ]]  <- as.numeric(diags$ess[[i]])
}

diags$psrf[[4]]

# init pdf 
pdf(file=file.path(input, "results", "PSRF_ESS_combined.pdf"), width=10, height=6)

par(mfrow=c(1,1))

# PSRF plot
print('psrf plot full')

all_psrf <- lapply(all_psrf, function(x) {
  summary(x)
  x[!is.finite(x)] <- 1
  if (length(unique(x)) == 1) {
    x[2] <- x[1] + 0.001
  }
  x
})

vioplot(all_psrf, col=cols, main="PSRF across parameters",
        ylim=c(0.99, 5),
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

