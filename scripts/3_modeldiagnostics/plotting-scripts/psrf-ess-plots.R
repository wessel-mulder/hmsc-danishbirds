# PLOTTING PARAMETERS 
# official names 
names <- c('beta','gamma','rho','V',
           'alpha_space','alpha_time',
           'omega_raw','omega_resid')

# plot names 
full_names <- c("Beta - Species x Environment",
                "Gamma - Traits x Environment",
                "Rho - Phylogeny",
                "V - Variation",
                "Alpha - Spatial",
                "Alpha - Time",
                "Omega1 - Raw associations",
                "Omega2 - Residual associations")

# plot colors 
cols <- c('green4','purple4','yellow4','orange4','ivory','ivory','blue4','blue4')

all_psrf <- list()
all_ess  <- list()

for(i in seq_along(params)){
  all_psrf[[ full_names[i] ]] <- as.numeric(diags$psrf[[i]])
  all_ess[[ full_names[i] ]]  <- as.numeric(diags$ess[[i]])
}

# init pdf 
if(!dir.exists(file.path(input,'results'))) {dir.create(file.path(input,'results'))}
pdf(file=file.path(input, "results", "PSRF_ESS_combined.pdf"), width=10, height=6)

par(mfrow=c(1,1))

# PSRF plot
print('psrf plot full')

vioplot(all_psrf, col=cols, main="PSRF across parameters",
        ylim=c(0.99, max(unlist(all_psrf), na.rm=TRUE)),
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

