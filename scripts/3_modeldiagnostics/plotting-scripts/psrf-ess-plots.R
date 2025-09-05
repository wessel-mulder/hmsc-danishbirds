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

