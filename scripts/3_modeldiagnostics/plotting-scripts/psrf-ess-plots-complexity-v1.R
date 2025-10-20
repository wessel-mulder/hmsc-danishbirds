# =============================
# PLOTTING PARAMETERS 
# =============================
megalist <- list(
  beta = list(
    obj   = mpost$Beta,
    name  = "Species x Environment",
    color = "green4"
  ),
  gamma = list(
    obj   = mpost$Gamma,
    name  = "Guilds/strategies x Environment",
    color = "purple4"
  ),
  rho = if (taxonomy_flag == 1) list(
    obj   = mpost$Rho,
    name  = "Taxonomy",
    color = "orange"
  ) else NULL,
  V = list(
    obj   = mpost$V,
    name  = "Residual covariance (species niches)",
    color = "firebrick"
  ),
  sigma = list(
    obj   = mpost$Sigma,
    name  = "Residual variance (species occurrences)",
    color = "firebrick1"
  ),
  eta = list(
    obj   = mpost$Eta[[1]],
    name  = "Site loadings 1",
    color = "cornsilk3"
  ),
  alpha = list(
    obj   = mpost$Alpha[[1]],
    name  = "Scale of latent factors 1",
    color = "ivory"
  ),
  omega = list(
    obj   = mpost$Omega[[1]],
    name  = "Species associations 1",
    color = "blue4"
  ),
  lambda = list(
    obj   = mpost$Lambda[[1]],
    name  = "Species loadings 1",
    color = "lightblue"
  ),
  psi = list(
    obj   = mpost$Psi[[1]],
    name  = "Local shrinkage - species loadings 1",
    color = "cyan3"
  ),
  delta = list(
    obj   = mpost$Delta[[1]],
    name  = "Global shrinkage - species loadings 1",
    color = "cyan4"
  )
)

# Drop NULL entries (so rho disappears cleanly if taxonomy_flag == 0)
megalist <- megalist[!sapply(megalist, is.null)]

# =============================
# GATHER DIAGNOSTICS
# =============================
all_psrf <- list()
all_ess  <- list()

for (i in names(diags$psrf)) {
  all_psrf[[i]] <- as.numeric(diags$psrf[[i]])
  all_ess[[i]]  <- as.numeric(diags$ess[[i]])
}

# Keep only diagnostics matching megalist parameters
all_psrf <- all_psrf[names(all_psrf) %in% names(megalist)]
all_ess  <- all_ess[names(all_ess)  %in% names(megalist)]

# =============================
# INIT PDF 
# =============================
pdf(file=file.path(input, "results", "PSRF_ESS_combined.pdf"), width=10, height=6)

par(mfrow=c(1,1))

# -----------------------------
# PSRF plots
# -----------------------------
cat("Plotting PSRF distributions...\n")

all_psrf <- lapply(all_psrf, function(x) {
  x[!is.finite(x)] <- 1
  if (length(unique(x)) == 1) x[2] <- x[1] + 0.001
  x
})

# Extract colors and names from megalist
cols <- sapply(megalist[names(all_psrf)], function(x) x$color)
names_for_plot <- sapply(megalist[names(all_psrf)], function(x) x$name)

# Full range PSRF plot
vioplot(all_psrf, col=cols, main="PSRF across parameters",
        ylim=c(0.99, 5), xaxt='n', ann=TRUE)
abline(h=1.1, lty=1, col="red")
legend('topleft', fill=cols, legend=names_for_plot, cex=0.8)

# Zoomed PSRF plot (<=1.1)
vioplot(all_psrf, col=cols, main="PSRF across parameters (zoomed)",
        ylim=c(0.99, 1.1), xaxt='n', ann=TRUE)
abline(h=1.1, lty=1, col="red")

# -----------------------------
# ESS plots
# -----------------------------
cat("Plotting ESS distributions...\n")

vioplot(all_ess, col=cols, main="ESS across parameters",
        xaxt='n', ann=TRUE)
abline(h=100, lty=1, col="red")

# Zoomed ESS
vioplot(all_ess, col=cols, main="ESS across parameters (zoomed)",
        ylim=c(0, 1000), xaxt='n', ann=TRUE)
abline(h=100, lty=1, col="red")

dev.off()
cat("✅ PSRF and ESS plots saved.\n")