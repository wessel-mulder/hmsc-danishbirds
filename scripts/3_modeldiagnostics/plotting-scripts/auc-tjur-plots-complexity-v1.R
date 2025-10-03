# get number of species occurrences
occs <- data.frame(occs = colSums(fitSepTF$Y,na.rm=T),
                   species = colnames(fitSepTF$Y))

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
# position, 1234 bltr
join$postauc <- 2
join$postauc[join$occs<1000] <- 4
# --- adjust for species that are "too close" ---
threshold <- 60  # occurrences distance threshold

# Sort by occs to make it easier
ord <- order(join$occs)
occs_sorted <- join$occs[ord]

# Sort by occs to make it easier
ord <- order(join$occs)
occs_sorted <- join$occs[ord]

# Loop through consecutive species
for (i in seq_along(occs_sorted)[-1]) {
  diff <- occs_sorted[i] - occs_sorted[i-1]
  if (diff <= threshold) {
    # Find their original indices
    idx1 <- ord[i-1]
    idx2 <- ord[i]
    
    # If only two within threshold, assign positions 1 and 3
    if(occs_sorted[i] > 500){
    join$postauc[idx1] <- 1
    join$postauc[idx2] <- 3
    }
    
    # Check if there's a third one also close
    if (i > 1 && i < length(occs_sorted)) {
      if ((occs_sorted[i+1] - occs_sorted[i]) <= threshold) {
        idx3 <- ord[i+1]
        join$postauc[idx1] <- 1
        join$postauc[idx2] <- 2
        join$postauc[idx3] <- 3
      }
    }
  }
}
# Add species labels with adjusted positions
text(join$occs, join$auc, labels = join$species_short,
     cex = 1, pos = join$postauc)
# Add mean line
text(0,0.99,
     paste0('Mean AUC = ', round(mean(join$auc), 2)),
     pos = 4, cex = 1.5)

## AND TJUR 
plot(tjur~occs,
     data = join,
     main = 'TjurR2 ~ nr of occurrences',
     xlab = 'Occurrences',
     ylab = 'TjurR2',
     ylim=c(0,0.5),
     xlim=c(0,2500))
join$postjur <- 2
join$postjur[join$occs<1000] <- 4
# --- adjust for species that are "too close" ---
threshold <- 60  # occurrences distance threshold

# Sort by occs to make it easier
ord <- order(join$occs)
occs_sorted <- join$occs[ord]

# Loop through consecutive species
for (i in seq_along(occs_sorted)[-1]) {
  diff <- occs_sorted[i] - occs_sorted[i-1]
  if (diff <= threshold) {
    # Find their original indices
    idx1 <- ord[i-1]
    idx2 <- ord[i]
    
    # If only two within threshold, assign positions 1 and 3
    # if their occs are below 1000
    if(occs_sorted[i] > 500){
    join$postjur[idx1] <- 1
    join$postjur[idx2] <- 3
    }
    # Check if there's a third one also close
    if (i > 1 && i < length(occs_sorted)) {
      if ((occs_sorted[i+1] - occs_sorted[i]) <= threshold) {
        idx3 <- ord[i+1]
        join$postjur[idx1] <- 1
        join$postjur[idx2] <- 2
        join$postjur[idx3] <- 3
      }
    }
  }
}
# Add species labels with adjusted positions
text(join$occs, join$tjur, labels = join$species_short,
     cex = 1, pos = join$postjur)
# Add mean line
text(0,0.5,
     paste0('Mean TjurR2 = ', round(mean(join$tjur), 2)),
     pos = 4, cex = 1.5)

dev.off()

