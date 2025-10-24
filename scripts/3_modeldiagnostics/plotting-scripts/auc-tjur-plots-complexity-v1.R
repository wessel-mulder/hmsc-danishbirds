# edit xlims for all species 
if(all_atlas==1){xlimits <- 6000}else if(all_atlas==0)(xlimits <- 2000)

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


# SPECIES  ----------------------------------------------------------------

plot.new()
# make pdf 
pdf(file=file.path(input,'results','AUC_TjurR2.pdf'),
    width = 8,
    height = 6)
par(mfrow=c(1,1))

plot(auc~occs,
     data = join,
     main = paste0('Mean AUC = ', round(mean(join$auc), 2)),
     xlab = 'Number of ccurrences',
     ylab = 'AUC',
     xlim=c(0,xlimits),
     ylim=c(0.5,1))
# # position, 1234 bltr
# join$postauc <- 2
# join$postauc[join$occs<1000] <- 4
# # --- adjust for species that are "too close" ---
# threshold <- 60  # occurrences distance threshold
# 
# # Sort by occs to make it easier
# ord <- order(join$occs)
# occs_sorted <- join$occs[ord]
# 
# # Sort by occs to make it easier
# ord <- order(join$occs)
# occs_sorted <- join$occs[ord]
# 
# # Loop through consecutive species
# for (i in seq_along(occs_sorted)[-1]) {
#   diff <- occs_sorted[i] - occs_sorted[i-1]
#   if (diff <= threshold) {
#     # Find their original indices
#     idx1 <- ord[i-1]
#     idx2 <- ord[i]
#     
#     # If only two within threshold, assign positions 1 and 3
#     if(occs_sorted[i] > 500){
#     join$postauc[idx1] <- 1
#     join$postauc[idx2] <- 3
#     }
#     
#     # Check if there's a third one also close
#     if (i > 1 && i < length(occs_sorted)) {
#       if ((occs_sorted[i+1] - occs_sorted[i]) <= threshold) {
#         idx3 <- ord[i+1]
#         join$postauc[idx1] <- 1
#         join$postauc[idx2] <- 2
#         join$postauc[idx3] <- 3
#       }
#     }
#   }
# }
# # Add species labels with adjusted positions
# text(join$occs, join$auc, labels = join$species_short,
#      cex = 1, pos = join$postauc)
# Add mean line
# text(0,0.99,
#      paste0('Mean AUC = ', round(mean(join$auc), 2)),
#      pos = 4, cex = 1.5)

## AND TJUR 
plot(tjur~occs,
     data = join,
     main = paste0('Mean TjurR2 = ', round(mean(join$tjur), 2)),
     xlab = 'Number of occurrences',
     ylab = 'TjurR2',
     ylim=c(0,0.6),
     xlim=c(0,xlimits))
# join$postjur <- 2
# join$postjur[join$occs<1000] <- 4
# # --- adjust for species that are "too close" ---
# threshold <- 60  # occurrences distance threshold
# 
# # Sort by occs to make it easier
# ord <- order(join$occs)
# occs_sorted <- join$occs[ord]

# # Loop through consecutive species
# for (i in seq_along(occs_sorted)[-1]) {
#   diff <- occs_sorted[i] - occs_sorted[i-1]
#   if (diff <= threshold) {
#     # Find their original indices
#     idx1 <- ord[i-1]
#     idx2 <- ord[i]
#     
#     # If only two within threshold, assign positions 1 and 3
#     # if their occs are below 1000
#     if(occs_sorted[i] > 500){
#     join$postjur[idx1] <- 1
#     join$postjur[idx2] <- 3
#     }
#     # Check if there's a third one also close
#     if (i > 1 && i < length(occs_sorted)) {
#       if ((occs_sorted[i+1] - occs_sorted[i]) <= threshold) {
#         idx3 <- ord[i+1]
#         join$postjur[idx1] <- 1
#         join$postjur[idx2] <- 2
#         join$postjur[idx3] <- 3
#       }
#     }
#   }
# }
# # Add species labels with adjusted positions
# text(join$occs, join$tjur, labels = join$species_short,
#      cex = 1, pos = join$postjur)
# # Add mean line
# text(0,0.6,
#       paste0('Mean TjurR2 = ', round(mean(join$tjur), 2)),
#       pos = 4, cex = 1.5)

dev.off()


# GUILD / STRATEGIES VIOLIN PLOTS  ----------------------------------------

head(join)
# make pdf 
pdf(file=file.path(input,'results','AUC_TjurR2_guilds-strategies.pdf'),
    width = 8,
    height = 10)
par(mar=c(3,15,3,3))

traits <- fitSepTF$TrData
rownames(join) <- join$Row.names
join_traits <- merge(join, traits, by = "row.names", all.x = TRUE)
rownames(join_traits) <- join_traits$Row.names
join_traits$Row.names <- NULL
head(join_traits)
specieschar <- 'foraging_guild_consensus'
for(specieschar in c('foraging_guild_consensus','Migration_a3_DOF')){
  for(speciesstat in c('tjur','auc')){
    
    if(specieschar == 'foraging_guild_consensus'){
      namechar <- 'guild'
    }else{
      namechar <- 'migratory strategy'
    }
    
    if(speciesstat == 'auc'){
      namestat <- 'AUC'
      xlim <- c(0.7,1)
    }else{
      namestat <- 'TjurR2'
      xlim <- c(0,0.6)
    }
    
    # Split your data into groups
    groups_start <- split(join_traits[[speciesstat]], join_traits[[specieschar]])
    
    # Compute median (or single value if only one observation)
    group_medians <- sapply(groups_start, function(x) {
      if (length(x) == 1) {
        x  # just return that single value
      } else {
        median(x, na.rm = TRUE)
      }
    })
    
    if(specieschar == 'foraging_guild_consensus' && speciesstat == 'tjur'){
    order_guild_tjur <- names(sort(group_medians))
    order <- order_guild_tjur
    }
    if(specieschar == 'foraging_guild_consensus' && speciesstat == 'auc'){
    order <- order_guild_tjur
    }
    if(specieschar == 'Migration_a3_DOF'){
    order <- c('long-distance',
               'short-and long-distance',
               'short-distance',
               'sedentary and short-distance',
               'sedentary')
    order <- intersect(order,names(groups_start))
    }

    groups <- groups_start[order]

    # Separate single-value and multi-value groups
    single_groups <- sapply(groups, length) == 1
    multi_groups  <- !single_groups
    
    #join_traits$foraging_guild_consensus <- as.factor(join_traits$foraging_guild_consensus)
    # Make the violin plot

    plot(0:1,0:1,type = 'n',ylim=c(0.5,length(groups)+0.5),xlim=xlim,
         axes=F,ann=F)
    vioplot(groups,add=T,horizontal=T,
            col=c('white'),colMed='grey20',rectCol='grey50',
            lineCol='grey50')
    
    # Add single-value points
    singles <- groups[single_groups]
    points(
      x = unlist(singles),
      y = which(single_groups),
      pch = 19,
      col = 'grey20',
      cex = 1.2
    )
    
    if(speciesstat == 'auc'){axis(side=1,at=seq(xlim[1],xlim[2],0.1))}else{axis(side=1,at=seq(xlim[1],xlim[2],0.2))}
    axis(side=2,at=1:length(groups),labels=names(groups),las=1)
    title(main = paste0(namestat,' by ',namechar))
  }
}
dev.off()

