# get number of species occurrences
occs <- data.frame(occs = colSums(m$Y,na.rm=T),
                   species = colnames(m$Y))


# get species and trait information to identify poorly modelled 
# species and traits 
sp <- fitSepTF$spNames
tr <- fitSepTF$TrData

sp_fits <- data.frame(
  auc = MF_10$AUC,
  tjur = MF_10$TjurR2
)
rownames(sp_fits) <- sp

semi <- merge(sp_fits,tr,by='row.names')
names(semi)[names(semi)=='Row.names'] <- 'species'
join <- merge(semi,occs,by='species')

par(mfrow=c(1,1),
    mar = c(10,4,4,2))

sort_desc_tjur <- join[order(join$tjur,decreasing = F),]


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
  mean <- mean(MF_10$AUC)
  hist(MF_10$AUC,main=paste0('Mean AUC across species = ',floor(mean*100)/100),
       ylab = 'Number of species',
       xlab = 'AUC')
  mean <- mean(MF_10$TjurR2)
  hist(MF_10$TjurR2,main=paste0('Mean TjurR2 across species = ',floor(mean*100)/100),
       ylab = 'Number of species',
       xlab = 'TjurR2')
  
  
  plot(auc~occs,
       data = join,
       main = 'AUC ~ nr of occurrences',
       xlab = 'Occurrences',
       ylab = 'AUC',
       ylim=c(0.75,1),
       xlim=c(0,6500))
  
  plot(tjur~occs,
       data = join,
       main = 'TjurR2 ~ nr of occurrences',
       xlab = 'Occurrences',
       ylab = 'TjurR2',
       ylim = c(0,0.65),
       xlim=c(0,6500))
  
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
