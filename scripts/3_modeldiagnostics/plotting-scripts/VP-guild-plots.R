VPs <- list(VP,VP_split,VP_season)
names <- c('VP_per_guild.pdf',
           'VP_split_per_guild.pdf',
           'VP_season_per_guild.pdf')

for(i in seq_along(VPs)){
  VP_cur <- VPs[[i]]
  pdf(file=file.path(input,'results',names[[i]]),
      width = 14,
      height = 6)
  
  par(mar = c(15,5,5,15))
  
  Tr <- m$TrData
  guilds <- unique(Tr$foraging_guild_consensus)
  
  #guilds <- c('Passerine seedeaters','Diurnal raptors')
  for(g in guilds){
    focal_guild <- g
    focal_species <- rownames(Tr[Tr$foraging_guild_consensus%in%g,])
    n_species <- length(focal_species)
    indices <-  which(colnames(VP_cur$vals)%in%focal_species)
    
    VP_guild <- VP_cur
    VP_guild$vals <- VP_cur$vals[,indices]
    if(n_species == 1){
      VP_guild$vals <- cbind(VP_guild$vals,NA)
    }
    
    if(i == 1){
    main <- paste0('Variance Partitioning - ',focal_guild)
    Hmsc::plotVariancePartitioning(m,VP_guild,
                                   cols = c('firebrick3',
                                            'firebrick2',
                                            'firebrick1',
                                            'dodgerblue3',
                                            'dodgerblue2',
                                            'dodgerblue1',
                                            'goldenrod2',
                                            'goldenrod1',
                                            colorRampPalette(c("springgreen4", "springgreen1"))(8),
                                            'cornsilk2',
                                            'cornsilk3'),
                                   main = main,
                                   las = 2,
                                   border = NA,
                                   space=0,
                                   axisnames=T,
                                   ann=T,
                                   args.legend = list(x = 'topright',
                                                      inset=c(-0.25,0)))
    }else if(i == 2){
      main <- paste0('Variance Partitioning - ',focal_guild,' - Grouped by category')
      Hmsc::plotVariancePartitioning(m,VP_guild,
                                     cols = c('firebrick3',
                                              'dodgerblue3',
                                              'goldenrod2',
                                              'springgreen4',
                                              'cornsilk2',
                                              'cornsilk3'),
                                     main = main,
                                     las = 2,
                                     border = NA,
                                     space=0,
                                     axisnames=T,
                                     ann=T,
                                     args.legend = list(x = 'topright',
                                                        inset=c(-0.25,0)))
      

    }else if(i == 3){
      main <- paste0('Variance Partitioning - ',focal_guild,' - Grouped by season')
      Hmsc::plotVariancePartitioning(m,VP_guild,
                                     cols = c('coral',
                                              'lightblue',
                                              'lightgreen',
                                              'goldenrod2',
                                              'springgreen4',
                                              'cornsilk2',
                                              'cornsilk3'),
                                     main = main,
                                     las = 2,
                                     border = NA,
                                     space=0,
                                     axisnames=T,
                                     ann=T,
                                     args.legend = list(x = 'topright',
                                                        inset=c(-0.25,0)))
    }
    
    
  }
  dev.off()
}

