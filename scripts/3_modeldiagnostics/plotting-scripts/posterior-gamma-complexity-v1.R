
# PLOTTING FUNCTION -------------------------------------------------------
m <- fitSepTF
if(!is.null(fitSepTF$TrData)){
parName = 'Gamma'
plotType = 'Sign'
characteristic = c('Migratory strategy / Guild')
newEnvNames <- c('Yearly temp.','Yearly prec.','Heterogeneity')

plotGammaSign <- function(m,
                             parName = c("Beta", "Gamma", "Omega"),
                              characteristic = c('Guild','Migratory'),
                             plotType = c("Sign", "Mean"),
                             returnPlot = TRUE,
                             supportLevel = 0.95,
                             newEnvNames = NULL,
                             newTraitNames = NULL) {

  # Identify the structure based on the parameter being plotted

    # For Beta, Gamma, use getPostEstimate
    postEstimate <- getPostEstimate(m, parName = parName)
    meanEst <- postEstimate$mean
    supportEst <- postEstimate$support

    rowNames <- m$covNames  # Covariates (use new names if provided)

    colNames <- abbreviate_genus(m$spNames)   # Species
    
    rowLabel <- "Environmental variable"
    colLabel <- characteristic
    
    # Format the data as sign or mean
    toPlot = sign(meanEst)
    toPlot = toPlot * ((supportEst > supportLevel) + (supportEst < (1 - supportLevel)) > 0)

    # Remove intercept for Beta (if it's present)
    toPlot <- toPlot[-1, ]  # Remove intercept row
    rowNames <- rowNames[-1]  # Adjust row names
    
    # get colnames 
    colnames <- m$trN
    
    # get columns sorted:
    categories <- colnames(m$TrData)
    cat <- 'Migration_a3_DOF'
    uniques_in_cat <- sapply(categories,function(cat){
      uniques <- unique(m$TrData[[cat]])
    })
    all_uniques <- unlist(uniques_in_cat, use.names = FALSE)    # Keep gsub ugly name out
    colnames <- gsub('Migration_a3_DOF',"",colnames)
    colnames <- gsub('foraging_guild_consensus',"",colnames)


      
      # find the intercept
      intercept <- setdiff(all_uniques,colnames)
      if(length(intercept)==1){
        colnames[1] <- paste0('Intercept - ',intercept)
      } else if(length(intercept)==2){
        colnames[1] <- paste0('Intercept - ',intercept[1])
        colnames[1] <- paste0(intercept[1],' - ',intercept[2])
        
      }

    #   print('test')
    #   for(cat in categories){
    #     uniques_in_cat <- unique(m$TrData[[cat]])
    #     
    #     # Keep gsub ugly name out
    #     colnames <- gsub(cat,"",colnames)
    #     
    #     # find the intercept
    #     intercept <- setdiff(uniques_in_cat,colnames)
    #     
    #     colnames[1] <- paste0('Intercept - ',intercept)
    #   }
    # }
    # 
    # Remove intercept and set row and column names just before plotting
    if (!is.null(newEnvNames)) {
      rowNames <- newEnvNames  # Use provided new names
    }
  
  colNames <- colnames
    
  # Format the data as matrix and add column and row names
  plotMat = matrix(toPlot, nrow = length(rowNames), ncol = length(colNames))
  plotMat
  
  colnames(plotMat) <- colNames
  rownames(plotMat) <- rowNames
  
  # Reorder matrix based on parameter type
  #plotMat <- plotMat[, rev(order(colnames(plotMat))),drop=F]  # Reorder columns
  
  # Prepare data for ggplot
  plotMatmelt <- as.data.frame(melt(as.matrix(plotMat)))
  
  # Define color scheme
  colors = colorRampPalette(c('springgreen4', '#fdfdf8', 'firebrick3'))
  colorLevels <- 3
  cols_to_use = colors(colorLevels)
  
  plotMatmelt$value <- factor(plotMatmelt$value,
                              levels=c(1,0,-1))
  # Plot based on Sign or Mean
  plot <- ggplot(plotMatmelt, aes(x = Var1, y = Var2, fill = value)) +
    labs(x = rowLabel, y = colLabel, fill = "Sign") +
    geom_tile(color = 'gray60') +
    scale_fill_manual(breaks = levels(plotMatmelt$value),
                      values = cols_to_use,
                      labels = c('+', '0', '-')) +
  
    theme(plot.background = element_rect(fill = "white"),
          legend.background = element_rect(fill = "white"),
          panel.border = element_rect(fill = NA, color = NA),
          legend.margin = margin(l = 1, unit = 'cm'),
          legend.title = element_text(hjust = 0.1),
          legend.key.width = unit(1, 'cm'),
          legend.key.height = unit(4, 'cm'),
          legend.text = element_text(size = 12),
          axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0)) +
    scale_x_discrete(expand = c(0, 0))+
    scale_y_discrete(expand = c(0, 0))
  plot
  
  return(plot)
  
  }

  

# PLOT --------------------------------------------------------------------
covariates <- colnames(fitSepTF$XData)

mar=c(1,20,1,1)
mar=c(0,10,0,0)

if(length(covariates) == 1){
  if(covariates == 'tmean_year'){
    alt_name <- c('Yearly temp.')
  }else if(covariates == 'prec_year'){
    alt_name <- c('Yearly prec.')
  }else if(covariates == 'hh'){
    alt_name <- c('Heterogeneity')
  }
}else if(length(covariates) == 3){
  if(identical(covariates, c('tmean_year','prec_year','hh'))){
    alt_name <- c('Yearly temp.','Yearly prec.','Heterogeneity')
  }
}

if(length(colnames(fitSepTF$TrData))==1){
  if(colnames(fitSepTF$TrData) == 'foraging_guild_consensus'){
    alt_char <- 'Guild'
  }else if(colnames(fitSepTF$TrData) == 'Migration_a3_DOF'){
    alt_char <- 'Migratory strategy'
  }
}else if (length(colnames(fitSepTF$TrData))==2){
  alt_char <- c('Migratory strategy / Guild')
}
  

pdf(file=file.path(input,'results','posterior-gamma.pdf'),
    width = 6,
    height = 6)

p<- plotGammaSign(fitSepTF,parName = 'Gamma',
                 plotType = 'Sign',
                 newEnvNames = alt_name,
                 characteristic = alt_char,
                 supportLevel = 0.9)
print(p)

dev.off()
}

