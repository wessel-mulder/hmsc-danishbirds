
# PLOTTING FUNCTION -------------------------------------------------------
m <- fitSepTF
parName = 'Beta'
plotType = 'Sign'
plotBetaSign <- function(m,
                             parName = c("Beta", "Gamma", "Omega"),
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
    # function to convert a name
    abbreviate_genus <- function(x) {
      sapply(x, function(nm) {
        parts <- unlist(strsplit(nm, "_"))  # split by underscore
        paste0(substr(parts[1], 1, 1), ". ", parts[2])
      })
    }
    
    colNames <- abbreviate_genus(m$spNames)   # Species
    
    rowLabel <- "Environmental variable"
    colLabel <- "Species"
    
    # Format the data as sign or mean
      toPlot = sign(meanEst)
      toPlot = toPlot * ((supportEst > supportLevel) + (supportEst < (1 - supportLevel)) > 0)

    
    # Remove intercept for Beta (if it's present)
    toPlot <- toPlot[-1, ]  # Remove intercept row
    rowNames <- rowNames[-1]  # Adjust row names
    
    # Remove intercept and set row and column names just before plotting
    if (!is.null(newEnvNames)) {
      rowNames <- newEnvNames  # Use provided new names
    }
  
  
  # Format the data as matrix and add column and row names
  plotMat = matrix(toPlot, nrow = length(rowNames), ncol = length(colNames))
  colnames(plotMat) <- colNames
  rownames(plotMat) <- rowNames
  
  # Reorder matrix based on parameter type
  plotMat <- plotMat[, rev(order(colnames(plotMat))),drop=F]  # Reorder columns
  
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

if(covariates == 'tmean_year'){
  alt_name <- c('Yearly temp.')
}else if(covariates == 'prec_year'){
  alt_name <- c('Yearly prec.')
}else if(covariates == 'hh'){
  alt_name <- c('Heterogeneity')
}

pdf(file=file.path(input,'results','posterior-beta.pdf'),
    width = 6,
    height = 6)

p<- plotBetaSign(fitSepTF,parName = 'Beta',
                 plotType = 'Sign',
                 newEnvNames = alt_name)
print(p)

dev.off()

