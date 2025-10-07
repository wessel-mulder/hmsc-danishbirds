# edit function
plotVariancePartitioning2 =
  function (hM, VP, cols=NULL, main = "Variance Partitioning", ...)
  {
    ng = dim(VP$vals)[1]
    if(is.null(cols)){
      cols = heat.colors(ng, alpha = 1)
    }
    leg = VP$groupnames
    for (r in seq_len(hM$nr)) {
      leg = c(leg, paste("Random: ", hM$rLNames[r], sep = ""))
    }
    means = round(100 * rowMeans(VP$vals), 1)
    for (i in 1:ng) {
      leg[i] = paste(leg[i], " (mean = ", toString(means[i]),
                     ")", sep = "")
    }
    
    VP$vals <- VP$vals[,rev(colnames(VP$vals))]
    barplot(VP$vals, main = main, xlab= " ", ylab = "",las = 1, horiz=T,
            legend = leg, col = cols,
            #las = 2,
            border = NA,
            space=0.1,
            ann=T,
            xlim=c(0,1),
            args.legend = list(x = 'topright',
                               inset = c(0.02,0.05),
                               bg = 'white'),...)
    #rotate 60 degrees (srt = 60)
  }


pdf(file=file.path(input,'results','VP.pdf'),
    width = 10,
    height = 15)
par(mar = c(5,10,2,2))

VP_2 <- VP

# set plotting parameters 
if(VP$groupnames == 'tmean_year'){
  alt_name <- c('Yearly temp.')
  colorway <- c('firebrick3','cornsilk2')
}else if(VP$groupnames == 'prec_year'){
  alt_name <- c('Yearly prec.')
  colorway <- c('dodgerblue3','cornsilk2')
}else if(VP$groupnames == 'hh'){
  alt_name <- c('Heterogeneity')
  colorway <- c('goldenrod2','cornsilk2')
}

VP_2$groupnames <- alt_name


colnames(VP_2$vals) <- sub("^([A-Za-z])[A-Za-z]+_([a-z]+)$", "\\1. \\2", colnames(VP_2$vals))
plotVariancePartitioning2(fitSepTF,
                          cols = colorway,
                          VP_2,
                          main = 'Proportion of variance explained')

# get MF vals 
MF <- readRDS(file.path(input,'model-outputs','model-fit.rds'))
tjurs <- MF$TjurR2
VP_2$vals <- sweep(VP_2$vals, 2, tjurs, "*")
plotVariancePartitioning2(fitSepTF,
                          cols = colorway,
                          VP_2,
                          main = 'Total variance explained')
dev.off()
