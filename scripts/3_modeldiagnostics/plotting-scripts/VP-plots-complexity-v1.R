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
    
    barplot(VP$vals, main = main, xlab= " ", ylab = "Total variance", las = 1,
            legend = leg, col = cols,...)
    #   mtext("Species", 1,line = 1)
  }


pdf(file=file.path(input,'results','VP.pdf'),
    width = 10,
    height = 5)
par(mar = c(10,5,5,2))

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
tjurs <- MF$TjurR2
VP_2$vals <- VP$vals * tjurs
colnames(VP_2$vals) <- sub("^([A-Za-z])[A-Za-z]+_([a-z]+)$", "\\1. \\2", colnames(VP_2$vals))


plotVariancePartitioning2(fitSepTF,
                          cols = colorway,
                          VP_2,
                          main = 'Total variance explained',
                          las = 2,
                          border = NA,
                          space=0.01,
                          ann=T,
                          ylim=c(0,1),
                          args.legend = list(x = 'topright'))
#inset=c(-0.25,0)))
dev.off()

