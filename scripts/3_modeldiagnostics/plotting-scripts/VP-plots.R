pdf(file=file.path(input,'results','VP.pdf'),
    width = 14,
    height = 6)

par(mar = c(5,5,5,15))
# ALL VARS 
main <- 'Variance Partitioning'
Hmsc::plotVariancePartitioning(m,VP,
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
                           axisnames=F,
                           ann=T,
                         args.legend = list(x = 'topright',
                                            inset=c(-0.25,0)))

main <- 'Variance Partitioning - Grouped by category'
Hmsc::plotVariancePartitioning(m,VP_split,
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
                         axisnames=F,
                         ann=T,
                         args.legend = list(x = 'topright',
                                            inset=c(-0.25,0)))

main <- 'Variance Partitioning - Grouped by season'
Hmsc::plotVariancePartitioning(m,VP_season,
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
                         axisnames=F,
                         ann=T,
                         args.legend = list(x = 'topright',
                                            inset=c(-0.25,0)))

dev.off()
