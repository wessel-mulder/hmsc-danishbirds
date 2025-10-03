
# COMPUTE ASSOCIATIONS  ---------------------------------------------------


OmegaCor <- computeAssociations(fitSepTF)
supportLevel <- 0.95

# effect of site
toPlot <- ((OmegaCor[[1]]$support > supportLevel) +
             (OmegaCor[[1]]$support < (1 - supportLevel)) > 0) *
  OmegaCor[[1]]$mean

# function to convert a name
abbreviate_genus <- function(x) {
  sapply(x, function(nm) {
    parts <- unlist(strsplit(nm, "_"))  # split by underscore
    paste0(substr(parts[1], 1, 1), ". ", parts[2])
  })
}

# apply to row and column names
rownames(toPlot) <- abbreviate_genus(rownames(toPlot))
colnames(toPlot) <- abbreviate_genus(colnames(toPlot))


# INIT PLOTTING -----------------------------------------------------------



mar=c(1,1,1,1)


pdf(file=file.path(input,'results','posterior-omega-coorplot.pdf'),
    width = 7,
    height = 7)


# reorder using eigenvectors
corrplot(toPlot,
         method = "color",
         type='lower',
         col = colorRampPalette(c("blue", "white", "red"))(200),
         tl.cex = 1,
         tl.col = "black",
         order = "hclust",      # <- this reorders nicely
         title = 'Residual associations',
         mar=mar)
dev.off()

