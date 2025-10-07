

# ORIGINAL DATA  ----------------------------------------------------------
# get environment and basic info
grid <- fitSepTF$XData
covariates <- colnames(grid)
nspecies <- ncol(preds)

# get design and xycoords of sites  
design <- fitSepTF$studyDesign
xycoords <- tryCatch(
  {
    # First attempt
    data.frame(fitSepTF$ranLevels$site$s@coords[drop = FALSE])
  },
  error = function(e) {
    # Fallback if the above errors
    data.frame(fitSepTF$ranLevels$site$s)
  }
)

head(xycoords)
# rename X and Y consistentely 
col_with_max <- colnames(xycoords)[which.max(lapply(xycoords, max))]
names(xycoords)[names(xycoords) == col_with_max] <- "Y"  # Rename it to "Y"
col_with_min <- colnames(xycoords)[which.min(lapply(xycoords, min))]
names(xycoords)[names(xycoords) == col_with_min] <- "X"  # Rename it to "Y"

# organize columns
xycoords <- xycoords[,c('X','Y')]

# merge by rownames  
merge <- merge(design,xycoords, by.x = 'site', by.y = 'row.names', all.x =T)
rownames(merge) <- rownames(design)


# OG RICHNESS  ------------------------------------------------------------
og_Y = fitSepTF$Y
og_S = data.frame(rich = rowSums(og_Y))
og_xyrich <- merge(merge,og_S,by='row.names')
rownames(og_xyrich) <- og_xyrich$Row.names


# PREDICTION DATA  --------------------------------------------------------
# get richnes 
S = data.frame(rich = rowSums(preds))
xyrich <- merge(merge,S,by='row.names')
rownames(xyrich) <- xyrich$Row.names


# DIFFERENCE  -------------------------------------------------------------
diff <- S - og_S
xydiff <- merge(merge,diff,by='row.names')
rownames(xydiff) <- xydiff$Row.names

table(og_xyrich$rich)

# PLOTTING ----------------------------------------------------------------
# make a continuous color palette
pal <- rev(brewer.pal(11,'RdYlBu'))
pal <- colorRampPalette(pal)

# scale richness to the palette
ncolz <- 100


# define legend range
zlim <- range(xyrich$rich, na.rm = TRUE)
zseq <- seq(zlim[1], zlim[2], length.out = ncolz)

pdf(file=file.path(input,'results','sp-preds-richness.pdf'),
    width = 10,
    height = 10)

par(mar=c(10,5,5,5))


# PREDICTIONS 
cols <- pal(ncolz)[as.numeric(cut(
  xyrich$rich, 
  breaks = seq(0, nspecies, length.out = ncolz),
  include.lowest = T
))]

# plot the points
plot(xyrich$X, xyrich$Y, col = cols, pch = 19,
     xlab = 'X',
     ylab = 'Y',
     main = 'Model predictions - richness')

# fixed legend from 0 to 12
image.plot(legend.only = TRUE,
           zlim = c(0, nspecies),      # force scale 0-12
           col = pal(ncolz),
           legend.lab = "Richness",
           horizontal = T)

# OG 
og_cols <- pal(ncolz)[as.numeric(cut(
  og_xyrich$rich, 
  breaks = seq(0, nspecies, length.out = ncolz),
  include.lowest = T
))]
table(is.na(og_cols))

# plot the points
plot(og_xyrich$X, og_xyrich$Y, col = og_cols, pch = 19,
     xlab = 'X',
     ylab = 'Y',
     main = 'Atlas 3 - richness')

# fixed legend from 0 to 12
image.plot(legend.only = TRUE,
           zlim = c(0, nspecies),      # force scale 0-12
           col = pal(ncolz),
           legend.lab = "Richness",
           horizontal = T)

# DIFFERENCE 
diff_cols <- pal(ncolz)[as.numeric(cut(
  xydiff$rich, 
  breaks = seq(scale_bar_richness[1], scale_bar_richness[2], length.out = ncolz),
  include.lowest=T
))]
if (any(is.na(diff_cols))) {
  stop("Error: There are NA values in diff_cols. Edit the scale bar")
}
plot(xydiff$X, xydiff$Y, col = diff_cols, pch = 19,
     xlab = 'X',
     ylab = 'Y',
     main = 'Predictions - observed')
image.plot(legend.only = TRUE,
           zlim = scale_bar_richness,      # force scale 
           col = pal(ncolz),
           legend.lab = "Difference",
           horizontal = T)


dev.off()
print('spatial preds finished')
