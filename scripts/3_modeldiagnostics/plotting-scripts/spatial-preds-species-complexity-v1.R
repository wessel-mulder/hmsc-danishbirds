

# ORIGINAL DATA  ----------------------------------------------------------
# get environment and basic info
grid <- fitSepTF$XData
covariates <- colnames(grid)
nspecies <- ncol(preds)
species <- colnames(preds)

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
og_xy_species <- merge(merge,og_Y,by='row.names')
rownames(og_xy_species) <- og_xy_species$Row.names
head(og_xy_species)

# PREDICTION DATA  --------------------------------------------------------
# get richnes 
xy_species <- merge(merge,preds,by='row.names')
rownames(xy_species) <- xy_species$Row.names
head(xy_species[,5:7])
head(og_xy_species[,5:7])

# PLOTTING ----------------------------------------------------------------
# make a continuous color palette
pal <- brewer.pal(9,'YlGn')
pal_raw <- pal
pal <- colorRampPalette(pal)

# and a diffpall
diffpal <- rev(brewer.pal(11,'RdYlBu'))
diffpal <- colorRampPalette(diffpal)


# scale richness to the palette
ncolz <- 100

head(xy_species)

# genus name change 
abbreviate_genus <- function(x) {
  sapply(x, function(nm) {
    parts <- unlist(strsplit(nm, "_"))  # split by underscore
    paste0(substr(parts[1], 1, 1), ". ", parts[2])
  })
}

# define legend range
#zlim <- range(xyrich$rich, na.rm = TRUE)
#zseq <- seq(zlim[1], zlim[2], length.out = ncolz)
if(!dir.exists(file.path(input,'results','sp-preds-singlespecies'))) {dir.create(file.path(input,'results','sp-preds-singlespecies'))}
specie <- species[1]
for(specie in species){
  
pdf(file=file.path(input,'results','sp-preds-singlespecies',paste0(specie,'.pdf')),
    width = 10,
    height = 10)

par(mar=c(10,5,5,5))

nicename <- abbreviate_genus(specie)
# PREDICTIONS 
cols <- pal(ncolz)[as.numeric(cut(
  xy_species[,specie], 
  breaks = seq(0, 1, length.out = ncolz),
  include.lowest = T
  ))]
# plot the points
plot(xy_species$X, xy_species$Y, col = cols, pch = 19,
     xlab = 'X',
     ylab = 'Y',
     main = paste0('Predicted occurrence probability - ',nicename))
# fixed legend from 0 to 12
image.plot(legend.only = TRUE,
           zlim = c(0, 1),      # force scale 0-12
           col = pal(ncolz),
           legend.lab = "Probability",
           horizontal = T)

# OG 
# define a binary palette: 0 = light red, 1 = dark red
binary_cols <- c("indianred1",pal_raw[9])

# assume og_xy_species[,specie] contains 0/1
og_cols <- binary_cols[og_xy_species[,specie] + 1]  # +1 because R indexes start at 1

# plot the points
plot(og_xy_species$X, og_xy_species$Y, col = og_cols, pch = 19,
     xlab = 'X',
     ylab = 'Y',
     main = 'Atlas 3 - Presence/Absence')
# add a binary legend
image.plot(legend.only = TRUE,
           zlim = c(0,1),                   # 0/1 scale
           col = binary_cols,                # binary colors
           legend.shrink = 0.5,
           axis.args = list(at = c(0,1), labels = c("Absent","Present")),
           legend.lab = "",
           horizontal = T)

# DIFFERENCE 
diff <- data.frame(xy_species[,specie] - og_xy_species[,specie])
rownames(diff) <- rownames(xy_species)
xydiff <- merge(merge,diff,by='row.names')
head(xydiff)
colnames(xydiff) <- c('row.names','site','X','Y','diff')
diff_cols <- diffpal(ncolz)[as.numeric(cut(
  xydiff$diff, 
  breaks = seq(-1, 1, length.out = ncolz),
  include.lowest = T,
  ))]

plot(xydiff$X, xydiff$Y, col = diff_cols, pch = 19,
     xlab = 'X',
     ylab = 'Y',
     main = 'Predicted probabilites vs Observed')
image.plot(legend.only = TRUE,
           zlim = c(-1,1),      # force scale 
           col = diffpal(ncolz),
           legend.lab = "Difference",
           horizontal = T)

dev.off()

}
 

