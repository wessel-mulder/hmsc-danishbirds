

# ORIGINAL DATA  ----------------------------------------------------------
# get environment and basic info
grid <- fitSepTF$XData
covariates <- colnames(grid)
nspecies <- ncol(preds)
species <- colnames(preds)

# get design and xycoords of sites  
# get design and xycoords of sites  
design <- fitSepTF$studyDesign
if(is.null(design)){
  # get design from a different model 
  temp_fitsep <- readRDS('tmp_rds/mods-complexity-v2/2025-10-14_16-33-15_threeenv_full/fitsepTF.rds')
  design <- temp_fitsep$studyDesign
  
  xycoords <- tryCatch(
    {
      # First attempt
      data.frame(temp_fitsep$ranLevels$site$s@coords[drop = FALSE])
    },
    error = function(e) {
      # Fallback if the above errors
      data.frame(temp_fitsep$ranLevels$site$s)
    }
  )
}else{
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
}


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

# PREDICTION DATA  --------------------------------------------------------
# get richnes 
xy_species <- merge(merge,preds,by='row.names')
rownames(xy_species) <- xy_species$Row.names

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

# get idea of how many atlases are in this set 
replicates <- sub("^.*_", "", rownames(xy_species))
replicates <- unique(replicates)
replicate <- 1
lapply(replicates,function(replicate){
  if(replicate=='1'){year <- '1970s'}
  if(replicate=='2'){year <- '1990s'}
  if(replicate=='3'){year <- '2010s'}
  
  
  print(replicate)
  singlespeciesatlasdir <- file.path(input,'results',paste0('sp-preds-singlespecies-atlas-',replicate))
  if(!dir.exists(singlespeciesatlasdir)){dir.create(singlespeciesatlasdir)}
  
  ### DATA SUBSETS FOR ATLAS 
  xy_species_sub <- xy_species[rownames(xy_species)[grep(paste0("_",replicate,"$"), rownames(xy_species))],,drop=F]
  og_xy_species_sub <- og_xy_species[rownames(og_xy_species)[grep(paste0("_",replicate,"$"), rownames(og_xy_species))],,drop=F]
  

  for(specie in species){
  
    pdf(file=file.path(singlespeciesatlasdir,
                       paste0(specie,'-atlas-',replicate,'.pdf')),
        width = 10,
        height = 10)

    par(mar=c(10,5,5,5))

    nicename <- abbreviate_genus(specie)
    # PREDICTIONS 
    cols <- pal(ncolz)[as.numeric(cut(
      xy_species_sub[,specie], 
      breaks = seq(0, 1, length.out = ncolz),
      include.lowest = T
      ))]
    # plot the points
    plot(xy_species_sub$X, xy_species_sub$Y, col = cols, pch = 19,
         xlab = 'X',
         ylab = 'Y',
         main = paste0('Predicted occurrence probability - ',nicename,' - ',year))
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
    og_cols <- binary_cols[og_xy_species_sub[,specie] + 1]  # +1 because R indexes start at 1
    
    # plot the points
    plot(og_xy_species_sub$X, og_xy_species_sub$Y, col = og_cols, pch = 19,
         xlab = 'X',
         ylab = 'Y',
         main = paste0(year,' - Presence/Absence'))
    # add a binary legend
    image.plot(legend.only = TRUE,
               zlim = c(0,1),                   # 0/1 scale
               col = binary_cols,                # binary colors
               legend.shrink = 0.5,
               axis.args = list(at = c(0,1), labels = c("Absent","Present")),
               legend.lab = "",
               horizontal = T)
    
    # DIFFERENCE 
    diff <- data.frame(diff = xy_species_sub[,specie] - og_xy_species_sub[,specie])
    rownames(diff) <- rownames(xy_species_sub)
    xydiff <- merge(merge,diff,by='row.names')
    head(xydiff)
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
}
)

