

# ORIGINAL DATA  ----------------------------------------------------------
# get environment and basic info
grid <- fitSepTF$XData
covariates <- colnames(grid)
nspecies <- ncol(preds)

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
og_S = data.frame(rich = rowSums(og_Y))
og_xyrich <- merge(merge,og_S,by='row.names')
rownames(og_xyrich) <- og_xyrich$Row.names

max(og_xyrich$rich) 

# PREDICTION DATA  --------------------------------------------------------
# get richnes 
S = data.frame(rich = rowSums(preds))
xyrich <- merge(merge,S,by='row.names')
rownames(xyrich) <- xyrich$Row.names

max(xyrich$rich) 

# DIFFERENCE  -------------------------------------------------------------
diff <- S - og_S
head(S)
head(og_S)
xydiff <- merge(merge,diff,by='row.names')
head(xydiff)
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

# get idea of how many atlases are in this set 
replicates <- sub("^.*_", "", rownames(xyrich))
replicates <- unique(replicates)
replicate <- 3
lapply(replicates,function(replicate){
  if(replicate=='1'){year <- '1970s'}
  if(replicate=='2'){year <- '1990s'}
  if(replicate=='3'){year <- '2010s'}

  pdf(file=file.path(input,'results',paste0('sp-preds-richness-atlas-',replicate,'.pdf')),
      width = 10,
      height = 10)
  
  par(mar=c(10,5,5,5))
  
  xyrich_sub <- xyrich[rownames(xyrich)[grep(paste0("_",replicate,"$"), rownames(xyrich))],,drop=F]
  og_xyrich_sub <- og_xyrich[rownames(og_xyrich)[grep(paste0("_",replicate,"$"), rownames(og_xyrich))],,drop=F]
  xydiff_sub <- xydiff[rownames(xydiff)[grep(paste0("_",replicate,"$"), rownames(xydiff))],,drop=F]
  print(head(xyrich_sub))
  
  # set number of species to go to 
  if(all_species==1){max_bar <- 120}else{max_bar <- nspecies}
  
  
  # PREDICTIONS 
  cols <- pal(ncolz)[as.numeric(cut(
    xyrich_sub$rich, 
    breaks = seq(0, max_bar, length.out = ncolz),
    include.lowest = T
  ))]
  
  # plot the points
  plot(xyrich_sub$X, xyrich_sub$Y, col = cols, pch = 19,
       xlab = 'X',
       ylab = 'Y',
       main = 'Model predictions - richness')
  
  #text(800000,6250000,sum(xyrich_sub$rich))
  
  # fixed legend from 0 to 12
  image.plot(legend.only = TRUE,
             zlim = c(0, max_bar),      # force scale 0-12
             col = pal(ncolz),
             legend.lab = "Richness",
             horizontal = T)
  
  # OG 
  og_cols <- pal(ncolz)[as.numeric(cut(
    og_xyrich_sub$rich, 
    breaks = seq(0, max_bar, length.out = ncolz),
    include.lowest = T
  ))]
  table(is.na(og_cols))
  
  # plot the points
  plot(og_xyrich_sub$X, og_xyrich_sub$Y, col = og_cols, pch = 19,
       xlab = 'X',
       ylab = 'Y',
       main = paste0('Species richness - ',year))
  
  # fixed legend from 0 to 12
  image.plot(legend.only = TRUE,
             zlim = c(0, max_bar),      # force scale 0-12
             col = pal(ncolz),
             legend.lab = "Richness",
             horizontal = T)
  
  # DIFFERENCE 
  # set scale bar depencing on model
  if(all_species==1){
    n <- 6
    if(all_atlas == 1){
      n <- 10
    }
  }else if(all_species == 0){
    n <- 2
    if(all_atlas == 1){
      n <- 4
    }
  }
  
    # set scale bars and aesthetics
  scale_bar_richness <- c(-n,n)
  plot_list <- list(at = c(-n,-n/2,0,n/2,n),
                    labels = c(paste0('< ',-n),
                               paste0(-n/2),
                               '0',
                               paste0(n/2),
                               paste0('> ',n)
                               )
  )
                    
  # Cap values below/above limits
  vals <- pmax(pmin(xydiff_sub$rich, scale_bar_richness[2]), scale_bar_richness[1])
  
  diff_cols <- pal(ncolz)[as.numeric(cut(
    vals, 
    breaks = seq(scale_bar_richness[1], scale_bar_richness[2], length.out = ncolz),
    include.lowest=T
  ))]
  
  if (any(is.na(diff_cols))) {
    stop("Error: There are NA values in diff_cols. Edit the scale bar")
  }
  
  plot(xydiff_sub$X, xydiff_sub$Y, col = diff_cols, pch = 19,
       xlab = 'X',
       ylab = 'Y',
       main = 'Predictions - observed')
  image.plot(legend.only = TRUE,
             zlim = scale_bar_richness,      # force scale 
             col = pal(ncolz),
             legend.lab = "Difference",
             axis.args = plot_list,
             horizontal = T)
  
  dev.off()
  print('spatial preds finished')
})
