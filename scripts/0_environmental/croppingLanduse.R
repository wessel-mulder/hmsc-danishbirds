rm(list=ls())
# GETTING STARTED ---------------------------------------------------------
library(terra)
library(devtools)
install_github('wessel-mulder/supportScripts')
library(supportScripts)
library(pbapply)

# LOADING FILES -----------------------------------------------------------
input <- list.files(path = '../data/environmental/hilda/hildap_vGLOB-1.0-f_tifs',
                    pattern = '\\.tif$',
                    full.names = T)

# get climate to crop to 
climate <- rast('../data/environmental/chelsaCRUTS/season/prec/breeding/CHELSAcruts_prec_breeding_1901.tif')

# get files to rasters
rasters <- lapply(input,rast)

# crop to extent of climate layer
rasters_c <- pblapply(rasters,crop,
                    y=climate)

# write output 
output_path <- '../data/environmental/hilda/hildap_vGLOB-1.0-f_tifs-cropped'
dir.create(output_path,recursive=T)



metaGenerator(output_path)


names <- c(basename(input))

pblapply(seq_along(rasters_c), function(i) {
  writeRaster(rasters_c[[i]],
              filename = file.path(output,names[i]),
              overwrite = TRUE)
  
}
)



