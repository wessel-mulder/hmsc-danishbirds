# GETTING STARTED ---------------------------------------------------------
library('supportScripts')

pkgs <- c('terra')
loadInstall(pkgs)

# LOADING DATA ------------------------------------------------------------
# define atlas years 
atlas <- list(
  'one' = c(1971:1974),
  'two' = c(1993:1996),
  'three' = c(2014:2017)
)

# define inputs and speces 
info <- list(
  'season' = list(
    name = 'chelsa',
    path = '../data/environmental/chelsaCRUTS/season',
    specs =  c('tmean','prec'),
    second_specs = c('year','breeding','winter')
  ),
  'hetero' = list(
    name = 'hilda',
    path = '../data/environmental/hilda/hildap_vGLOB-1.0-f_heterogeneity',
    specs = c('heterogeneity','proportion'),
    second_specs = c("")
    
  )
)

output_path <- 'data/0_environmental/'

# AGGREGATE INFORMATION ---------------------------------------------------
i <- info$hetero
spec <- 'heterogeneity'
second <- ""

for(i in info){
  data <- i
  path <- data$path
  specs <- data$specs
  second_specs <- data$second_specs
  
  for(spec in specs){
    
    for(second in second_specs){
      # grab files from path
      if(second == ""){
        file_path <- file.path(path,spec)
      }else{
        file_path <- file.path(path,spec,second)
      }
      files <- list.files(file_path,
                          "\\.tif$",
                          full.names = T)
      
      for(atlas_seq in seq_along(atlas)){
        a <- atlas[[atlas_seq]]
        # loop over different atlases 
        years <- as.character(a)
        
        print(years)
        # collaps years into pattern 
        pattern <- paste(years,collapse = '|')
        # grab the files that match the pattern
        subset <- files[grepl(pattern, files)]  
        
        # Load rasters individually to a list
        rasters <- lapply(subset, rast)
        
        # Stack them into one SpatRaster
        r_all <- rast(rasters)
        
        # Number of layers in each input raster
        n_layers <- nlyr(rasters[[1]])
        
        # Group matching layers together (e.g., layer 1 of each, layer 2 of each, ...)
        group_index <- rep(1:n_layers, times = length(rasters))
        
        # landuse already aggregated, but climate not

        if (spec %in% c('tmean','prec')){
          if(spec == 'tmean'){
            r_all <- aggregate(r_all, fact = 6,fun = mean) # mean for temperature
          }else if (spec == 'prec'){
            r_all <- aggregate(r_all, fact = 6,fun = mean) # also mean for precipitation
          }
        }

        
        # Apply mean per group
        agg <- tapp(r_all, index = group_index, fun = mean)
        
        # assign raster names 
        names <- names(rasters[[1]])
        
        if(spec == 'proportion'){
          names(agg) <- paste0('prop_',names) # name original names of landuse types 
        }else if (spec == 'heterogeneity'){
          names(agg) <- spec #name heterogeneity
        }else if(spec == 'tmean'){
          if(second == 'winter'){
            names(agg) <- paste0('tmp_min_',second)
            
          }else{
            names(agg) <- paste0('tmp_avg_',second)
            
          }
        }else if(spec == 'prec'){
          names(agg) <- paste0('prc_sum_',second)
        }
        
        # assign yearly values 
        time(agg,tstep='years') <- rep(a[1],nlyr(agg))
        
        if(second == ""){
          filename<-paste0(data$name,'_',spec,'_atlas',atlas_seq,'.tif')
        }else{
          filename<-paste0(data$name,'_',spec,'_',second,'_atlas',atlas_seq,'.tif')
        }
        writeRaster(agg,
                    file.path(output_path,paste0('atlas',atlas_seq),filename),
                    gdal = c("COMPRESS=DEFLATE"),filetype = 'GTiff',overwrite = TRUE)
        
        }
        
    }
  }
  metaGenerator(output_path)
  }


