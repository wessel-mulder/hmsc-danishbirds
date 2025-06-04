# GETTING STARTED ---------------------------------------------------------
library(terra)
library(sf)
library(dplyr)
library(supportScripts)
library(ggplot2)
library(pbapply)
rm(list=ls())
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
    second_specs = c('year','breeding','winter'),
    project_method = 'bilinear'
  ),
  'hetero' = list(
    name = 'hilda',
    path = '../data/environmental/hilda/hildap_vGLOB-1.0-f_tifs-cropped/',
    specs = c("heterogeneiety","number","proportions"),
    second_specs = c(""),
    project_method = 'near'
  )
)

grid <- st_read('../data/distributions/dof_atlas/DK5km_ED50grid_approx_kvadrkod_DOF/DK5km_ED50grid_approx_kvadrkod_DOF.shp')

output_path <- 'data/0_environmental'
dir.create(output_path,recursive=T)

classes_labels <- c('Ocean_00','Urban_11','Cropland_22','Pasture_33',
                    'Forest_44','GrassShrub_55','Other_66','Water_77')

# FUNCTIONS  --------------------------------------------------------------
get_years <- function(years = years,files = files){
  # convert to character 
  years <- as.character(years)

  # collaps years into pattern 
  pattern <- paste(years,collapse = '|')
  # grab the files that match the pattern
  subset <- files[grepl(pattern, files)]  
  
  # Load rasters individually to a list
  rasters <- lapply(subset, rast)
  return(rasters)
}

grid_means_matrix <- function(rast, grid_data) {
  # Extract raster values
  ex <- terra::extract(rast, grid_data$cell)[, 1]
  
  # Keep only rows where ex is not NA
  valid_idx <- !is.na(ex)
  ex <- ex[valid_idx]
  grid_data <- grid_data[valid_idx, ]
  
  # Calculate the weighted sum and the sum of fractions for each polygon ID
  weighted_sum <- tapply(grid_data[,2] * grid_data$fraction, grid_data$ID, sum, na.rm = TRUE)
  fraction_sum <- tapply(grid_data$fraction, grid_data$ID, sum, na.rm = TRUE)
  
  # Calculate weighted means
  weighted_means <- weighted_sum / fraction_sum

    # Convert to matrix (choose number of rows and columns as needed)
  grid_means_df <- data.frame(ID = names(weighted_means),
                              grid_means = weighted_means)  # Change ncol to adjust matrix dimensions
  
  return(grid_means_df)
}

grid_heterogeneiety_matrix <- function(rast, grid_data) {
  # Extract raster values
  ex <- terra::extract(rast, grid_data$cell)[, 1]
  
  # Keep only rows where ex is not NA
  valid_idx <- !is.na(ex)
  ex <- ex[valid_idx]
  grid_data <- grid_data[valid_idx, ]
  
  head(grid_data)
  
  # get individual classes 
  classes <- grid_data %>%
    group_by(ID,LULC_states) %>%
    summarise(fraction = sum(fraction), .groups='drop')
  
  sum <- grid_data %>%
    group_by(ID) %>%
    summarise(sum = sum(fraction), .groups='drop')
  
  join <- left_join(classes,sum,by='ID')
  join$proportion <- join$fraction / join$sum
  
  # get heterogeneity index 
  heterogeneiety <- join %>%
    group_by(ID) %>%
    summarise(hh = 1-sum(proportion^2), .groups='drop')
  
  # get the unique classes 
  unique <- grid_data %>%
    group_by(ID) %>%
    summarise(unique = length(unique(LULC_states)), .groups='drop')
  
  # get the proportional classes 
  join_wide <- join %>%
    select(ID, LULC_states, proportion) %>%  # Keep only the necessary columns
    tidyr::pivot_wider(
    names_from = LULC_states,
    values_from = proportion,
    values_fill = 0,
    names_prefix = "LULC_"
  )
  
  joined <- heterogeneiety %>%
    left_join(unique, by = 'ID') %>%
    left_join(join_wide,by='ID')
  
  joined <- as.data.frame(joined)

  return(joined)
}


                      






# CLIMATE -----------------------------------------------------------------

data <- info$season
path <- data$path
specs <- data$specs
second_specs <- data$second_specs

for(spec in specs){
  
  for(second in second_specs){
    # grab files from path
    file_path <- file.path(path,spec,second)
    
    files <- list.files(file_path,
                        "\\.tif$",
                        full.names = T)
    
    for(atlas_seq in seq_along(atlas)){
      atlas_path <- file.path(output_path,paste0('atlas',atlas_seq))
      dir.create(atlas_path)    
      
      ### DEFINE NAMES SO WE CAN CHECK IF THEY EXIST
      # define name of the csv file
      name <- paste(c(data$name,spec,second,paste0('atlas',atlas_seq)),collapse = '_')
      name_csv <- paste0(name,'.csv')
      # define the name of the shapefiles 
      name_shp <- paste0(name,'.shp')
      name_shp_folder <- paste0(name,'_shapefile')
      
      if(file.exists(file.path(atlas_path,name_shp_folder,name_shp))) {
        next
      }
      

      
      #if(file.exists(paste0('figs/atlas',atlas_seq,'_',spec,'_',second,'.pdf'))){
      #  print(paste('skipping',paste0('figs/atlas',atlas_seq,'_',spec,'_',second,'.pdf')))
      #  next
      #}
      
      # get the right years 
      years <- atlas[[atlas_seq]]
      print(years)
      
      # get the appropriate atlas_sequence
      rasters <- get_years(years = years,files = files)
      
      # project rasters to the same projection as the grid 
      rasters_projected <- pblapply(rasters,project,
                                    y=grid, # specificy the DOF grid 
                                    method = data$project_method) # choose appropriate method based on data (bilinear for temp/prec, nearest neightbour for landuse) 
      
      print(paste0('method used to project: ',data$project_method))
      
      # extract values 
      cells <- pblapply(rasters_projected,extract,
                        y=grid,
                        cells=T,
                        exact=T)
      

      weighted_means <- pblapply(1:length(rasters_projected), function(myrast) {
          # Extract the appropriate raster and cells for this iteration
          rast <- rasters_projected[[myrast]]
          grid_data <- cells[[myrast]]
          
          # Call the grid_means_matrix function with the current raster and cell data
          grid_means_matrix(rast, grid_data)
        })
      


      # get the right IDs
      IDs <- data.frame(site = grid$kvadratkod,
                        ID = as.character(seq_along(grid$kvadratkod)))
      
      # join the values with the appropriate grid codes 
      weighted_means_merged <- pblapply(weighted_means,left_join,
                                        y=IDs,
                                        by='ID')
      
      # all rasters into one dataframe
      stacked <- do.call(rbind,weighted_means_merged)
      
      # get the average per site over the years 
      atlas_means <- pbtapply(stacked$grid_means,stacked$site,mean)
      
      # make a new dataframe with the columns I'm interested in 
      df <- data.frame(kvadratkod = names(atlas_means),
                       survey = paste0(names(atlas_means),'_',atlas_seq),
                       val = atlas_means)
      
      # rename the columns 
      names(df) <- c('kvadratkod','survey',paste(spec,second,sep='_'))

      # merge original grid with data for shapefile 
      grid_merged <- left_join(grid,df,by='kvadratkod')
      
      # write the csv file 
      write.csv(df,
                file.path(atlas_path,name_csv))
      
      # write the shapefile 
      dir.create(file.path(atlas_path,name_shp_folder))
      st_write(grid_merged,
               dsn = file.path(atlas_path,name_shp_folder,name_shp),
               append=F)        }
      
    }
  metaGenerator(output_path)
}



# LANDUSE -----------------------------------------------------------------
rm(list=ls())
data <- info$hetero
path <- data$path
spec <- data$specs
second <- data$second_specs

file_path <- file.path(path)

files <- list.files(file_path,
                    "\\.tif$",
                    full.names = T)

atlas_seq <- 1
spec <- 'proportions'
for(atlas_seq in seq_along(atlas)){

  #if(file.exists(paste0('figs/atlas',atlas_seq,'_',spec,'_',second,'.pdf'))){
  #  print(paste('skipping',paste0('figs/atlas',atlas_seq,'_',spec,'_',second,'.pdf')))
  #  next
  #}
  
  # get the right years 
  years <- atlas[[atlas_seq]]
  print(years)
  
  # get the appropriate atlas_sequence
  rasters <- get_years(years = years,files = files)
  
  # project rasters to the same projection as the grid 
  rasters_projected <- pblapply(rasters,project,
                                y=grid, # specificy the DOF grid 
                                method = data$project_method) # choose appropriate method based on data (bilinear for temp/prec, nearest neightbour for landuse) 
  
  print(paste0('method used to project: ',data$project_method))
  
  # if land-use map turn all different types of forests into one class 

  rasters_projected <- pblapply(rasters_projected, function(r) {
      r[r >= 40 & r <= 45] <- 44
      return(r)
    })
    # and rename the layers to get consistent column names later on
  rasters_projected <- pblapply(rasters_projected, function(x) {
      names(x) <- 'LULC_states'
      x
    })
  
  # extract values 
  cells <- pblapply(rasters_projected,extract,
                    y=grid,
                    cells=T,
                    exact=T)
  
    weighted_means <- pblapply(1:length(rasters_projected), function(myrast) {
      # Extract the appropriate raster and cells for this iteration
      rast <- rasters_projected[[myrast]]
      grid_data <- cells[[myrast]]
      
      # Call the grid_means_matrix function with the current raster and cell data
      grid_heterogeneiety_matrix(rast, grid_data)
    })
    
  # get the right IDs
  IDs <- data.frame(site = grid$kvadratkod,
                    ID = as.character(seq_along(grid$kvadratkod)))
  
  # join the values with the appropriate grid codes 
  weighted_means_merged <- pblapply(weighted_means, function(df) {
    df$ID <- as.character(df$ID)
    left_join(df, IDs, by = "ID")
  })
  
  # all rasters into one dataframe
  stacked <- do.call(rbind,weighted_means_merged)
  
  # get the average per site over the years 
  hh_means <- pbtapply(stacked$hh,stacked$site,mean)
  unique_modes <- pbtapply(stacked$unique,stacked$site,function(x){
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  })
  cols <- paste0('LULC_',c(0,11,22,33,44,55,66,77))
  
  props_means <- as.data.frame(lapply(stacked[cols], function(col) {
      tapply(col, stacked$site, mean)
    }))
  

  # make a new dataframe with the columns I'm interested in 
  df <- data.frame(kvadratkod = names(hh_means),
                   survey = paste0(names(hh_means),'_',atlas_seq),
                   hh = hh_means,
                   unique = unique_modes)
  
  df_joined <- cbind(df,props_means)
  
  ### DEFINE NAMES SO WE CAN CHECK IF THEY EXIST
  # define name of the csv file
  name <- paste(c(data$name,paste0('atlas',atlas_seq)),collapse = '_')
  name_csv <- paste0(name,'.csv')
  
  # write the csv file 
  write.csv(df_joined,
            file.path(atlas_path,name_csv))
  
  names <- c('heterogeneiety','unique','proportions')
  col_choose <- list('hh','unique',cols)
  
  for(i in seq_along(col_choose)){
  col <- col_choose[[i]]
  name_spec <- names[[i]]
  
  # merge original grid with data for shapefile 
  grid_merged <- left_join(grid,df_joined[c(col,'kvadratkod')],by='kvadratkod')
  # define the name of the shapefiles 
  name_shp <- paste0(name,'_',name_spec,'.shp')
  name_shp_folder <- paste0(name,'_',name_spec,'_shapefile')
  # write the shapefile 
  dir.create(file.path(atlas_path,name_shp_folder))
  st_write(grid_merged,
           dsn = file.path(atlas_path,name_shp_folder,name_shp),
           append=F)
  
  }
  
  metaGenerator(output_path)
  
}



# PLOTTING TEMPERATURES ---------------------------------------------------
a1 <- st_read('data/0_environmental/atlas1/chelsa_year_tmean_atlas1_shapefile/chelsa_year_tmean_atlas1.shp') 
a2 <- st_read('data/0_environmental/atlas2/chelsa_year_tmean_atlas2_shapefile/chelsa_year_tmean_atlas2.shp') 
a3 <- st_read('data/0_environmental/atlas3/chelsa_year_tmean_atlas3_shapefile/chelsa_year_tmean_atlas3.shp') 

range(a1$prc_wnt[!is.na(a1$prc_wnt)])

range(a2$val[!is.na(a2$val)])
range(a3$val[!is.na(a3$val)])

plotting_atlases <- list(
  atlas1 = list(
    atlas = a1,
    name = 'Atlas 1',
    years = '(1971-1974)'
  ),
  atlas2 = list(
    atlas = a2,
    name = 'Atlas 2',
    years = '(1993-1996)'
  ),
  atlas3 = list(
    atlas = a3,
    name = 'Atlas 3',
    years = '(2014-2017)'
  )
)

plotting_info <- plotting_atlases[[1]]
lapply(plotting_atlases,function(plotting_info){
  
  data <- plotting_info$atlas
  
  ggplot(data)+
    geom_sf(aes(fill=val)) + 
    scale_fill_gradientn(
      colours = c("#313695", "#74add1", 
                  "#ffffbf", 
                  "#f46d43", "#a50026"),
      #values = c(7, 8.5, 9.5, 10.5, 12),  # Centered near the midpoint (~9.5)
      limits = c(7, 11),
      name = "°C"
    ) +
    theme_minimal() +
    labs(title = paste0('Average Annual Temperature - ',plotting_info$name),
         subtitle = plotting_info$years)
  
  
  #ggsave(paste0('figs/atlas',atlas_seq,'_',spec,'_',second,'.pdf'),
     #    plot=last_plot())
})


  


# PLOTTING PRECIPITATION --------------------------------------------------
a3 <- st_read('data/0_environmental/atlas1/chelsa_prec_winter_atlas3_shapefile/chelsa_prec_winter_atlas3.shp')

plotting_atlases <- list(
  atlas1 = list(
    atlas = a1,
    name = 'Atlas 1',
    years = '(1971-1974)'
  ),
  atlas2 = list(
    atlas = a2,
    name = 'Atlas 2',
    years = '(1993-1996)'
  ),
  atlas3 = list(
    atlas = a3,
    name = 'Atlas 3',
    years = '(2014-2017)'
  )
)

range(a3$prc_wnt[!is.na(a3$prc_wnt)])

plotting_info <- plotting_atlases[[3]]

lapply(plotting_atlases,function(plotting_info){
  
  data <- plotting_info$atlas
  
  ggplot(data)+
    geom_sf(aes(fill=prc_wnt)) + 
    scale_fill_gradientn(
      colours = c("#f7fbff", "#deebf7", "#9ecae1", "#3182bd", "#08519c"),
      #values = c(7, 8.5, 9.5, 10.5, 12),  # Centered near the midpoint (~9.5)
      limits = c(90, 190),
      name = "Precipitation (mm)"
    ) +
    theme_minimal() +
    labs(title = paste0('Total Winter Precipitation - ',plotting_info$name),
         subtitle = plotting_info$years)
  
  
  #ggsave(paste0('figs/atlas',atlas_seq,'_',spec,'_',second,'.pdf'),
  #    plot=last_plot())
})


