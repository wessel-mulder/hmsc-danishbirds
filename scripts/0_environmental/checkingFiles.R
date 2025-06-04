# GETTING STARTED ---------------------------------------------------------
library(terra)
library(sf)
library(dplyr)
library(supportScripts)
library(ggplot2)
library(pbapply)

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



# PLOTTING LANDUSE ATLAS 1 --------------------------------------------------------
# load data - heterogeneiety 
a1_hh <- st_read('data/0_environmental/atlas1/hilda_atlas1_heterogeneiety_shapefile/')
a1_hh

ggplot(a1_hh)+
  geom_sf(aes(fill=hh)) + 
  scale_fill_gradientn(
    colours = c("#1a9850", "#ffffbf", "#d73027"),
    limits = c(0,0.875), # max hetero number
    name = "Heterogeneity") +
  theme_minimal() +
  labs(title = paste0('Heterogeneity - Atlas 1'),
       subtitle = '1971-1974')

# load data - unique 
a1_unique <- st_read('data/0_environmental/atlas1/hilda_atlas1_unique_shapefile/')
a1_unique

ggplot(a1_unique)+
  geom_sf(aes(fill=unique)) + 
  scale_fill_gradientn(
    colours = c("#08306B", "#ffffbf","#d73027"),
    limits = c(1,8),
    name = "Unique") +
  theme_minimal() +
  labs(title = paste0('Unique - Atlas 1'),
       subtitle = '1971-1974')

# load data - proportions 
a1_proportions <- st_read('data/0_environmental/atlas1/hilda_atlas1_proportions_shapefile/hilda_atlas1_proportions.shp')
a1_proportions

# define colors
colors <- c("#468FAF",  # 00 ocean (blue)
            "#D00000",  # 11 urban (red)
            "#FCE762",  # 22 cropland (yellow)
            "#DD7230",  # 33 pasture/rangeland (orange)
            "#679436",  # 44 forest (green)
            "#7DD181",  # 55 unmanaged grass/shrubland (light green)
            "#DCEED1",  # 66 sparse/no vegetation (gray)
            "#AFD0D6")  # 77 water (cyan)

colors_dark <- c(
  "#2C6D84",  # Ocean (was #468FAF)
  "#9B0000",  # Urban (was #D00000)
  "#D6C200",  # Cropland (was #FCE762)
  "#A54F00",  # Pasture/rangeland (was #DD7230)
  "#4B6F29",  # Forest (was #679436)
  "#4D9B5F",  # Grass/shrubland (was #7DD181)
  "#A8BBA4",  # Other / Sparse veg (was #DCEED1)
  "#6A9DA6"   # Water (was #AFD0D6)
)

landuse_labels <- c("Ocean", "Urban", "Cropland", "Pasture/rangeland",
                    "Forest", "Grass/shrubland", "Other", "Water")
landuse_vals <- paste0('LULC_',c(0,11,22,33,44,55,66,77))
plots_a1 <- list()
for(i in seq_along(landuse_labels)){
  plots_a1[[i]] <- 
    ggplot(a1_proportions)+
    geom_sf(aes(fill=!!sym(landuse_vals[i]))) + 
    scale_fill_gradientn(
      colours = c("#F7F7F7", colors_dark[i]),
      limits = c(0,1),
      name = landuse_labels[i]) +
    theme_minimal() +
    labs(title = paste0(landuse_labels[i],' - Atlas 1'),
         subtitle = '1971-1974')
}

for (i in seq_along(plots_a1)) {
  print(plots_a1[[i]])
}

# PLOTTING LANDUSE ATLAS 3 --------------------------------------------------------
# load data - heterogeneiety 
a3_hh <- st_read('data/0_environmental/atlas3/hilda_atlas3_heterogeneiety_shapefile/')
a3_hh

ggplot(a3_hh)+
  geom_sf(aes(fill=hh)) + 
  scale_fill_gradientn(
    colours = c("#1a9850", "#ffffbf", "#d73027"),
    limits = c(0,0.875), # max hetero number
    name = "Heterogeneity") +
  theme_minimal() +
  labs(title = paste0('Heterogeneity - Atlas 3'),
       subtitle = '2014-2017')

# load data - unique 
a3_unique <- st_read('data/0_environmental/atlas3/hilda_atlas3_unique_shapefile/')
a3_unique

ggplot(a3_unique)+
  geom_sf(aes(fill=unique)) + 
  scale_fill_gradientn(
    colours = c("#08306B", "#ffffbf","#d73027"),
    limits = c(1,8),
    name = "Unique") +
  theme_minimal() +
  labs(title = paste0('Unique - Atlas 3'),
       subtitle = '2014-2017')

# load data - proportions 
a3_proportions <- st_read('data/0_environmental/atlas3/hilda_atlas3_proportions_shapefile/hilda_atlas3_proportions.shp')
a3_proportions

# define colors
colors <- c("#468FAF",  # 00 ocean (blue)
            "#D00000",  # 11 urban (red)
            "#FCE762",  # 22 cropland (yellow)
            "#DD7230",  # 33 pasture/rangeland (orange)
            "#679436",  # 44 forest (green)
            "#7DD181",  # 55 unmanaged grass/shrubland (light green)
            "#DCEED1",  # 66 sparse/no vegetation (gray)
            "#AFD0D6")  # 77 water (cyan)

colors_dark <- c(
  "#2C6D84",  # Ocean (was #468FAF)
  "#9B0000",  # Urban (was #D00000)
  "#D6C200",  # Cropland (was #FCE762)
  "#A54F00",  # Pasture/rangeland (was #DD7230)
  "#4B6F29",  # Forest (was #679436)
  "#4D9B5F",  # Grass/shrubland (was #7DD181)
  "#A8BBA4",  # Other / Sparse veg (was #DCEED1)
  "#6A9DA6"   # Water (was #AFD0D6)
)

landuse_labels <- c("Ocean", "Urban", "Cropland", "Pasture/rangeland",
                    "Forest", "Grass/shrubland", "Other", "Water")
landuse_vals <- paste0('LULC_',c(0,11,22,33,44,55,66,77))
plots_a3 <- list()
for(i in seq_along(landuse_labels)){
  plots_a3[[i]] <- 
    ggplot(a3_proportions)+
    geom_sf(aes(fill=!!sym(landuse_vals[i]))) + 
    scale_fill_gradientn(
      colours = c("#F7F7F7", colors_dark[i]),
      limits = c(0,1),
      name = landuse_labels[i]) +
    theme_minimal() +
    labs(title = paste0(landuse_labels[i],' - Atlas 3'),
         subtitle = '2014-2017')
}

for (i in 1:8) {
  paired_plot <- plots_a1[[i]] | plots_a3[[i]]
  print(paired_plot)
  diff <- plots_a3[[i]] - plots_a1[[i]]
  print(diff)
}





