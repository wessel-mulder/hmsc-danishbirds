
# GETTING STARTED ---------------------------------------------------------
message("Running in RStudio")
library(sp)
library(sf)
library(terra)
input <- '.'

### filter out ocean plots 
coast <- st_read(file.path(input,'data/1_preprocessing/coastline/EEA_Coastline_Polygon_Shape/EEA_Coastline_20170228.shp'))
grids <- st_read(file.path(input, 'data/1_preprocessing/atlas-grids/DOF_Shapefiles_/DK5km_ED50grid_approx_kvadrkod_DOF.shp'))

coast_proj <- st_transform(coast,st_crs(grids))

intersections <- st_intersection(grids, coast_proj)
intersections$land_area <- st_area(intersections)
grids$grid_area <- st_area(grids)

hist(intersections$land_area)

plot(intersections)

# sum land area within each grid
land_fraction <- intersections %>%
  st_set_geometry(NULL) %>%  # drop geometry for aggregation
  group_by(kvadratkod) %>%      # replace with your grid ID column
  summarise(land_area_sum = sum(as.numeric(land_area)))

#join back with original grids
grids_area <- grids %>%
  left_join(land_fraction, by = "kvadratkod") %>%
  mutate(
    land_area_sum = ifelse(is.na(land_area_sum), 0, land_area_sum), # grids with no land
    pct_land = 100 * land_area_sum / as.numeric(grid_area)
  )

library(ggplot2)

ggplot(grids_area) +
  geom_sf(aes(fill = pct_land)) +
  scale_fill_viridis_c(limits = c(0, 100)) +  # force 0–100%
  labs(title = paste0(length(grids_area$pct_land),
                      ' grids')) +
  theme_minimal()

ggplot(grids_area[grids_area$pct_land>25,]) +
  geom_sf(aes(fill = pct_land)) +
  scale_fill_viridis_c(limits = c(0, 100)) +  # force 0–100%
  labs(title = paste0(length(grids_area$pct_land[grids_area$pct_land>25]),
                      ' grids')) +
  theme_minimal()

st_write(grids_area,file.path(input,'data/1_preprocessing/atlas-grids/grids-ocean-thresholds/grids_ocean_thresholds.shp'))

plot(test)
