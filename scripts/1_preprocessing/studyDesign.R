# GETTING STARTED ---------------------------------------------------------
library('sf')

# get information on the sites 
grid <- vect('../data/distributions/dof_atlas/DK5km_ED50grid_approx_kvadrkod_DOF/DK5km_ED50grid_approx_kvadrkod_DOF.shp')
plot(grid)

# initiate a new document 
df <- data.frame(matrix(nrow=length(grid),
                        ncol = 0))

# assign site 
df$site <- grid$kvadratkod

# assign site lat and longs 
df$lat <- grid$center_lat
df$lon <- grid$center_lng

df_tripled <- do.call(rbind, lapply(1:3, function(i) {
  df_copy <- df
  df_copy$atlas <- i
  df_copy$survey <- paste0(df_copy$site,'_',i)
  return(df_copy)
}))

write.csv(df_tripled,
          'data/0_data/design/studyDesign.csv',
          row.names = F)

