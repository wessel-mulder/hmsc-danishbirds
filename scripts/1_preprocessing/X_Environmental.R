rm(list = ls())
# GETTING STARTED ---------------------------------------------------------
library(dplyr)

# ATLAS MERGING -----------------------------------------------------------------
paths <- paste0('data/0_environmental/atlas',c(1,2,3))

# load all atlas information 
dfs <- lapply(paths,function(atlas_path){
  # grab atlases 
  atlas <- list.files(path = atlas_path,
                   pattern = '\\.csv$',
                   full.names = T)
  # load csvs 
  atlas_csv <- lapply(atlas,read.csv)
  
  # left join sequentially 
  atlas_merged <- Reduce(function(x, y) left_join(x, y, by = "survey"), atlas_csv)
  
  # select appropriate columns 
  atlas_selected <- atlas_merged %>%
    select(survey,
           tmean_year,tmean_winter,tmean_breeding,
           prec_year,prec_winter,prec_breeding,
           hh,unique,
           LULC_0,LULC_11,LULC_22,LULC_33,LULC_44,LULC_55,LULC_66,LULC_77)
}
)

# merge into one 
dfs_merged <- do.call(rbind,dfs)

# write csvs 
write.csv(dfs_merged,
          'data/1_preprocessing/X_environmental/X_Environmental.csv',
          row.names = F)


