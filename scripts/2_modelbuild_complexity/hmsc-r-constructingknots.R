rm(list = ls())
n_cores=4

# GETTING STARTED ---------------------------------------------------------
if (interactive() && Sys.getenv("RSTUDIO") == "1") {
  message("Running in RStudio")
  library(Hmsc)
  library(jsonify)
  library(knitr)
  library(corrplot)
  library(ape)
  library(dplyr)
  library(MASS)
  library(vioplot)
  library(ggplot2)
  
  start <- '.'
  python <- file.path(getwd(), 'hmsc-hpc-main',"hmsc-venv", "bin", "python3.11")
  
} 

# make dirs 
mod <- 'v1_tmean'
input <- file.path('.','tmp_rds','mods-complexity','v1_tmean')
if(!dir.exists(input)) {dir.create(input)}
if(!dir.exists(file.path(input,'model-outputs'))){dir.create(file.path(input,'model-outputs'))}
if(!dir.exists(file.path(input,'results'))){dir.create(file.path(input,'results'))}

# check if accurately installed
summary(TD)


# LOADING BASIC DATA - TEMPERATURE MODEL   ------------------------------------------------------------
### ENVIRONMENT
X <- read.csv(file.path(start,'data/1_preprocessing/X_environmental/X_Environmental.csv'),row.names=1)
X <- X[sort(row.names(X)),]

# keep only sites with data  
X <- na.omit(X)

# keep only atlas 3
atlas3 <- X[rownames(X)[grep("_3$", rownames(X))],,drop=F]

# grab tmean_year
atlas3_temp <- atlas3[,'tmean_year',drop=F]

sites_actual <- row.names(atlas3_temp)

# --> LOAD STUDY DESIGN 
Design <- read.csv(file.path(start,"data/1_preprocessing/design/studyDesign.csv"),row.names=5)

# remove sites without data 
Design <- Design[row.names(Design) %in% sites_actual,]

# sort
Design <- Design[sort(row.names(Design)), ]
Design3 <- Design[rownames(Design)[grep("_3$", rownames(Design))],,drop=F]
Design3 <- Design3[,c('site','lat','lon')]

head(Design3)

# convert to factors
Design3$site <- as.factor(Design3$site)

library(terra)

xycoords <- data.frame(lon = Design3$lon, lat = Design3$lat)
rownames(xycoords) <- Design3$site

v <- vect(xycoords, geom = c("lon","lat"), crs = "EPSG:4326")
v_proj <- project(v, "EPSG:23032")

proj_xycoords <- crds(v_proj)
rownames(proj_xycoords) <- rownames(xycoords)
write.csv(proj_xycoords,
          '~/Desktop/csv-coords.csv')
head(proj_xycoords)
plot(proj_xycoords)


library(igraph)

# Compute distance matrix
dmat <- as.matrix(distance(v_proj))

threshold <- 5050
print(t)
# Build adjacency matrix: 1 if distance <= threshold, else 0
adj <- dmat <= threshold
diag(adj) <- 0  # remove self-links

g <- graph_from_adjacency_matrix(adj, mode = "undirected")
clusters <- components(g)

# Cluster ID for each site
cluster_ids <- clusters$membership

# Add to original data
v_proj$cluster <- cluster_ids
print(table(v_proj$cluster))

# Define a color palette
n_clusters <- length(unique(v_proj$cluster))
colors <- c(
  "#E41A1C", # red
  "#377EB8", # blue
  "#4DAF4A", # green
  "#984EA3", # purple
  "#FF7F00", # orange
  "#FFFF33", # yellow
  "#A65628", # brown
  "#F781BF", # pink
  "#999999", # gray
  "#66C2A5", # teal
  "#FC8D62"  # coral
)

# Plot points, colored by cluster
plot(v_proj, col=colors[v_proj$cluster], pch=16, main="Clusters of Sites")

