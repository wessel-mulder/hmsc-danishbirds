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

### OCCURRENCES 
Y <- read.csv(file.path(start,'data/1_preprocessing/Y_occurrences/Y_occurrences.csv'),row.names=1)

# remove sites without data 
Y <- Y[row.names(Y) %in% sites_actual,]

# grab 13 warblers
genera <- c('Phylloscopus','Curruca','Sylvia','Acrocephalus')
keep <- sapply(strsplit(colnames(Y),'_'),head,1) %in% genera
Y_warblers <- Y[,keep]

# check if they have enough records 
colSums(Y_warblers, na.rm =T)
# barred warbler is absent so we remove it 
Y_warblers <- Y_warblers[colnames(Y_warblers) != "Curruca_nisoria"]
# fixed
colSums(Y_warblers, na.rm =T)


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

# get xycoords and project to UTM 
xycoords <- data.frame(lon = Design3$lon,lat = Design3$lat)
rownames(xycoords) <- Design3$site
library(sf)
sf_coords <- st_as_sf(xycoords, coords=c('lon','lat'),crs=4326)

sf_coords_proj <- st_transform(sf_coords,crs="EPSG:23032")
proj_xycoords <- st_coordinates(sf_coords_proj)
rownames(proj_xycoords) <- rownames(xycoords)
head(proj_xycoords)


# make SpatialPointsDataFrame
xycoords <- data.frame(lon = Design3$lon, lat = Design3$lat, row.names = Design3$site)
sp_coords <- SpatialPointsDataFrame(
  coords = xycoords,
  data   = data.frame(site = rownames(xycoords)),
  proj4string = CRS("+proj=longlat +datum=WGS84") # WGS84
)

plot(sp_coords)
# project to UTM zone 32 (ETRS89 / UTM zone 32N, EPSG:23032)
sp_coords_proj <- spTransform(sp_coords, CRS("+init=epsg:23032"))
plot(sp_coords_proj)
# extract coordinates
proj_xycoords <- coordinates(sp_coords_proj)
rownames(proj_xycoords) <- rownames(xycoords)
colnames(proj_xycoords) <- c('X','Y')

head(proj_xycoords)

library(terra)

xycoords <- data.frame(lon = Design3$lon, lat = Design3$lat)
rownames(xycoords) <- Design3$site

v <- vect(xycoords, geom = c("lon","lat"), crs = "EPSG:4326")
v_proj <- project(v, "EPSG:23032")

proj_xycoords <- crds(v_proj)
rownames(proj_xycoords) <- rownames(xycoords)
head(proj_xycoords)


# SANITY CHECK  -----------------------------------------------------------
### XY DISTRIBUTION OF ENVIRONMENT
tmean <- data.frame(temp=atlas3_temp[,'tmean_year',drop=F])
# strip the "_3" from the rownames
rownames(tmean) <- sub(paste0("_",3,"$"), "", rownames(tmean))
#mergs
tmean_space <- merge(proj_xycoords,tmean,by='row.names')

# plot 
ggplot(tmean_space,
       aes(x=X,y=Y,color=tmean_year))+
  geom_point()+
  scale_color_gradientn(,colors=c('blue','red'))


# TEMPERATURE MODEL  ------------------------------------------------------
Formula_temp <- as.formula(paste("~", paste(colnames(atlas3_temp), collapse = "+"), sep = " "))

knot_dist <- c(25000,50000)
names <- c('25km','50km')
knotmods <- list()
par(mfrow=c(1,1))
for(j in seq_along(knot_dist)){
  print(j)
  i <- knot_dist[j]
  name <- names[j]
  i <- 100000
  name <- '100km'
  
  i <- 25000
  i <- 10000
  
  t <- list(
    test = 'hello'
  )
  
  t[['knotdist']] <- 'bug'
  
  xyKnots = constructKnots(proj_xycoords,knotDist = i,minKnotDist =  i)
  nKnots <- nrow(xyKnots)
  plot(proj_xycoords[,1],proj_xycoords[,2],pch=18, asp=1,
       main=i)
  points(xyKnots[,1],xyKnots[,2],col='red',pch=18)
  print(nKnots)
  
  struc_space <- HmscRandomLevel(sData = proj_xycoords, sMethod = "GPP",
                                 sKnot = xyKnots)
  struc_space$alphapw
  
  dists <- dist(proj_xycoords)
  min(dists)
  max(dists)
  
  freq <- c(0.5,rep(0.005,100))
  samples <- c(0,seq(from = 4999, to = 477312, length.out = 100))
  small <- cbind(samples,freq)
  

  # define m 
  m.temp <-Hmsc(Y = Y_warblers, 
                XData = atlas3_temp, 
                XFormula = Formula_temp,
                studyDesign = Design3[,c('site'),drop=F], 
                ranLevels = list('site'=struc_space),
                distr='probit')
  
  m.temp$ranLevels[[1]]$alphapw
  
  nChains = 4
  test.run = T
  if(test.run){
    thin = 4
    samples = 25
    transient = 100
    verbose =1
  }else{
    thin = 10
    samples = 100
    transient = 1
    verbose = 10
  }
  
  
  m.temp.sampled = sampleMcmc(m.temp, thin = thin, samples = samples, transient = transient,
                                nChains = nChains, nParallel = nChains, verbose = verbose)
  knotmods[[name]] <- m.temp.sampled
}

m1<-knotmods[["25km"]]
m2<-knotmods[["50km"]]
m3 <- knotmods[['100km']]


mpost1<-convertToCodaObject(m1)
summary(mpost1$Alpha[[1]])
plot(mpost1$Alpha[[1]][,1:3])


mpost2<-convertToCodaObject(m2)
summary(mpost2$Alpha[[1]])
plot(mpost2$Alpha[[1]][,1:3])

mpost3<-convertToCodaObject(m3)
summary(mpost3$Alpha[[1]])
plot(mpost3$Alpha[[1]][,1:3])

mpost1$Eta[[1]]

# get param estimates -----------------------------------------------------
mpost <- mpost2
params <- list(mpost$Beta,mpost$Gamma,mpost$V,mpost$Sigma,
               mpost$Eta[[1]],
               mpost$Alpha[[1]],
               mpost$Omega[[1]],
               mpost$Lambda[[1]],
               mpost$Psi[[1]],
               mpost$Delta[[1]])

diags <- list(psrf = list(), ess = list())
chunk_size <- 10

for(j in seq_along(params)){
  print(j)
  mat <- params[[j]]
  ncols <- ncol(mat[[1]]) # get ncol 
  
  # process in chunks if ncols > 1
  if(length(ncols)){
    idx_chunks <- split(1:ncols, ceiling(seq_along(1:ncols)/chunk_size))
    
    psrf_list <- list()
    ess_list  <- list()
    
    # parallelize
    results <- parallel::mclapply(idx_chunks, function(cols) {
      list(
        psrf = gelman.diag(mat[,cols], multivariate=FALSE,autoburnin=F)$psrf[,1],
        ess  = effectiveSize(mat[,cols])
      )
    }, mc.cores = n_cores)  # adjust cores to liking
    
    # Combine
    psrf_list <- lapply(results, `[[`, "psrf")
    ess_list  <- lapply(results, `[[`, "ess")
    
    diags$psrf[[j]] <- unlist(psrf_list)
    diags$ess[[j]]  <- unlist(ess_list)
  }else{
    diags$psrf[[j]] <- gelman.diag(mat, multivariate=FALSE,autoburnin=F)$psrf[,1]
    diags$ess[[j]] <- effectiveSize(mat)
  }
  
}

# plot names 
full_names <- c("Species x Environment",
                "Gamma intercept",
                "Residual covariance",
                'Residual variance in species occurrences',
                'Site loadings',
                'Decay factor of latent, spatial process',
                'Species associations',
                "Species loadings",
                'Local shrinkage - species loadings',
                'Global shrinkage - species loadings')

cols <- c('green4','purple4','orange4','firebrick',
          'lightgreen','ivory','blue4',
          'cornsilk2','cyan1','cyan3')

all_psrf <- list()
all_ess  <- list()

for(i in seq_along(full_names)){
  print(i)
  print(as.numeric(diags$psrf[[i]]))
  all_psrf[[ full_names[i] ]] <- as.numeric(diags$psrf[[i]])
  all_ess[[ full_names[i] ]]  <- as.numeric(diags$ess[[i]])
}

# process a bit to make plotting visible

all_psrf <- lapply(all_psrf, function(x) {
  summary(x)
  x[!is.finite(x)] <- 10
  x[is.na(x)] <- 10
  if(!length(x)){
    x[2] <- x[1] + 0.001
  } else if (length(unique(x)) == 1) {
    # All values are identical
    x[1] <- x[1] + 0.001
    x[2] <- x[2] - 0.001
  }
  x
})

all_ess <- lapply(all_ess, function(x) {
  summary(x)
  x[!is.finite(x)] <- 1
  x[is.na(x)] <- 1
  if(!length(x)){
    x[2] <- x[1] + -1
  } else if (length(unique(x)) == 1) {
    # All values are identical
    x[1] <- x[1] + 1
    x[2] <- x[2] - 2
  }
  x
})

vioplot(all_psrf, col=cols, main="PSRF across parameters",
        ylim=c(0.99, 11),
        xaxt = 'n',ann =  T)
abline(h=1.1, lty=1, col="red") # warning threshold
#axis(1, at=seq_along(all_psrf), labels=names(all_psrf), las=2, cex.axis=0.7)
legend('topleft',
       fill = cols,
       full_names)
print('psrf plot threshold')
vioplot(all_psrf, col=cols, main="PSRF across parameters",
        ylim=c(0.99, 1.1),
        xaxt = 'n',ann =  T)
abline(h=1.1, lty=1, col="red") # warning threshold
#axis(1, at=seq_along(all_psrf), labels=names(all_psrf), las=2, cex.axis=0.7)


