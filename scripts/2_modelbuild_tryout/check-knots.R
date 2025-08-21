input <- '.'

# --> LOAD ENVIRONMENT
X <- read.csv(file.path(input,'data/1_preprocessing/X_environmental/X_Environmental.csv'),row.names=1)
X <- X[sort(row.names(X)),]

# keep only sites with data  
X <- na.omit(X)
sites_actual <- row.names(X)


# --> LOAD STUDY DESIGN 
Design <- read.csv(file.path(input,"data/1_preprocessing/design/studyDesign.csv"),row.names=5)

# remove sites without data 
Design <- Design[row.names(Design) %in% sites_actual,]

# sort
Design <- Design[sort(row.names(Design)), ]

# convert to factors
Design$site <- as.factor(Design$site)
Design$atlas <- as.factor(Design$atlas)
Design$year[Design$atlas == '1'] <- 1971
Design$year[Design$atlas == '2'] <- 1992
Design$year[Design$atlas == '3'] <- 2014

xycoords <- as.matrix(xycoords)
plot.new()
par(mfrow=c(2,2))
dist <- c(0.2,0.3,0.4,0.5)

pdf(file=file.path(input,'figs','1_preprocessing','knots.pdf'),
    width = 8,height = 8)
par(mfrow=c(2,2))
for(i in dist){
  #for(j in dist){

    xyKnots = constructKnots(xycoords,minKnotDist = i,knotDist = i)
    n_knots <- nrow(xyKnots)
    plot(xycoords[,2],xycoords[,1],pch=18, asp=1,
         main = paste0('knotDist = ',i,
                       ', nknots = ',n_knots))
    points(xyKnots[,2],xyKnots[,1],col='red',pch=18)
  #}
}
dev.off()


xyKnots = constructKnots(xycoords,knotDist = 0.5, minKnotDist = 0.3)
#pdf(file=file.path(input,'figs','1_preprocessing',paste0(n_knots,'_knots.pdf')))
plot(xycoords[,2],xycoords[,1],pch=18, asp=1)
points(xyKnots[,2],xyKnots[,1],col='red',pch=18)
struc_space <- HmscRandomLevel(sData = xycoords, sMethod = "GPP",
                               sKnot = xyKnots)