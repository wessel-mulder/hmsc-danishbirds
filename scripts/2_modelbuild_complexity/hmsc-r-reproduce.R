rm(list = ls())

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
  input <- '.'
  python <- file.path(getwd(), 'hmsc-hpc-main',"hmsc-venv", "bin", "python3.11")
  flagInit = 0
  flagFitR = 1
} 

# check if accurately installed
summary(TD)

# VIGNETTE 4 --------------------------------------------------------------
set.seed(6)

# records 
n = 100 # number of observations
ns = 5 # number of species 

# parameters 
beta1 = c(-2,-1,0,1,2) #beta parameters diffferent species 
alpha = rep(0,ns) # intercepts for all species 

beta = cbind(alpha,beta1) 

# response
x = cbind(rep(1,n),rnorm(n))
#linear predictors 
Lf = x%*%t(beta)
# xy coords 
xycoords = matrix(runif(2*n),ncol=2)
colnames(xycoords) = c("x-coordinate","y-coordinate")
rownames(xycoords) = 1:n

# spatial variance and decay parameters 
sigma.spatial = c(2)
alpha.spatial = c(0.35)

Sigma = sigma.spatial^2*exp(-as.matrix(dist(xycoords))/alpha.spatial)

# site loadings 
eta1 = mvrnorm(mu=rep(0,n), Sigma=Sigma)

# species loadings
lambda1 = c(1,2,-2,-1,0)

# full linear predictor 
Lr = eta1%*%t(lambda1)
L = Lf + Lr

# generate response 
y = as.matrix(L + matrix(rnorm(n*ns),ncol=ns))

# generate binary response 
yprob = 1*((L +matrix(rnorm(n*ns),ncol=ns))>0)
#store predictors 
XData = data.frame(x1=x[,2])

# visualize based on proximity
rbPal = colorRampPalette(c('cyan','red'))
par(mfrow=c(2,3))
Col = rbPal(10)[as.numeric(cut(x[,2],breaks = 10))]
plot(xycoords[,2],xycoords[,1],pch = 20,col = Col,main=paste('x'), asp=1)
for(s in 1:ns){
  Col = rbPal(10)[as.numeric(cut(y[,s],breaks = 10))]
  plot(xycoords[,2],xycoords[,1],pch = 20,col = Col,main=paste('Species',s), asp=1)
}

# spatially explicit model 
studyDesign = data.frame(sample = as.factor(1:n))
rL.spatial = HmscRandomLevel(sData = xycoords)
rL.spatial = setPriors(rL.spatial,nfMin=1,nfMax=1) #We limit the model to one latent variables for visua
m.spatial = Hmsc(Y=yprob, XData=XData, XFormula=~x1,
                 studyDesign=studyDesign, ranLevels=list("sample"=rL.spatial),distr="probit")

nChains=2
thin = 10
samples = 1000
transient = 1000
verbose = 100

m.spatial = sampleMcmc(m.spatial, thin = thin, samples = samples, transient = transient,
                       nChains = nChains, nParallel = nChains, verbose = verbose,
                       updater=list(GammaEta=FALSE))

#Explanatory power
preds.spatial = computePredictedValues(m.spatial)
MF.spatial = evaluateModelFit(hM=m.spatial, predY=preds.spatial)
MF.spatial

#Predictive power 
#Predictive power
partition = createPartition(m.spatial, nfolds = 2, column = "sample")
cvpreds.spatial = computePredictedValues(m.spatial, partition=partition,
                                         nParallel = nChains, updater=list(GammaEta=FALSE))

cvMF.spatial = evaluateModelFit(hM=m.spatial, predY=cvpreds.spatial)
cvMF.spatial

# looking at alpha parameter: 
mpost.spatial = convertToCodaObject(m.spatial)
plot(mpost.spatial$Alpha[[1]])
summary(mpost.spatial$Alpha[[1]])

# what happens when we fit a non-spatial model 
m = Hmsc(Y=yprob, XData=XData, XFormula=~x1, studyDesign = studyDesign, distr="probit")
m = sampleMcmc(m, thin = thin, samples = samples, transient = transient,
               nChains = nChains, nParallel = nChains, verbose = verbose)

# model fit for non-spatial 
preds = computePredictedValues(m)
MF = evaluateModelFit(hM=m, predY=preds)
MF

partition = createPartition(m, nfolds = 2, column = "sample")
preds = computePredictedValues(m, partition=partition, nParallel = nChains)

MF = evaluateModelFit(hM=m, predY=preds)
MF


# MAKING A SPATIAL MODEL WITH KNOTS ---------------------------------------
par(mfrow=c(1,1))
Knots = constructKnots(xycoords, knotDist = 0.2, minKnotDist = 0.4)

plot(xycoords[,1],xycoords[,2],pch=18, asp=1)
points(Knots[,1],Knots[,2],col='red'
       ,pch=18)

rL.gpp = HmscRandomLevel(sData = xycoords, sMethod = 'GPP'
                         , sKnot = Knots)
rL.gpp = setPriors(rL.gpp,nfMin=1,nfMax=1)
m.gpp = Hmsc(Y=yprob, XData=XData, XFormula=~x1,
             studyDesign=studyDesign, ranLevels=list("sample"=rL.gpp),distr="probit")
m.gpp = sampleMcmc(m.gpp, thin = thin, samples = samples, transient = transient,
                   nChains = nChains, nParallel = nChains, verbose = verbose,
                   updater=list(GammaEta=FALSE))
preds.gpp = computePredictedValues(m.gpp, updater=list(GammaEta=FALSE))
MF.gpp = evaluateModelFit(hM=m.gpp, predY=preds.gpp)
MF.gpp

# cross validation 
cvpreds.gpp = computePredictedValues(m.gpp, partition=partition,
                                     nParallel = nChains, updater=list(GammaEta=FALSE))

cvMF.gpp = evaluateModelFit(hM=m.gpp, predY=cvpreds.gpp)
cvMF.gpp

mpost.gpp = convertToCodaObject(m.gpp)
plot(mpost.gpp$Alpha[[1]])

summary(mpost.gpp$Alpha[[1]])

