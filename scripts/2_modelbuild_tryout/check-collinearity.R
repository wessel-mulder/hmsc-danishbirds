# GETTING STARTED  --------------------------------------------------------
library(dplyr)

# 1. CHECK COLLINEARITY 
input <- '.'
# --> LOAD ENVIRONMENT
X <- read.csv(file.path(input,'data/1_preprocessing/X_environmental/X_Environmental.csv'),row.names=1)
X <- X[sort(row.names(X)),]

# keep only sites with data  
X <- na.omit(X)

# --> LOAD SITES
sites <- read.csv(file.path(input,'data/1_preprocessing/design/studyDesign.csv'))
rownames(sites) <- sites$survey
sites <- select(sites,
                c(lat,lon))
head(sites)

# GET COR MATS 
cor_mat <- cor(X, use = "pairwise.complete.obs")
p.mat <- cor_pmat(cor_mat)

# plot 
if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggcorrplot")
library(ggcorrplot)
ggcorrplot(cor_mat,lab=T,
           type='lower',
           p.mat = p.mat
           ) 

# 2. VIF check for all predictors
#    Trick: fit a "fake" model with all variables
#    (response doesn't matter for vif, so just use the first one)
hist(X$LULC_77)
#X$LULC_77 <- NULL

library(car)
formula_all <- as.formula(
  paste(names(X)[1], "~ .")
)

lm_fit <- lm(formula_all, data = X)
vif_vals <- vif(lm_fit)
vif_vals

# 3. Combine into a quick summary table
summary_tab <- data.frame(
  Variable = names(vif_vals),
  VIF = round(vif_vals, 2)
)
print(summary_tab)

library(performance)
check_collinearity(fitSepTF)

library(usdm)
v <- vifstep(X,th=10)
ex <- exclude(X,v)
colnames(ex)
v <- vifstep(X,th=5)
ex <- exclude(X,v)
colnames(ex)

v <- vifstep(X,th=2.5)
ex <- exclude(X,v)
colnames(ex)

cor_mat <- cor(ex, use = "pairwise.complete.obs")
p.mat <- cor_pmat(cor_mat)
ggcorrplot(cor_mat,lab=T,
           type='full',
           insig = 'blank',
           p.mat = p.mat
) 

cor_mat <- cor(X, use = "pairwise.complete.obs")
p.mat <- cor_pmat(cor_mat)
ggcorrplot(cor_mat,lab=T,
           type='full',
           insig = 'blank',
           p.mat = p.mat
) 

X_simple <- X[,c(1,4,7)]
head(X_simple)
cor_mat <- cor(X_simple, use = "pairwise.complete.obs")
p.mat <- cor_pmat(cor_mat)
ggcorrplot(cor_mat,lab=T,
           type='full',
           insig = 'blank',
           p.mat = p.mat
) 

X_simple <- X[,c(2,5,7)]
head(X_simple)
cor_mat <- cor(X_simple, use = "pairwise.complete.obs")
p.mat <- cor_pmat(cor_mat)
ggcorrplot(cor_mat,lab=T,
           type='full',
           insig = 'blank',
           p.mat = p.mat
) 

X_simple <- X[,c(3,6,7)]
head(X_simple)
cor_mat <- cor(X_simple, use = "pairwise.complete.obs")
p.mat <- cor_pmat(cor_mat)
ggcorrplot(cor_mat,lab=T,
           type='full',
           insig = 'blank',
           p.mat = p.mat
) 


# SPATIAL CORRELATION -----------------------------------------------------
join <- merge(X,sites,by=0) # by rownames
rownames(join) <- join$Row.names
join <- join[,-1]

library('SpatialPack')
# get coords
coords <- cbind(join$lon,join$lat)
# get vars
vars <- join[, !(names(join) %in% c("lon", "lat"))]
# construct mats 
mat <- matrix(nrow=ncol(vars),
              ncol = ncol(vars))
p.mat <- matrix(nrow=ncol(vars),
              ncol = ncol(vars))

# fill in matrix with t_test results 
for(i in seq_along(colnames(vars))){
  for(j in i:length(colnames(vars))){  # j starts at i
    t <- modified.ttest(vars[,i], vars[,j], coords)
    mat[i,j] <- t$corr
    mat[j,i] <- t$corr
    p.mat[i,j] <- t$p.value
    p.mat[j,i] <- t$p.value
  }
  print(i)
}

colnames(mat) <- colnames(vars)
colnames(p.mat) <- colnames(vars)

library(ggcorrplot)
ggcorrplot(mat,lab=T,
           type='full',
           insig='blank',
           p.mat = p.mat) 


# MORAN  ------------------------------------------------------------------
library('spdep')
library('caret')
coords <- cbind(join$lon,join$lat)

# Identify variables to detrend (exclude lon and lat)
#nearZeroVar(join)
sub <- join[,!names(join) %in% c('LULC_66','lat','lon')]
sub <- join[,1:8]
vars_to_detrend <- names(sub)

# Initialize a list to store residuals
detrended_list <- vector("list", length(vars_to_detrend))
names(detrended_list) <- vars_to_detrend

# Loop through each variable and fit lm(~lon*lat)
for(v in vars_to_detrend){
  lm_model <- lm(as.formula(paste(v, "~ lon + lat")), data = join)
  detrended_list[[v]] <- resid(lm_model)
}

# Combine residuals into a dataframe
detrended_df <- as.data.frame(detrended_list)
cor.mat <- cor(detrended_df, use = "pairwise.complete.obs")
p.mat <- cor_pmat(cor.mat)
ggcorrplot(cor.mat,lab=T,
           type='full',
           insig='blank',
           p.mat = p.mat) 

# Fit a model with all three as predictors for one “dummy” response
# This is just to compute VIF; the response can be any variable
library(usdm)
v <- vifstep(detrended_df,th=5)
ex <- exclude(detrended_df,v)
cor.mat <- cor(ex, use = "pairwise.complete.obs")
p.mat <- cor_pmat(cor.mat)
ggcorrplot(cor.mat,lab=T,
           type='full',
           insig='blank',
           p.mat = p.mat) 


