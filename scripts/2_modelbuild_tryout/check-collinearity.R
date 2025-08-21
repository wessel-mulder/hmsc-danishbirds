# 1. CHECK COLLINEARITY 
input <- '.'
# --> LOAD ENVIRONMENT
X <- read.csv(file.path(input,'data/1_preprocessing/X_environmental/X_Environmental.csv'),row.names=1)
X <- X[sort(row.names(X)),]

# keep only sites with data  
X <- na.omit(X)

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
