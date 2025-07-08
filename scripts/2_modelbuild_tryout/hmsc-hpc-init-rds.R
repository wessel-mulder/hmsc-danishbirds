rm(list = ls())
input <- '~/home/projects/hmsc-danishbirds'
python <- '~/opt/software/anaconda3/2021.05/bin/python'
# Define MCMC settings
nChains <- 4
thin <- 10
nSamples <- 2000
transient <- thin*nSamples
verbose <- 100
params <- list(
    nChains = nChains,
    thin = thin,
    nSamples = nSamples,
    transient = transient,
    verbose = verbose
)
# GETTING STARTED ---------------------------------------------------------
if (interactive() && Sys.getenv("RSTUDIO") == "1") {
  message("Running in RStudio")
  library(Hmsc)
  library(jsonify)
  library(knitr)
  library(corrplot)
} else {
  message("Running from terminal or non-interactive environment")
  library(RColorBrewer,lib="~/Rlibs")
  library(farver,lib="~/Rlibs")
  library(scales,lib="~/Rlibs")
  library(jsonify,lib="~/Rlibs")
  library(ape,lib="~/Rlibs")
  library(dplyr,lib="~/Rlibs")
  library(Hmsc,lib="~/Rlibs")
}

summary(TD)

# LOADING DATA -----------------------------------------------------

# --> LOAD ENVIRONMENT
X <- read.csv(file.path(input,'data/1_preprocessing/X_environmental/X_Environmental.csv'),row.names=1)
X <- X[sort(row.names(X)),]

# keep only sites with data  
X <- na.omit(X)
sites_actual <- row.names(X)

# --> LOAD OCCURRENCES 
Y <- read.csv(file.path(input,'data/1_preprocessing/Y_occurrences/Y_occurrences.csv'),row.names=1)
# merge subspecies for now 
# Merge Motacilla_alba_yarrellii into Motacilla_alba
Y$Motacilla_alba <- pmax(Y$Motacilla_alba, Y$Motacilla_alba_yarrellii, na.rm = TRUE)
Y$Motacilla_alba_yarrellii <- NULL  # Remove the subspecies column

# Merge Motacilla_flava_flavissima into Motacilla_flava
Y$Motacilla_flava <- pmax(Y$Motacilla_flava, Y$Motacilla_flava_flavissima, na.rm = TRUE)
Y$Motacilla_flava_flavissima <- NULL  # Remove the subspecies column

# rename acanthis to species name
names(Y)[names(Y) == "Acanthis_flammea_cabaret"] <- "Acanthis_flammea"
Y <- Y[sort(row.names(Y)), sort(colnames(Y))]

# remove sites without data 
Y <- Y[row.names(Y) %in% sites_actual,]

# --> LOAD TRAITS 
Tr <- read.csv(file.path(input,"data/1_preprocessing/Tr_aits/traits-guild_migration.csv"),row.names = 2)[,c(2,3)]
Tr <- Tr[sort(row.names(Tr)),]
Tr['Motacilla_alba_yarrellii',] == Tr['Motacilla_alba',] # TRUE 
Tr['Motacilla_flava_flavissima',] == Tr['Motacilla_flava',] # MIGRATION FALSE 
Tr <- Tr[!(rownames(Tr) %in% c('Motacilla_alba_yarrellii','Motacilla_flava_flavissima')), ]
rownames(Tr)[rownames(Tr) == "Acanthis_flammea_cabaret"] <- "Acanthis_flammea"

# --> LOAD STUDY DESIGN 
Design <- read.csv(file.path(input,"data/1_preprocessing/design/studyDesign.csv"),row.names=5)

# remove sites without data 
Design <- Design[row.names(Design) %in% sites_actual,]

# sort
Design <- Design[sort(row.names(Design)), ]

# convert to factors
Design$site <- as.factor(Design$site)
Design$atlas <- as.factor(Design$atlas)

# phylogeny
phy <- read.tree(file.path(input,'data/1_preprocessing/Taxonomy/tree_fromPD.tre'))
# check pd
pd_matrix <- cophenetic.phylo(phy)
pd_matrix <- pd_matrix[sort(rownames(pd_matrix)), sort(colnames(pd_matrix))]

# check if species-lists are identical
table(rownames(pd_matrix) == names(Y))
# stunning 



# get xycoords
xycoords <- Design[,colnames(Design) %in% c('lat','lon','site')]
xycoords <- distinct(xycoords)
rownames(xycoords) <- xycoords$site
xycoords <- xycoords[,colnames(xycoords) %in% c('lat','lon')]

# PREPARING MODEL BUILD ---------------------------------------------------
# Define random effects
xycoords <- as.matrix(xycoords)
rand_site <- HmscRandomLevel(sData = xycoords)
rand_atlas <- HmscRandomLevel(units = levels(Design$atlas))

# Define model formulas for environmental and trait data
XFormula <- as.formula(paste("~", paste(colnames(X), collapse = "+"), sep = " "))
TrFormula <- as.formula(paste("~", paste(colnames(Tr), collapse = "+"), sep = " "))

# define m 
m <-Hmsc(Y = Y, XData = X, XFormula = XFormula, TrData = Tr,
         TrFormula = TrFormula, phyloTree = phy,
         studyDesign = Design[,c(1,4)], 
         ranLevels = list('site'=rand_site,'atlas'=rand_atlas),
         distr='probit')

# initiate mcmc sampling 
init_obj <- sampleMcmc(m, samples=nSamples, thin=thin,
                      transient=transient, nChains=nChains,
                      verbose = verbose,
                      engine="HPC")

dir_name <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
dir.create(file.path(input,'tmp_rds',dir_name))

init_file_path = file.path(input,'tmp_rds',dir_name, "init_file.rds")
m_file_path = file.path(input,'tmp_rds',dir_name, "m_object.rds")
param_file_path = file.path(input,'tmp_rds',dir_name, "params.rds")
lines <- paste(names(params), unlist(params), sep = "=")
writeLines(lines, file.path(input,'tmp_rds',dir_name, "params.txt"))


# save as json 
saveRDS(to_json(init_obj), file=init_file_path)
saveRDS(m,file=m_file_path)
saveRDS(params,file=param_file_path)


# operates in python, so formulate the required call 
post_file_path = file.path(input,'tmp_rds',dir_name, "post_file.rds")
python_cmd_args = paste("-m hmsc.run_gibbs_sampler",
                        "--input", shQuote(init_file_path),
                        "--output", shQuote(post_file_path),
                        "--samples", nSamples,
                        "--transient", transient,
                        "--thin", thin,
                        "--verbose", verbose)
cat(paste(shQuote(python), python_cmd_args), "\n")

# running this actuall script in python
#system2(python, python_cmd_args)


# GETTING POSTERIORS ------------------------------------------------------
#importFromHPC = from_json(readRDS(file = post_file_path)[[1]])
##postList = importFromHPC[1:nChains]
#cat(sprintf("fitting time %.1f sec\n", importFromHPC[[nChains+1]]))

#fitTF = importPosteriorFromHPC(m, postList, nSamples, thin, transient)
#plotVariancePartitioning(fitTF, computeVariancePartitioning(fitTF), args.legend=list(x="bottomright"))




