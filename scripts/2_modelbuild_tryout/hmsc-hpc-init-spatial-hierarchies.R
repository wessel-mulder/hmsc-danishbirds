rm(list = ls())

# Define MCMC settings
subset_env_vars <- 0 # flip to 0 to reverse subset
knotDistance = 0.1 # knotdistances 
nChains <- 4
verbose <- 100

thin <- c(10)
nSamples <- c(250)
transient <- 5000

# GETTING STARTED ---------------------------------------------------------
if (interactive() && Sys.getenv("RSTUDIO") == "1") {
  message("Running in RStudio")
  library(Hmsc)
  library(jsonify)
  library(knitr)
  library(corrplot)
  library(ape)
  library(dplyr)
  input <- '.'
  python <- file.path(getwd(), 'hmsc-hpc-main',"hmsc-venv", "bin", "python3.11")
  flagInit = 0
  flagFitR = 1
} else {
  message("Running from terminal or non-interactive environment")
  library(RColorBrewer,lib="~/Rlibs")
  library(farver,lib="~/Rlibs")
  library(scales,lib="~/Rlibs")
  library(jsonify,lib="~/Rlibs")
  library(ape,lib="~/Rlibs")
  library(dplyr,lib="~/Rlibs")
  library(Hmsc,lib="~/Rlibs")
  library(dplyr,lib="~/Rlibs")
  library(withr,lib="~/Rlibs")
  input <- '~/home/projects/hmsc-danishbirds'
  python <- '~/opt/software/anaconda3/2021.05/bin/python'
  flagInit = 1
  flagFitR = 0
}

# check if accurately installed
summary(TD)

# LOADING DATA -----------------------------------------------------
# --> LOAD ENVIRONMENT
X <- read.csv(file.path(input,'data/1_preprocessing/X_environmental/X_Environmental.csv'),row.names=1)
X <- X[sort(row.names(X)),]

# keep only sites with data  
sites_with_NA <- rownames(X)[apply(X, 1, function(x) any(is.na(x)))]
X[rownames(X) %in% sites_with_NA,]

# rows <- X[rownames(X) %in% sites_with_NA,]
# 
# # for each row, get the column names where value == 1
# cols_by_row <- lapply(sites_with_NA, function(r) {
#   colnames(Y)[Y[r, ] == 1]
# })
# 
# # name the list by row
# names(cols_by_row) <- sites_with_NA
# 
# cols_by_row
# cols_by_row$
# 
X <- na.omit(X)
sites_actual <- row.names(X)

# --> LOAD OCCURRENCES 
Y <- read.csv(file.path(input,'data/1_preprocessing/Y_occurrences/Y_occurrences.csv'),row.names=1)

# remove sites without data 
Y <- Y[row.names(Y) %in% sites_actual,]

# --> LOAD TRAITS 
Tr <- read.csv(file.path(input,"data/1_preprocessing/Tr_aits/traits-guild_migration.csv"),row.names = 2)[,c(2,3)]
# sort by Y
Tr <- Tr[colnames(Y), , drop = FALSE]

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


# phylogeny
phy <- read.tree(file.path(input,'data/1_preprocessing/Taxonomy/tree_fromPD.tre'))
# check pd
pd_matrix <- cophenetic.phylo(phy)
pd_matrix <- pd_matrix[sort(rownames(pd_matrix)), sort(colnames(pd_matrix))]

# check if species-lists are identical
setdiff(rownames(pd_matrix),names(Y))
# stunning 

# number of occurrences, make sure there are no below 5 
info <- apply(Y,2,sum)
table(info)

# make sure they line up
all(rownames(Tr) == colnames(Y))



# PREPARING MODEL BUILD ---------------------------------------------------
# Define model types: 

date <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
for(i in thin){
  print(i)
  
  for(j in nSamples){
    print(j)

    params <- list(
    nChains = nChains,
    thin = i,
    nSamples = j,
    transient = transient,
    verbose = verbose,
    knotDist = knotDistance
    )



  XFormula <- as.formula(paste("~", paste(colnames(X), collapse = "+"), sep = " "))
  TrFormula <- as.formula(paste("~", paste(colnames(Tr), collapse = "+"), sep = " "))
  
  # get year sorted 
  Design$year <- as.numeric(as.character(Design$year))
  time <- data.frame(year = sort(unique(Design$year)))
  rownames(time) <- sort(unique(Design$year))
  struc_time <- HmscRandomLevel(sData = time)
  Design$year <- as.factor(Design$year)
  
  # get space
  xycoords <- Design[,colnames(Design) %in% c('lat','lon','site')]
  xycoords <- distinct(xycoords)
  rownames(xycoords) <- xycoords$site
  xycoords <- xycoords[,colnames(xycoords) %in% c('lat','lon')]

  # get blocks and midpoints 
Design$block <- substr(Design$site, 1, 2)
Design$block <- as.factor(Design$block)

midpoints <- lapply(levels(Design$block), function(block) {
  x <- Design$lon[Design$block == block]
  y <- Design$lat[Design$block == block]
  
  xmid <- (max(x) + min(x)) / 2
  ymid <- (max(y) + min(y)) / 2
  
  data.frame(block = block, lon = xmid, lat = ymid)
})

midpoints <- do.call(rbind, midpoints)
rownames(midpoints) <- midpoints$block
midpoints <- midpoints[,colnames(midpoints) %in% c('lat','lon')]
  
  # construct nots 
  xycoords <- as.matrix(xycoords)
  xyKnots = constructKnots(xycoords,knotDist = knotDistance, minKnotDist = knotDistance)
  nKnots <- nrow(xyKnots)
  plot(xycoords[,2],xycoords[,1],pch=18, asp=1)
  points(xyKnots[,2],xyKnots[,1],col='red',pch=18)
  
  # smallest level 
  struc_space <- HmscRandomLevel(sData = xycoords, sMethod = "GPP",
                                 sKnot = xyKnots)
  struc_block <- HmscRandomLevel(sData = midpoints, sMethod = "GPP",
                                 sKnot = xyKnots)
  struc_space$alphapw
  # standard frequencies 
  freq <- c(0.5,rep(0.005,100))
  #under 1 degree
  samples <- c(0,seq(from = 0.05, to = 0.1, length.out = 100))
  smallest <- cbind(samples,freq)

  samples <- c(0,seq(from = 0.1, to = 0.5, length.out = 100))
  small <- cbind(samples,freq)
  
  #under 1-3 degree
  samples <- c(0,seq(from = 0.5, to = 1, length.out = 100))
  medium <- cbind(samples,freq)
  
  # 3 - max degree
  samples <- c(0,seq(from = 1, to = 7.7, length.out = 100))
  large <- cbind(samples,freq)
  
  struc_space_smallest <- setPriors(struc_space,alphapw=smallest,nfMax=5)
  struc_space_small <- setPriors(struc_space,alphapw=small,nfMax=5)

  struc_space_medium <- setPriors(struc_block,alphapw=medium,nfMax=5)
  struc_space_large <- setPriors(struc_block,alphapw=large,nfMax=5)
  
  Design$blockregion <- Design$block
  Design$siteregion <- Design$site

  # define m 
  m <-Hmsc(Y = Y, XData = X, XFormula = XFormula, TrData = Tr,
           TrFormula = TrFormula, phyloTree = phy,
           studyDesign = Design[,c('year','site','siteregion','block','blockregion')], 
           ranLevels = list('site'=struc_space_smallest,'siteregion'=struc_space_small,'block'=struc_space_medium,'blockregion'=struc_space_large,
                            'year'=struc_time),
           distr='probit')
  
  
  rownames(phy)
  
  ### IN RSTUDIO START SAMPLING 
  if(flagFitR){
    print('Rstudio stuff executed')
  }
  ### IN HPC ENVIORNMENT SET UP INIT
  if(flagInit){
  # initiate mcmc sampling 
  init_obj <- sampleMcmc(m, samples=j, thin=i, # flag 
                         transient=transient, nChains=nChains,
                         verbose = verbose,
                         engine="HPC")
  
  dir_name <- paste0(date,'_spatialhierarchies')
  print(dir_name)
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
  print('Init files created')
  }
  }
}

