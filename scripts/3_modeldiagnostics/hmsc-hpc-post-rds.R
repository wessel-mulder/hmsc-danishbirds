input <- '/home/bhr597/home/projects/hmsc-danishbirds/tmp_rds/2025-05-22_16-39-28'

# GETTING STARTED ---------------------------------------------------------
library(RColorBrewer,lib="~/Rlibs")
library(farver,lib="~/Rlibs")
library(scales,lib="~/Rlibs")
library(jsonify,lib="~/Rlibs")
library(ape,lib="~/Rlibs")
library(dplyr,lib="~/Rlibs")
library(Hmsc,lib="~/Rlibs")

# GETTING STARTED --------------------------------------------------------
# load unfitted object
m <- readRDS(file.path(input,'m_object.rds'))
params <- readRDS(file.path(input,'params.rds'))


# loading the chains 
chainList = vector("list", nChains)
for(cInd in 1:nChains){
    chain_file_path = file.path(input, sprintf("post_chain%.2d_file.rds", cInd-1))
    chainList[[cInd]] = from_json(readRDS(file = chain_file_path)[[1]])[[1]]
}

fitSepTF = importPosteriorFromHPC(m, chainList, 100, 10, 1000)


