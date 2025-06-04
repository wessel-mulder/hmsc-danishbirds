rm(list = ls())
# GETTING STARTED ---------------------------------------------------------
python = file.path(getwd(), 'hmsc-hpc-main',"hmsc-venv", "bin", "python3.11")

# get python running 
package_path = file.path('hmsc-hpc-main/')
system2(python, "-m pip install --upgrade pip")
system2(python, paste("-m pip install", shQuote(package_path)))

# check tensorflow 
Sys.setenv(TF_CPP_MIN_LOG_LEVEL=3)
system2(python, "-c \"import tensorflow as tf; print(tf.constant(1))\"")
system2(python, "-c \"import hmsc\"")

library(jsonify)

# Or from a local source folder
library(devtools)
devtools::install_local("HMSC-master/",force=T)
library(Hmsc)

# DEFINING PARAMETERS -----------------------------------------------------
nSamples=100
thin=2
nChains=4
verbose=100
transient=nSamples*thin
summary(TD)

# init hmsc object
m = Hmsc(Y=TD$Y,
         XData=TD$X, XFormula=~.,
         TrData=TD$Tr[,-1], TrFormula=~.,
         phyloTree=TD$phy,
         studyDesign=TD$studyDesign,
         ranLevels=list(plot=TD$rL1, sample=TD$rL2))

# initiate mcmc sampling 
init_obj = sampleMcmc(m, samples=nSamples, thin=thin,
                      transient=transient, nChains=nChains,
                      verbose=verbose, engine="HPC")
# write to path
init_file_path = file.path(getwd(),'hmsc-hpc-main', "init_file.rds")
# save as json 
saveRDS(to_json(init_obj), file=init_file_path)

# operates in python, so formulate the required call 
post_file_path = file.path(getwd(),'hmsc-hpc-main', "post_file.rds")
python_cmd_args = paste("-m hmsc.run_gibbs_sampler",
                        "--input", shQuote(init_file_path),
                        "--output", shQuote(post_file_path),
                        "--samples", nSamples,
                        "--transient", transient,
                        "--thin", thin,
                        "--verbose", verbose)
cat(paste(shQuote(python), python_cmd_args), "\n")

# running this actuall script in python
system2(python, python_cmd_args)


# POSTERIORS --------------------------------------------------------------
importFromHPC = from_json(readRDS(file = post_file_path)[[1]])
postList = importFromHPC[1:nChains]
cat(sprintf("fitting time %.1f sec\n", importFromHPC[[nChains+1]]))

fitTF = importPosteriorFromHPC(m, postList, nSamples, thin, transient)
plotVariancePartitioning(fitTF, computeVariancePartitioning(fitTF), args.legend=list(x="bottomright"))


# HPC FOR PARALLEL --------------------------------------------------------
chain_cmd_args_list = vector("list", nChains)
for(cInd in 1:nChains){
  chain_file_path = file.path(getwd(),'hmsc-hpc-main', sprintf("post_chain%.2d_file.rds", cInd-1))
  chain_cmd_args = paste("-m hmsc.run_gibbs_sampler",
                         "--input", shQuote(init_file_path),
                         "--output", shQuote(chain_file_path),
                         "--samples", nSamples,
                         "--transient", transient,
                         "--thin", thin,
                         "--verbose", verbose,
                         "--chain", cInd-1)
  cat(paste(shQuote(python), chain_cmd_args), "\n")
  chain_cmd_args_list[[cInd]] = chain_cmd_args
}

for(cInd in 1:nChains){
  system2(python, chain_cmd_args_list[[cInd]])
}
# load posteriors into R
chainList = vector("list", nChains)
for(cInd in 1:nChains){
  chain_file_path = file.path(getwd(),'hmsc-hpc-main', sprintf("post_chain%.2d_file.rds", cInd-1))
  chainList[[cInd]] = from_json(readRDS(file = chain_file_path)[[1]])[[1]]
}

fitSepTF = importPosteriorFromHPC(m, chainList, nSamples, thin, transient)
plotVariancePartitioning(fitSepTF, computeVariancePartitioning(fitSepTF), args.legend=list(x="bottomright"))
