for (para in 1:5){
  parameter = c("beta","gamma","rho","V","omega")[para]
  if (parameter=="beta"){
    getpar = function(a)
      return(a$Beta)
    true = all.parameters$beta
    prior = list()
    for(i in 1:npost){
      gamma = matrix(MASS::mvrnorm(n=1, mu=m$mGamma, Sigma = m$UGamma), ncol=m$nt, nrow=m$nc)
      mu = tcrossprod(gamma,m$Tr)
      rho = sample(size=1,x=m$rhopw[,1],prob=m$rhopw[,2])
      V = solve(MCMCpack::rwish(m$f0, solve(m$V0)))
      Si = kronecker(V,rho*all.data$C + (1-rho)*diag(m$ns))
      prior[[i]] = matrix(MASS::mvrnorm(n=1, mu=as.vector(mu), Sigma=Si), ncol=m$ns, nrow=m$nc)
    }
  }
  if (parameter=="rho"){
    getpar = function(a)
      return(a$rho)
    true = all.parameters$rho
    prior = list()
    for(i in 1:npri){
      prior[[i]] = sample(size=1,x=m$rhopw[,1],prob=m$rhopw[,2])
    }
  }
  if (parameter=="gamma"){
    getpar = function(a)
      return(a$Gamma)
    true = all.parameters$gamma
    prior = list()
    for(i in 1:npri){
      prior[[i]] = matrix(MASS::mvrnorm(n=1, mu=m$mGamma, Sigma = m$UGamma), ncol=m$nt, nrow=m$nc)
    }
  }
  
  if (parameter=="V"){
    getpar = function(a)
      return(a$V)
    true = all.parameters$V
    prior = list()
    for(i in 1:npost){
      prior[[i]] = solve(MCMCpack::rwish(m$f0, solve(m$V0)))
    }
  }
  if (parameter=="omega"){
    getpar = function(a)
      return(t(a$Lambda[[1]])%*%a$Lambda[[1]])
    true = t(all.parameters$lambda)%*%all.parameters$lambda
    nu = m$rL[[1]]$nu
    a1 = m$rL[[1]]$a1
    a2 = m$rL[[1]]$a2
    b1 = m$rL[[1]]$b1
    b2 = m$rL[[1]]$b2
    prior = list()
    for(i in 1:npost){
      nf = 10
      delta = rep(NA,nf)
      delta[1] = rgamma(1,a1,b1)
      for (j in 2:nf){
        delta[j] = rgamma(1,a2,b2)
      }
      tau = cumprod(delta)
      psi = matrix(rgamma(nf*m$ns, nu/2, nu/2), nf, m$ns)
      sd = sqrt(1/(psi*matrix(rep(tau,m$ns),nr=nf, nc=m$ns)))
      la = rnorm(length(sd), mean=0, sd=sd)
      dim(la) = dim(sd)
      prior[[i]] = t(la)%*%la
    }
  }
  post = lapply(postList, getpar)
  post.mean = Reduce("+", post) / length(post)
  cors = rep(0,npost)
  for (i in 1:npost){
    cors[i] = cor(as.vector(true),as.vector(post[[i]]))
  }
  se = function(a)
    return((a-true)^2)
  post.se = lapply(post, se)
  prior.se = lapply(prior, se)
  post.rmse = sqrt(Reduce("+", post.se) / length(post.se))
  prior.rmse = sqrt(Reduce("+", prior.se) / length(prior.se))
  nrmse = post.rmse / prior.rmse
  
  if (parameter=="beta"){
    estimates$cor.beta = cors
    estimates$nrmse.beta = nrmse
    
  }
  if (parameter=="gamma"){
    estimates$cor.gamma = cors
    estimates$nrmse.gamma = nrmse
  }
  if (parameter=="rho"){
    estimates$post.rho = post
    estimates$nrmse.rho = nrmse
  }
  if (parameter=="V"){
    estimates$cor.V = cors
    estimates$nrmse.V = nrmse
  }
  if (parameter=="omega"){
    estimates$cor.omega = cors
    estimates$nrmse.omega = nrmse
  }
}

save(file=file.path(localDir, "performance",paste("Case",toString(case),".R",sep="")),
     computational.time, mixing, estimates, MF)