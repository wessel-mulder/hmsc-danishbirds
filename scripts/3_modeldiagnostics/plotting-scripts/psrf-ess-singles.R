# init output documents 


for(i in seq_along(params)){
  ### PSRFS 
  if(i == 1){
    if(!dir.exists(file.path(input,'results'))) {dir.create(file.path(input,'results'))}
    pdf(file=file.path(input,'results','PSRF_ESS.pdf'),
        width = 8,
        height = 6)
  }
  
  # plot PSRF
  par(mfrow=c(1,2))
  diags$psrf[[i]][diags$psrf[[i]] > 10] <- 10
  vioplot(diags$psrf[[i]],col=cols[i],ylim=c(1,max(diags$psrf[[i]],na.rm = T)),
          main=paste0('PSRF: ',full_names[i]))
  abline(h=1.1,col='red')
  vioplot(diags$psrf[[i]],col=cols[i],ylim=c(1,1.1),main=paste0('Mean = ',round(mean(diags$psrf[[i]],na.rm = T),2)))
  abline(h=1.1,col='red')
  par(mfrow=c(1,1))
  hist(diags$ess[[i]],breaks = seq(0,max(diags$ess[[i]])+100,by=100),
       xlab = 'Number of effective samples',
       main=paste0('ESS: ',full_names[i]))
  abline(v=100,col='red',lty=1)
  abline(v=200,col='red',lty=2)
  abline(v=1000,col='red',lty=3)
  legend(
    "topright",
    legend = c(
      paste0("Proportion of effective samples >100 = ",
             round(length(diags$ess[[i]][diags$ess[[i]] > 100]) / length(diags$ess[[i]]), 2)),
      paste0("Proportion of effective samples >200 = ",
             round(length(diags$ess[[i]][diags$ess[[i]] > 200]) / length(diags$ess[[i]]), 2)),
      paste0("Proportion of effective samples >1000 = ",
             round(length(diags$ess[[i]][diags$ess[[i]] > 1000]) / length(diags$ess[[i]]), 2))
    ),
    col = "red",
    lty = c(1, 2, 3)
  )
  
  if(i == length(params)){dev.off()}
}
