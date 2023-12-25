library('TOSTER')
# calculating the equivalence bounds by power and sample size
# powerTOSTone(alpha=0.05,N=200,statistical_power=.8) # d=+-0.21
# powerTOSTone(alpha=0.05,N=29,statistical_power=.8) # d=+-0.54

nsub = 29

setwd('D:\\Matlab_Workspace\\tom_NewStory\\modelling\\revision_cogpsy\\model_run')
for (irun in c(1,2)){
  tmp = vector(length=9801)
  for (ifile in seq(1,9801)){
    f = paste('.\\data_mi_run',as.character(irun),'_newexp\\data',as.character(ifile),'.txt',sep='')
    d = read.table(f)
    delta = c(d[1:nsub,1] - d[1:nsub,2])
    m_d = mean(delta)
    s_d = sd(delta)
    obj = TOSTone(m=m_d,mu=0,sd=s_d,n=nsub,low_eqbound_d=-0.54,high_eqbound_d=0.54,alpha=0.05)
    tmp[ifile] = obj$TOST_p2
    dev.off()
    cat('\014')
  }
  write.table(tmp,file=paste('.\\R_run',as.character(irun),'_',as.character(nsub),'.txt',sep=''),sep=' ',row.names=FALSE,col.names=FALSE)
}
