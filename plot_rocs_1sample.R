library(pROC)

createPlot=function(rocs,cbPalette=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),filename=NA,main=""){
  if(!is.na(filename)){
    pdf(filename)
  }
  par(family="serif",font=2)
  plot(rocs$R0s,t='l',main=main,lwd=2,col=cbPalette[1],cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
  lines(rocs$SIGN,t='l',xlim=c(0,1),ylim=c(0,1),col=cbPalette[2],lwd=2)
  lines(rocs$snrs,t='l',xlim=c(0,1),ylim=c(0,1),col=cbPalette[3],lwd=2)
  lines(rocs$wilcox,t='l',xlim=c(0,1),ylim=c(0,1),col=cbPalette[4],lwd=2)
  legend("bottomright",c("r","sign","t","Wilcoxon"),col=cbPalette,lwd=2,title="statistics",cex=1.5)
  if(!is.na(filename)){
    dev.off()
  }
}

load('precomputed/rocs/roc_SNR0.11_sigma1_N100_Navg10000_nPerm10000_nSampH050000_Gauss_nreps1_2sampFALSE_volR1_lenghR1_diffT:records.rda')
whats=c("R0s","snrs","wilcox","SIGN")
rocs=lapply(toList(whats),function(what){
  print(what)
  roc(c(rep(0,length(allres[["0"]][[what]])),rep(1,length(allres[["1"]][[what]]))),(-1)^(what=="downs"||what=="wilcox")*c(allres[["0"]][[what]],allres[["1"]][[what]]))
})
