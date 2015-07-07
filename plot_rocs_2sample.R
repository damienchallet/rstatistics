library(pROC)
source('libdc.R')

createPlot=function(rocs,cbPalette=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),filename=NA,main=""){
  if(!is.na(filename)){
    pdf(filename)
  }
  par(family="serif",font=2)
  pal=cbPalette
  plot(rocs$snrs,t='l',main=main,col=pal[1],lwd=2,cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
  lines(rocs$Rd,t='l',xlim=c(0,1),ylim=c(0,1),col=pal[2],lwd=2)
  lines(rocs$wilcox,t='l',xlim=c(0,1),ylim=c(0,1),col=pal[3],lwd=2)
  lines(rocs$ups,t='l',xlim=c(0,1),ylim=c(0,1),col=pal[4],lwd=2)
  lines(rocs$downs,t='l',xlim=c(0,1),ylim=c(0,1),col=pal[5],lwd=2)
  lines(rocs$Rz,t='l',xlim=c(0,1),ylim=c(0,1),col=pal[6],lwd=2)
  
  legend_text=c("Welch",expression(R[d]),"U",expression(R['+']^'(2)'),expression(R['-']^'(2)'),expression(R[z]))
  legend("bottomright",legend_text,col=pal[1:6],lwd=2,cex=1.5)
  
  if(!is.na(filename)){
    dev.off()
  }
}

load('precomputed/rocs/roc_SNR0.11_sigma1_N100_Navg10000_nPerm10000_nSampH010000_Exponential_nreps1_2sampTRUE_volR1_lenghR1_diffT:records.rda')
allres[[1]]$Rd=allres[[1]]$ups-allres[[1]]$downs
allres[[2]]$Rd=allres[[2]]$ups-allres[[2]]$downs

allres[[1]]$Rz=allres[[1]]$ups_v-allres[[1]]$downs_v
allres[[2]]$Rz=allres[[2]]$ups_v-allres[[2]]$downs_v

whats=c("snrs","wilcox","ups","downs","Rd","Rz")
rocs=lapply(toList(whats),function(what){
  print(what)
  roc(c(rep(0,length(allres[["0"]][[what]])),rep(1,length(allres[["1"]][[what]]))),(-1)^(what=="downs"||what=="wilcox")*c(allres[["0"]][[what]],allres[["1"]][[what]]))
})


createPlot(rocs)
