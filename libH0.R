source('lib.R')
source('libpermutations.R')
library('pROC')

getTFPN=function(ecdfs,paramRange=seq(-5,5,0.01)){
  FPR=ecdfs[["1"]](paramRange)
  TPR=ecdfs[["0"]](paramRange)
  return(list(FPR=FPR,TPR=TPR))
}

runSampleEstimates=function(N=100,Navg=10000,numSamplesH0=NA,numPerm=1000,nu=3,SNR=0.1,sigma=1,topdf=TRUE,dtype="Gauss",dirSave="precomputed/rocs",mylapply=mclapply,forceRun=FALSE,skipComputationMissing=FALSE,twoSamples=FALSE,nreps=1,volRatio=1,lengthRatio=1,diffType="records",withBEST=FALSE,forceRocRUN=FALSE,withBoot=FALSE){
  
  mymain=dtype
  if(dtype=="Student"){
    mymain=paste(dtype," nu=",nu,sep="")
  }
  filename=paste(dirSave,"/",paste("roc_SNR",SNR,"_sigma",sigma,"_N",N,"_Navg",Navg,"_nPerm",numPerm,"_nSampH0",numSamplesH0,"_",mymain,"_nreps",nreps,"_2samp",twoSamples,"_volR",volRatio,"_lenghR",lengthRatio,"_diffT:",diffType,".rda",sep=""),sep="")
  #   filename_old=paste(dirSave,"/",paste("roc_SNR",SNR,"_sigma",sigma,"_N",N,"_Navg",Navg,"_Nsamples",numPerm,"_",mymain,"_nreps",nreps,"_twosamples",twoSamples,"_volRatio",volRatio,"_lengthRatio",lengthRatio,"_diffType:",diffType,".rda",sep=""),sep="")
  print(filename)
  #   print(filename_old)
  #     if(file.exists(filename_old)){
  #     file.rename(filename_old,filename)
  #   }
  if(!forceRun && file.exists(filename) && file.info(filename)$size>1000){
    load(filename)
  }else{
    if(skipComputationMissing){
      return(NULL)
    }
    factors=c(0,1)
    if(twoSamples){
      allres=lapply(toList(factors),function(myfactor) twosampleEstimates(N=N, Navg=Navg, numPerm=numPerm, numSamplesH0=numSamplesH0, SNR=SNR*myfactor, dtype=dtype,nu=nu,mylapply=mylapply,diffType=diffType,lengthRatio=lengthRatio))       
    }else{
      allres=lapply(toList(factors),function(myfactor) sampleEstimates(N=N, Navg=Navg, numPerm=numPerm, SNR=SNR*myfactor, dtype=dtype,nu=nu,mylapply=mylapply,withBEST=withBEST,withBoot=withBoot))       
    }
    if(!file.exists(dirSave)){
      dir.create(dirSave,recursive = TRUE)
    }
    save(allres,file=filename)
  }
  
  
  filename=paste(dirSave,"/",paste("roc_pROC_SNR",SNR,"_sigma",sigma,"_N",N,"_Navg",Navg,"_nPerm",numPerm,"_nSampH0",numSamplesH0,"_",mymain,"_nreps",nreps,"_2samp",twoSamples,"volR",volRatio,"_lengthR",lengthRatio,"_diffT:",diffType,"_BEST.rda",sep=""),sep="")
  if(!forceRocRUN && file.exists(filename) && file.info(filename)$size>4000){
    load(filename)
  }else{
    if(twoSamples){
      allres[[1]]$Rd=allres[[1]]$ups-allres[[1]]$downs
      allres[[2]]$Rd=allres[[2]]$ups-allres[[2]]$downs

      allres[[1]]$Rz=allres[[1]]$ups_v-allres[[1]]$downs_v
      allres[[2]]$Rz=allres[[2]]$ups_v-allres[[2]]$downs_v
      
      whats=c("rtest2","snrs","wilcox","ups","downs","Rd","Rz")
    }else{
      whats=c("R0s","snrs","wilcox","snr_boot","RT","SIGN")
      if(withBEST){
        whats=c(whats,"BEST")
      }
    }
    rocs=lapply(toList(whats),function(what){
      print(what)
      roc(c(rep(0,length(allres[["0"]][[what]])),rep(1,length(allres[["1"]][[what]]))),(-1)^(what=="downs"||what=="wilcox")*c(allres[["0"]][[what]],allres[["1"]][[what]]))
    })
    
    #    rocs$Rpm= roc(c(rep(0,allres[["0"]]$Navg),rep(1,allres[["1"]]$Navg)),c(allres[["0"]]$ups-allres[["0"]]$downs,allres[["1"]]$ups-allres[["1"]]$downs))
    
    save(rocs,file=filename)
  }
  
  if(topdf){
    if(!file.exists("pdfs/rocs/")){
      dir.create("pdfs/rocs",recursive = TRUE)
    }
    filename=paste("pdfs/rocs/roc_SNR",SNR,"_sigma",sigma,"_N",N,"_Navg",Navg,"_numPerm",numPerm,"_numSamplesH0",numSamplesH0,"_",mymain,"_nreps",nreps,"_twosamples",twoSamples,"volRatio",volRatio,"_lengthRatio",lengthRatio,"_diffType:",diffType,".pdf",sep="")
    pdf(filename)
  }
  if(!twoSamples){
    plot(rocs$snrs,t='l',main=mymain,lwd=2)
    lines(rocs$R0s,t='l',xlim=c(0,1),ylim=c(0,1),col=2,lwd=2)
    lines(rocs$wilcox,t='l',xlim=c(0,1),ylim=c(0,1),col=3,lwd=2)
    #   lines(rocs$rbis,t='l',xlim=c(0,1),ylim=c(0,1),col=6)
    #   browser()
    #  lines(rocs$snr_boot,t='l',xlim=c(0,1),ylim=c(0,1),col=4)
    lines(rocs$RT,t='l',xlim=c(0,1),ylim=c(0,1),col=5,lwd=2)
    lines(rocs$SIGN,t='l',xlim=c(0,1),ylim=c(0,1),col=6,lwd=2)
    nameT="t-stat"
    nameW="Wilcoxon stat"
    legend_text=c(nameT,"r-stat",nameW,"rT-stat","sign stat")
    if(withBEST){
      lines(rocs$BEST,xlim=c(0,1),ylim=c(0,1),col=6)
      legend_text=c(legend_text,"BEST")
      colors=c(colors,6)
    }
    legend("bottomright",legend_text,col=colors,lwd=2)
  }else{
    plot(rocs$snrs,t='l',main=mymain,lwd=2)
    lines(rocs$rtest2,t='l',xlim=c(0,1),ylim=c(0,1),col=2,lwd=2)
    lines(rocs$Rd,t='l',xlim=c(0,1),ylim=c(0,1),col=2,lwd=2)
    lines(rocs$wilcox,t='l',xlim=c(0,1),ylim=c(0,1),col=3,lwd=2)
    lines(rocs$ups,t='l',xlim=c(0,1),ylim=c(0,1),col=4,lwd=2)
    lines(rocs$downs,t='l',xlim=c(0,1),ylim=c(0,1),col=5,lwd=2)
    lines(rocs$Rz,t='l',xlim=c(0,1),ylim=c(0,1),col=6,lwd=2)
    
    legend_text=c("Welch","Rp","Rd","U-test","R+","R_","Rz")
    legend("bottomright",legend_text,col=1:7,lwd=2)
  }    
  
  colors=c(1,2,3,5,6)
  if(twoSamples){
    lines(rocs$Rpm,t='l',xlim=c(0,1),ylim=c(0,1),col=7)
    lines(rocs$rtest2,col=4)
    legend_text=c(legend_text,"r-stat diff","r-stat np")
    colors=c(colors,7,4)
  }  
  legend("bottomright",legend_text,col=colors,lwd=2)
  if(topdf){
    dev.off()
  }
  return(rocs)
}


