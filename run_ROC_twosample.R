library(data.table)
library(ggplot2)
source('libH0.R')
source('ggplot_theme_bw_latex.R')

Ns=c(50,100,200)
Ns=100
SNRs=seq(0.11,0.21,0.05)
nus=seq(2.5,5,0.5)
Navg=1000
numSamplesH0=10000
numPerm=100
volRatio=1
lengthRatio=1
diffType="records"
mylapply=mclapply


if(exists("aucs")){
  rm(aucs)
}


for(dtype in c("Student","Exponential","Uniform","Gauss")){
  for(N in Ns){
    for(nu in nus){
      for(SNR in SNRs){
        print(paste(dtype,T,nu,SNR))
        myres=runSampleEstimates(N=N,Navg = Navg,numPerm = numPerm,numSamplesH0 = numSamplesH0,nu=nu,SNR=SNR,dtype=dtype,mylapply=mylapply,forceRun=TRUE,topdf=TRUE,skipComputationMissing=FALSE,twoSamples=TRUE)   
        if(is.null(myres)){
          print("skipping")
          next
        }
        newlines=lapply(names(myres),function(mystat){
          auc.ci=ci(myres[[mystat]])
          newline=data.frame(N=N,dtype=dtype,nu=nu,SNR=SNR,stat=mystat,auc=as.numeric(myres[[mystat]]$auc),auc.ci_upper=auc.ci[3],auc.ci_lower=auc.ci[1],youden=max(myres[[mystat]]$sensitivities+myres[[mystat]]$specificities)-1,Navg=Navg,volRatio=volRatio,diffType=diffType)
          return(newline)
          #          setkeyv(newline,names(newline))
        })
        newlines=as.data.frame(do.call(rbind,newlines))
        if(!exists("aucs")){
          aucs=newlines
        }else{
          aucs=rbind(aucs,newlines)
        }
        if(nrow(aucs)>1){
          save(aucs,file="aucs_twosided.rda")
        }
      }
    }
  }
  if(dtype=="Student"){
    nus=""
  }
}

aucs=as.data.table(aucs)

#Student plots
setkeyv(aucs,c("dtype","N","SNR"))
#SNRs=seq(0.01,0.31,0.05)
for(mySNR in SNRs){
  pdf(paste('auc_vs_nu_Student_N100_SNR',mySNR,'.pdf',sep=""))
  qplot(data=aucs[J("Student",100,mySNR),],nu,auc,color=stat,geom="line",group=stat,main=paste("SNR=",mySNR,sep=""),ylab="AUC",xlab=expression(nu))+theme_bw_latex+scale_color_hue("statistics",labels = c("R+", "R-", "t", "Wilcox","R0"))+geom_linerange(aes(ymin=auc.ci_lower,ymax=auc.ci_upper))
  dev.off()
}

#Gauss plots
setkeyv(aucs,c("dtype","N"))
SNR=0.11
distrs=c("Gauss","Exponential","Uniform")
for(distr in distrs){
  pdf(paste('auc_vs_SNR_',distr,'_N100.pdf',sep=""))
  print(distr)
  qplot(data=aucs[J(distr,100),],SNR,auc,color=stat,geom="line",group=stat,main=distr,ylab="AUC")+theme_bw_latex+scale_color_hue("statistics",labels = c("R+", "R-", "t", "Wilcox","R0"))+geom_linerange(aes(ymin=auc.ci_lower,ymax=auc.ci_upper))
  dev.off()
}
