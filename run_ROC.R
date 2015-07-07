library(data.table)
library(ggplot2)
library(pROC)

library(doMC)
registerDoMC(cores = 4)

source('libH0.R')
source('ggplot_theme_bw_latex.R')

#Ns=c(10,20,30,50,100,200)
Ns=100
numSamplesH0=50000
numPerm=1000
Navg=1000
SNRs=seq(0.01,0.21,0.025)
nus=seq(2,5,0.5)
forceRun=FALSE
topdf=FALSE
skipComputationMissing=FALSE
nreps=1
twoSamples=FALSE
withBoot=TRUE
mylapply=mclapply

if(exists("aucs")){
  rm(aucs)
}

for(dtype in c("Student","Exponential","Uniform","Gauss","Cauchy")){
  for(volRatio in 1){ #c(0.1,0.3,1,3,10)){
    for(lengthRatio in 1){  # c(0.1,0.3,1,3,10)){
      for(diffType in "records"){  #},"vectors")){
        for(Nsamples in 1000){      #c(10,31,100,310,1000,3100)){
          for(N in Ns){
            for(nu in nus){
              for(SNR in SNRs){
                print(paste(dtype,N,nu,SNR))
                myres=runSampleEstimates(N=N,Navg = Navg,numSamplesH0=numSamplesH0,numPerm = numPerm,nu=nu,SNR=SNR,dtype=dtype,mylapply=mylapply,forceRun=forceRun,topdf=topdf,skipComputationMissing=skipComputationMissing,nreps=nreps,twoSamples = twoSamples,volRatio=volRatio,lengthRatio=lengthRatio,diffType=diffType,withBEST=FALSE,forceRocRUN=TRUE,withBoot=withBoot)   
                if(is.null(myres)){
                  print("  skipping")
                  next
                }
                if(!is.null(myres$up) && !is.null(myres$wilcox)){
                  mytest=roc.test(myres$ups,myres$wilcox)
                }else{
                  mytest=list(Z=NA,p.value=NA)
                }
                newlines=mylapply(names(myres),function(mystat){
                  auc.ci=ci(myres[[mystat]])
                  newline=data.frame(N=N,dtype=dtype,nu=nu,SNR=SNR,stat=mystat,auc=as.numeric(myres[[mystat]]$auc),auc.ci_upper=auc.ci[3],auc.ci_lower=auc.ci[1],youden=max(myres[[mystat]]$sensitivities+myres[[mystat]]$specificities)-1,Nsamples=Nsamples,volRatio=volRatio,diffType=diffType,Zdelong=mytest$Z,pval.delong=mytest$p.value)
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
                  save(aucs,file="aucs_onesided.rda")
                }
                
              }
            }
          }
        }
      }
    }
    if(dtype=="Student"){
      nus=""
    }
  }
}

aucs=as.data.table(aucs)
save(aucs,file="aucs.rda")

browser()

toPdf=TRUE
toPdf=FALSE

#Student plots
setkeyv(aucs,c("dtype","N","SNR","Nsamples"))
SNRs=sort(unique(aucs$SNR))
for(mySNR in SNRs){
  if(toPdf){
    pdf(paste('auc_vs_nu_Student_N',N,"_Navg",Navg,"Nsamples",Nsamples,"_twoSided=",twoSided,"_volRatio",volRatio,"_diffType",diffType,'_SNR',mySNR,'.pdf',sep=""))
  }
  print(qplot(data=aucs[J("Student",100,mySNR,1000),],nu,auc,color=stat,geom="line",group=stat,main=paste("SNR=",mySNR,sep=""),ylab="AUC",xlab=expression(nu))+theme_bw_latex+scale_color_hue("statistics",labels = c("R+", "R-", "t", "Wilcox","R0"))+geom_linerange(aes(ymin=auc.ci_lower,ymax=auc.ci_upper)))
  if(toPdf){
    dev.off()
  }
}

#Student plots
setkeyv(aucs,c("dtype","N","SNR","diffType"))
SNRs=sort(unique(aucs$SNR))
for(mySNR in SNRs){
  if(toPdf){
    mySNRstr=sub("\\.","",as.character(mySNR),perl=TRUE)
    pdf(paste('youden_vs_nu_Student_N',N,"_Navg",Navg,"Nsamples",Nsamples,"_twoSided=",twoSided,"_volRatio",volRatio,"_diffType",diffType,'_SNR',mySNRstr,'.pdf',sep=""))
  }
  print(qplot(data=aucs[J("Student",100,mySNR,1000),],nu,youden,color=stat,geom="line",group=stat,main=paste("SNR=",mySNR,sep=""),ylab="AUC",xlab=expression(nu))+theme_bw_latex+scale_color_hue("statistics",labels = c("R+", "R-", "t", "Wilcox","R0")))
  if(toPdf){
    dev.off()
  }
}



#Gauss plots
setkeyv(aucs,c("dtype","N","Nsamples"))
distrs=c("Gauss","Exponential","Uniform")
for(distr in distrs){
  if(toPdf){
    mySNRstr=sub("\\.","",as.character(mySNR),perl=TRUE)
    pdf(paste('auc_vs_SNR_',distr,'_N',N,"_Navg",Navg,"Nsamples",Nsamples,"_twoSided=",twoSided,"_volRatio",volRatio,"_diffType",diffType,'_SNR',mySNRstr,'.pdf',sep=""))
  }
  print(distr)
  print(qplot(data=aucs[J(distr,100),],SNR,auc,color=stat,geom="line",group=stat,main=distr,ylab="AUC")+theme_bw_latex+scale_color_hue("statistics",labels = c("R+", "R-", "t", "Wilcox","R0"))+geom_linerange(aes(ymin=auc.ci_lower,ymax=auc.ci_upper)))
  if(toPdf){
    dev.off()
  }
}

#Gauss plots
setkeyv(aucs,c("dtype","N","Nsamples"))
distrs=c("Gauss","Exponential","Uniform")
for(distr in distrs){
  if(toPdf){
    mySNRstr=sub("\\.","",as.character(mySNR),perl=T)
    pdf(paste('youden_vs_SNR_',distr,'_N',N,"_Navg",Navg,"Nsamples",Nsamples,"_twoSided=",twoSided,"_volRatio",volRatio,"_diffType",diffType,'_SNR',mySNRstr,'.pdf',sep=""))
  }
  print(distr)
  print(qplot(data=aucs[J(distr,100),],SNR,youden,color=stat,geom="line",group=stat,main=distr,ylab="AUC")+theme_bw_latex+scale_color_hue("statistics",labels = c("R+", "R-", "t", "Wilcox","R0")))
  if(toPdf){
    dev.off()
  }
}



setkeyv(aucs,c("dtype","N","nu"))
distr="Student"
nu="3"
if(toPdf){
  pdf(paste('auc_vs_SNR_',distr,'_nu',nu,'_N',N,"_Navg",Navg,"Nsamples",Nsamples,"_twoSided=",twoSided,"_volRatio",volRatio,"_diffType",diffType,'.pdf',sep=""))
}
print(qplot(data=aucs[J(distr,100,"3"),],SNR,auc,color=stat,geom="line",group=stat,main=paste(distr," nu=",nu,sep=""),se,ylab="AUC")+theme_bw_latex+scale_color_hue("statistics",labels = c("R+", "R-", "t", "Wilcox","R0"))+geom_linerange(aes(ymin=auc.ci_lower,ymax=auc.ci_upper)))
if(toPdf){
  dev.off()
}

setkeyv(aucs,c("dtype","N","nu"))
distr="Student"
nu="3"
if(toPdf){
  pdf(paste('youden_vs_SNR_',distr,'_nu',nu,'_N',N,"_Navg",Navg,"Nsamples",Nsamples,"_twoSided=",twoSided,"_volRatio",volRatio,"_diffType",diffType,'_SNR',mySNR,'.pdf',sep=""))
}
print(qplot(data=aucs[J(distr,100,"3"),],SNR,youden,color=stat,geom="line",group=stat,main=paste(distr," nu=",nu,sep=""),se,ylab="AUC")+theme_bw_latex+scale_color_hue("statistics",labels = c("R+", "R-", "t", "Wilcox","R0")))
if(toPdf){
  dev.off()
}


setkeyv(aucs,c("dtype","N","nu"))
distr="Student"
nu="3"
if(toPdf){
  pdf(paste('youden_vs_SNR_',distr,'_nu',nu,'_N',N,"_Navg",Navg,"Nsamples",Nsamples,"_twoSided=",twoSided,"_volRatio",volRatio,"_diffType",diffType,'_SNR',mySNR,'.pdf',sep=""))
}
print(qplot(data=aucs[J(distr,100,"3"),],SNR,youden,color=stat,geom="line",group=stat,main=paste(distr," nu=",nu,sep=""),se,ylab="AUC")+theme_bw_latex+scale_color_hue("statistics",labels = c("R+", "R-", "t", "Wilcox","R0")))
if(toPdf){
  dev.off()
}


