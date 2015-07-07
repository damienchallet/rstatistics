source('libdc.R')
source('libpermutations.R')

N=100
Navg=10000
Nsamples=1000
SNRmin=0.01
SNRs=SNRmin*(10^(1/15))^c(0:15)  #from 0.01 to 0.1
GaussS=rev(c(TRUE,FALSE))

forceRUN=TRUE
dirSave="precomputed/efficiency"
dir.create(dirSave,recursive = TRUE,showWarnings = FALSE)

for(Gauss in GaussS){
  if(Gauss){
    what="Gauss"
  }else{
    what="Student"
  }
  for(i in 1:1){
    filename=paste(dirSave,"/allres_","T",N,"_Navg",Navg,"_Nsamples",Nsamples,"_Gauss",Gauss,"_",i,"_SNRmin",SNRmin,".rda",sep="")
    print(filename)
    if(file.exists(filename) && file.info(filename)$size>4000 && !forceRUN){
        print("  already computed, skipping")
        next
    }
    allres=lapply(toList(SNRs),function(SNR){
      print(paste(what,SNR))
      return(sampleEstimates(N=N, Navg=Navg, numPerm=Nsamples, SNR=SNR, dtype=what,mylapply=mclapply))
    })
    save(allres,file=filename)
    rm(allres)
    gc()
  }
}
