library(ggplot2)
library(Rcpp)
sourceCpp("num_records.cpp")

plotCumRetsCols=function(rets,cols,lwd=2,filename=NA,colMin=1,colMax=1,ylab="cumulative sum"){
  srets=c(0,cumsum(rets))
  
  x=1:length(srets)-1
  x=gl(length(x),2 ,labels=x)
  x=as.numeric(as.vector(x))
  x=head(x,-1)
  x=tail(x,-1)
  
  y=srets
  y=gl(length(y),2 ,labels=y)
  y=as.numeric(as.vector(y))
  y=head(y,-1)
  y=tail(y,-1)
  
  
  cols=gl(length(rets),2 ,labels=cols)
  cols=as.factor(as.vector(cols))
  cols=as.numeric(as.vector(cols))
  
  ym=cummax(y)
  myDF=data.frame(x,y,ym,cols)  
  
  cols_unique=unique(cols)
  
  if(!is.na(filename)){
    pdf(filename)
  }
  
  i=1
  par(family="serif",font=2)
  plot(myDF[cols==cols_unique[i],]$x,myDF[cols==cols_unique[i],]$y,t='l',col=cols[i],lwd=lwd,ylab=ylab,xlab="t",xlim=c(0,length(rets)),ylim=range(srets),cex.lab=1.7,cex.axis=1.5)
  for(i in 2:length(rets)){
    lines(myDF[cols==cols_unique[i],]$x,myDF[cols==cols_unique[i],]$y,t='l',col=cols_unique[i],lwd=5)
  }
  lines(tail(unique(x),-1),cummax(tail(srets,-1)),lty=1,lwd=lwd,t='s',col=colMax)
  lines(tail(unique(x),-1),cummin(tail(srets,-1)),lty=1,lwd=lwd,t='s',col=colMin)
  
  if(!is.na(filename)){
    dev.off() 
  }
  #  print(qplot(data=myDF,x,y,group=cols,col=as.numeric(cols),geom="line",lwd=2),ylab="price",xlab="time")+geom_line(data=myDF,aes(x,ym),lty=2)
  
}

plotRetsCols=function(rets,cols,filename=NA){
  
  x=1:length(rets)
  x=rep(x,each=2)
  
  
  y=data.frame(rep(0,length(rets)),rets)
  y=as.vector(t(y))
  cols=rep(cols,each=2)
  
  myDF=data.frame(x,y,cols)  
  

  if(!is.na(filename)){
    pdf(filename)
  }
  
  cols_unique=unique(cols)
  
  i=1
  plot(myDF[cols==cols_unique[i],]$x,myDF[cols==cols_unique[i],]$y,t='l',col=cols_unique[i],lwd=5,ylab="value",xlab="sample number",xlim=c(1,length(rets)),ylim=range(y),cex.lab=1.7,cex.axis=1.5)
  for(i in 2:length(rets)){
    lines(myDF[cols==cols_unique[i],]$x,myDF[cols==cols_unique[i],]$y,t='l',col=cols_unique[i],lwd=5)
  }
  abline(h=0,lty=2,lwd=2)
  
  if(!is.na(filename)){
    dev.off() 
  }
  #  print(qplot(data=myDF,x,y,group=cols,col=as.numeric(cols),geom="line",lwd=2),ylab="price",xlab="time")+geom_line(data=myDF,aes(x,ym),lty=2)
  
}


rets=c(2,1,-4,7,-2,-1)
cols=seq_along(rets)

myorder=seq_along(rets)

Nplots=5

for(i in 1:Nplots){
  filename=paste0("plot",i,".pdf")
  filename=NA
  plotCumRetsCols(rets[myorder],cols[myorder],filename=filename)
  print(num_records_up(cumsum(rets[myorder])))
  print(num_records_down(cumsum(rets[myorder])))
  filename=paste0("rets_",filename)
  plotRetsCols(rets[myorder],cols[myorder],filename = filename)
  myorder=sample(myorder)
}
  
  

