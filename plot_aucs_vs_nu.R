library(xts)
library(data.table)
library(ggplot2)

source('ggplot_theme_bw_latex.R')
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

load('aucs_onesided_20150622.rda')
aucs=as.data.table(aucs)
mystats=aucs$stat
mystats=as.vector(mystats)
mystats[mystats=="R0s"]="r"
mystats[mystats=="snrs"]="t"
mystats[mystats=="SIGN"]="sign"
mystats[mystats=="wilcox"]="Wilcox"
mystats=as.factor(mystats)
aucs$stat=mystats
setnames(aucs,old="stat",new="statistics")

pdf('auc_Student_vs_nu_onesided.pdf')
qplot(data=aucs[dtype=="Student" & SNR==0.11 & statistics %in% c("r","t","Wilcox","sign"),],nu,auc,geom="line",group=statistics,color=statistics,main="Student",ylab="AUC",lwd=I(1.1))+geom_linerange(aes(ymin=auc.ci_lower,ymax=auc.ci_upper))+ theme_bw_latex+scale_colour_manual(values=cbPalette)
dev.off()

pdf('auc_Gauss_vs_SNR_onesided.pdf')
qplot(data=aucs[dtype=="Student" &nu==4 & statistics %in% c("r","t","Wilcox","sign"),],SNR,auc,geom="line",group=statistics,color=statistics,main="Gauss",ylab="AUC",lwd=I(1.1))+geom_linerange(aes(ymin=auc.ci_lower,ymax=auc.ci_upper))+ theme_bw_latex+scale_colour_manual(values=cbPalette)
dev.off()


pdf('auc_Exponential_vs_SNR_onesided.pdf')
qplot(data=aucs[dtype=="Exponential" & statistics %in% c("r","t","Wilcox","sign"),],SNR,auc,geom="line",group=statistics,color=statistics,main="Exponential",ylab="AUC",lwd=I(1.1))+geom_linerange(aes(ymin=auc.ci_lower,ymax=auc.ci_upper))+ theme_bw_latex+scale_colour_manual(values=cbPalette)
dev.off()


pdf('auc_Student_vs_SNR_nu3_5_onesided.pdf')
qplot(data=aucs[dtype=="Student" &nu==3.5 & statistics %in% c("r","t","Wilcox","sign"),],SNR,auc,geom="line",group=statistics,color=statistics,main="Student nu=3.5",ylab="AUC",lwd=I(1.1))+geom_linerange(aes(ymin=auc.ci_lower,ymax=auc.ci_upper))+ theme_bw_latex+scale_colour_manual(values=cbPalette)
dev.off()
