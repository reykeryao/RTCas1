rm(list=ls())
setwd("/stor/work/Lambowitz/yaojun/Work/RTCas/Documents/Script/")
library(RColorBrewer)
library(stringi)
col<-brewer.pal(6,"Paired")[c(1,3,5,2,4,6)]

#function to make RC sequences
rev.comp<-function(x,rev=TRUE)
{
  x<-toupper(x)
  y<-rep("N",nchar(x))
  xx<-unlist(strsplit(x,NULL))
  for (bbb in 1:nchar(x))
  {
    if(xx[bbb]=="A") y[bbb]<-"T"        
    if(xx[bbb]=="C") y[bbb]<-"G"        
    if(xx[bbb]=="G") y[bbb]<-"C"        
    if(xx[bbb]=="T") y[bbb]<-"A"
  }
  if(rev==FALSE) 
  {
    for(ccc in (1:nchar(x)))
    {
      if(ccc==1) yy<-y[ccc] else yy<-paste(yy,y[ccc],sep="")
    }
  }
  if(rev==T)
  {
    zz<-rep(NA,nchar(x))
    for(ccc in (1:nchar(x)))
    {
      zz[ccc]<-y[nchar(x)+1-ccc]
      if(ccc==1) yy<-zz[ccc] else yy<-paste(yy,zz[ccc],sep="")
    }
  }
  return(yy)    
}


WT<-read.delim("MMB1_WT_sp_detail.info")
RTD<-read.delim("MMB1_RTdelta_sp_detail.info")
WT$tail<-0
WT$END<-""
WT$SpacerSeq<-""
for (i in 1:dim(WT)[1]){
  if (WT[i,15]!="U"){
    WT[i,17]<-WT[i,4]+WT[i,5]
  }
  if (WT[i,15]=="T"){
    WT[i,18]<-WT[i,7]
    WT[i,19]<-substr(WT[i,3],1,nchar(WT[i,3])-WT[i,5])
  } else if (WT[i,15]=="B"){
    WT[i,18]<-WT[i,8]
    WT[i,19]<-rev.comp(substr(WT[i,3],WT[i,4]+1,nchar(WT[i,3])))
  }
}

RTD$tail<-0
RTD$END<-""
RTD$SpacerSeq<-""
for (i in 1:dim(RTD)[1]){
  if (RTD[i,15]!="U"){
    RTD[i,17]<-RTD[i,4]+RTD[i,5]
  }
  if (RTD[i,15]=="T"){
    RTD[i,18]<-RTD[i,7]
    RTD[i,19]<-substr(RTD[i,3],1,nchar(RTD[i,3])-RTD[i,5])
  } else if (RTD[i,15]=="B"){
    RTD[i,18]<-RTD[i,9]
    RTD[i,19]<-rev.comp(substr(RTD[i,3],RTD[i,4]+1,nchar(RTD[i,3])))
  }
}


## Fig1D
pdf("Fig1D.pdf",width=8,height=3)
par(mfrow=c(1,2),cex=0.7)
#from the start of the soft-clipped seq (+1 to +5 from the end of the sapcer)
plotspcar_stack<-function(data_df,title){
  tmp<-data.frame("Nucl"=c("A","G","C","T"))
  for (i in 1:5){
    tmp<-merge(tmp,data.frame(table(substr(data_df$END,i,i))),by=1,all.x=T)
  }
  y_max<-round(max(colSums(tmp[,2:6]))/500)*500
  colnames(tmp)[2:6]<-c("1","2","3","4","5")
  df.bar <-barplot(as.matrix(tmp[,2:6]),ylab="Spacers", main=title,width = 0.1, ylim=c(0,y_max),
                   beside = F,col=c("green","blue","chocolate","red"),
                   args.legend = list(x="right",bty="n"),legend.text = tmp$Nucl,
                   space=0.1,xlim=c(0,0.7),xlab="Position from the end of spacer")
}
df1<-WT[WT$Sp_type!="U",]
plotspcar_stack(df1,"WT")

df1<-RTD[RTD$Sp_type!="U",]
plotspcar_stack(df1,"RTdelta")
dev.off()

Full<-data.frame(table(nchar(WT$Seq)))
Full<-merge(Full,data.frame(table(nchar(RTD$Seq))),by=1,all=T)
colnames(Full)<-c("Len","Mm-WT","Mm-RTD")
Full[is.na(Full)]<-0
Full$Len<-as.integer(as.character(Full$Len))
Len<-data.frame("Len"=20:50)
Full<-merge(Len,Full,by=1,all=T)
Full[is.na(Full)]<-0
Full[,2:3]<-100*prop.table(as.matrix(Full[,2:3]),margin =2)
WoSS<-data.frame(table(nchar(WT$Seq)-WT$SS5-WT$SS3))
WoSS<-merge(WoSS,data.frame(table(nchar(RTD$Seq)-RTD$SS5-RTD$SS3)),by=1,all=T)
colnames(WoSS)<-c("Len","Mm-WT","Mm-RTD")
WoSS[is.na(WoSS)]<-0
WoSS$Len<-as.integer(as.character(WoSS$Len))
Len<-data.frame("Len"=20:50)
WoSS<-merge(Len,WoSS,by=1,all=T)
WoSS[is.na(WoSS)]<-0
WoSS[,2:3]<-100*prop.table(as.matrix(WoSS[,2:3]),margin =2)

pdf("Fig1E.pdf",height=3,width=8)
par(mfrow=c(1,2),bty="n",cex=0.7)
Full<-Full[Full$Len<=50,]
plot(Full$`Mm-WT`~Full$Len,type="l",col="black",xlim=c(20,50),yaxt="n",
     bty="n",ylim=c(0,50),
     xlab="Length (nt)",ylab="Frequency (%)")
title(main="Spacer length\n(full sequence)")
lines(Full$`Mm-RTD`~Full$Len,col="red")
#lines(Full$Fs~Full$Len,col="green")
#lines(Full$Vv~Full$Len,col="goldenrod1")
#lines(Full$TT~Full$Len,col="purple")
points(Full$`Mm-WT`~Full$Len,pch=0)
points(Full$`Mm-RTD`~Full$Len,pch=1,col="red")
axis(2,at=seq(0,50,10),labels = seq(0,50,10),las=2)
#points(Full$Fs~Full$Len,pch=2,col="green")
#points(Full$Vv~Full$Len,pch=3,col="goldenrod1")
#points(Full$TT~Full$Len,pch=4,col="purple")
legend("topright",legend = c("Mm: WT",expression(paste("Mm: RT",Delta))),lty=1,col=c("black","red"),
       bty="n",pch=c(0:1))
WoSS<-WoSS[WoSS$Len<=50,]
plot(WoSS$`Mm-WT`~WoSS$Len,type="l",col="black",xlim=c(20,50),ylim=c(0,50),
     xlab="Length (nt)",ylab="Frequency (%)",bty="n",yaxt="n")
title(main="Spacer length\n(minus non-encoded nucleotides)")
lines(WoSS$`Mm-RTD`~WoSS$Len,col="red")
points(WoSS$`Mm-WT`~WoSS$Len,pch=0)
points(WoSS$`Mm-RTD`~WoSS$Len,pch=1,col="red")
axis(2,at=seq(0,50,10),labels = seq(0,50,10),las=2)
dev.off()
