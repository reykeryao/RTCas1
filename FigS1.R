rm(list=ls())
setwd("/stor/work/Lambowitz/yaojun/Work/RTCas1/Documents/Script/")
library(RColorBrewer)
library(stringi)
library(grDevices)
library(showtext)
set.seed(740714)
font_add_google("Arimo","Arial")
### TA ana
### WT
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

WT_NTA<-table(substr(WT$END[nchar(WT$END)>1],1,2))
RTD_NTA<-table(substr(RTD$END[nchar(RTD$END)>1],1,2))
NTA2<-data.frame("WT"=WT_NTA[grep("N",WT_NTA,invert=T)],
                 "RTD"=RTD_NTA[grep("N",names(RTD_NTA),invert=T)])
rownames(NTA2)<-NTA2[,1]
NTA2<-NTA2[,c(-1,-3)]
NTAsum<-colSums(NTA2)
NTA2$pvalue<-apply(NTA2,1,function(x){chisq.test(rbind(x,NTAsum))$p.value})
NTA2$sign<-""
NTA2[,1:2]<-100*prop.table(as.matrix(NTA2[,1:2]),2)

NTA2$sign[NTA2$pvalue<0.001 & NTA2$WT.Freq>NTA2$RTD.Freq]<-"*"
row.names(NTA2)<-paste0(row.names(NTA2),NTA2$sign)
NTA2<-NTA2[,1:2]



pdf("FigS1B.pdf")
scol <- c(brewer.pal(8, "Set3"),brewer.pal(8, "Set1"))
mp<-barplot(as.matrix(NTA2),col=scol,axes=F,
            names.arg=NTAsum,
            width=0.15,
            space=0.3,
            xlim=c(0,.8),
            legend.text=rownames(NTA2),
            adj=0,args.legend=list(x=0.4,xjust=0,y=50,yjust=0.5,bty="n",cex=0.7))
mtext("Spacers with non-encoded nucleotides
      (5' end of the NEU)",line=2,at=median(mp),cex=0.7)
mtext(c("WT","RTD"),line=0.25,at=mp,cex=0.7)
axis(1,labels=NA,at=c(0,0.4),lwd=1,lwd.ticks=0)
axis(1,labels=NA,at=mp,lwd=0,lwd.ticks=1)
axis(2,labels=c(0,25,50,75,100),las=1,at=c(0,25,50,75,100),pos=0,cex.axis=0.7)
par(xpd=T)
text(x=-0.1,y=50,"Spacers (%)",srt=90,adj=0.5)
dev.off()

##geno cov
WT_gen<-read.delim("WT_Spacer_intersect.bed",header=F)
RTD_gen<-read.delim("RTD_Spacer_intersect.bed",header=F)
colnames(WT_gen)[4]<-colnames(RTD_gen)[4]<-"ID"
###
WT_gen<-WT_gen[,c(-1,-5,-7,-11,-15)]
colnames(WT_gen)<-c("SP_start","Sp_end","ID","Sp_strand","Host_start","Host_end",
                           "Host_ID","Host_strand","Host_name","Host_type")
WT_gen$Gene_pos<-0
for (i in 1:dim(WT_gen)[1]){
  G_strand<-WT_gen$Host_strand[i]
  if (G_strand=="+"){
    Start<-WT_gen$Host_start[i]
    End<-WT_gen$Host_end[i]
    Len<-abs(Start-End+1)
    Sp_mid<-mean(WT_gen$SP_start[i],WT_gen$Sp_end[i])
    Sp_pos<-abs(Sp_mid-Start+1)
    WT_gen$Gene_pos[i]<-100*Sp_pos/Len
  }else if (G_strand=="-"){
    Start<-WT_gen$Host_end[i]
    End<-WT_gen$Host_start[i]
    Len<-abs(Start-End+1)
    Sp_mid<-mean(WT_gen$SP_start[i],WT_gen$Sp_end[i])
    Sp_pos<-abs(Sp_mid-Start)
    WT_gen$Gene_pos[i]<-100*Sp_pos/Len
  }
}

RTD_gen<-RTD_gen[,c(-1,-5,-7,-11,-15)]
colnames(RTD_gen)<-c("SP_start","Sp_end","ID","Sp_strand","Host_start","Host_end",
                    "Host_ID","Host_strand","Host_name","Host_type")
RTD_gen$Gene_pos<-0
for (i in 1:dim(RTD_gen)[1]){
  G_strand<-RTD_gen$Host_strand[i]
  if (G_strand=="+"){
    Start<-RTD_gen$Host_start[i]
    End<-RTD_gen$Host_end[i]
    Len<-abs(Start-End+1)
    Sp_mid<-mean(RTD_gen$SP_start[i],RTD_gen$Sp_end[i])
    Sp_pos<-abs(Sp_mid-Start+1)
    RTD_gen$Gene_pos[i]<-100*Sp_pos/Len
  }else if (G_strand=="-"){
    Start<-RTD_gen$Host_end[i]
    End<-RTD_gen$Host_start[i]
    Len<-abs(Start-End+1)
    Sp_mid<-mean(RTD_gen$SP_start[i],RTD_gen$Sp_end[i])
    Sp_pos<-abs(Sp_mid-Start)
    RTD_gen$Gene_pos[i]<-100*Sp_pos/Len
  }
}

pdf("FigS1C.pdf",family = "ArialMT",width=10,height=5)
plot(density(WT_gen$Gene_pos),main="Relative position of acquired spacers",
     xlab=c("Genebody (5' to 3', %)"),ylab="Density")
lines(density(RTD_gen$Gene_pos),col="red")
sums<-c(dim(WT_gen)[1],dim(RTD_gen)[1])
legend(50,0.004,legend = paste0(c("WT (","RTD ("),sums," spacers)"),col=c("black","red"),lty=1,bty="n")
dev.off()

