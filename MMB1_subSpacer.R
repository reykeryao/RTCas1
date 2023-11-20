rm(list=ls())
setwd("/stor/work/Lambowitz/yaojun/Work/RTCas1/Documents/Script/")
library(RColorBrewer)
library(stringi)
library(stringr)
library(showtext)
set.seed(740714)
font_add_google("Arimo","Arial")
library(Polychrome)
P64 = createPalette(64,  c("#ff0000", "#00ff00", "#0000ff"))
#swatch(P48)
P16 = createPalette(16,  c("#ff0000", "#00ff00", "#0000ff"))
#swatch(P22)

setwd("/stor/work/Lambowitz/yaojun/Work/RTCas1/Documents/Script/")

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

### isolate spacers mapped to know gene (can define sense/antisense) 
### and with SS in only one end
WT<-read.delim("MMB1_WT_sp_detail.info")
RTD<-read.delim("MMB1_RTdelta_sp_detail.info")

WT<-WT[WT$Temp_type!="U",]
WT<-WT[WT$SS5==0 | WT$SS3==0,]
RTD<-RTD[RTD$Temp_type!="U",]
RTD<-RTD[RTD$SS5==0 | RTD$SS3==0,]

### rearrange spacers into sense only, 2 type2 sense-SS, or SS-sense
WT$SSeq5<-""
WT$SpacerSeq<-""
WT$SSeq3<-""
for (i in 1:dim(WT)[1]){
  if (WT$Temp_type[i]=="S"){
    if (WT$SS5[i]>0){
      WT$SSeq5[i]<-WT$SS5_seq[i]
      WT$SpacerSeq[i]<-substr(WT$Seq[i],WT$SS5[i]+1,nchar(WT$Seq[i]))
    } else if (WT$SS3[i]>0) {
      WT$SpacerSeq[i]<-substr(WT$Seq[i],1,nchar(WT$Seq[i])-WT$SS3[i])
      WT$SSeq3[i]<-WT$SS3_seq[i]
    } else {
      WT$SpacerSeq[i]<-WT$Seq[i]
    }
  } else if (WT$Temp_type[i]=="A"){
    if (WT$SS5[i]>0){
      WT$SpacerSeq[i]<-rev.comp(substr(WT$Seq[i],WT$SS5[i]+1,nchar(WT$Seq[i])))
      WT$SSeq3[i]<-WT$SS5_rev[i]
    } else if (WT$SS3[i]>0) {
      WT$SSeq5[i]<-WT$SS3_rev[i]
      WT$SpacerSeq[i]<-rev.comp(substr(WT$Seq[i],1,nchar(WT$Seq[i])-WT$SS3[i]))
    } else {
      WT$SpacerSeq[i]<-rev.comp(WT$Seq[i])
    }
  }
}

RTD$SSeq5<-""
RTD$SpacerSeq<-""
RTD$SSeq3<-""
for (i in 1:dim(RTD)[1]){
  if (RTD$Temp_type[i]=="S"){
    if (RTD$SS5[i]>0){
      RTD$SSeq5[i]<-RTD$SS5_seq[i]
      RTD$SpacerSeq[i]<-substr(RTD$Seq[i],RTD$SS5[i]+1,nchar(RTD$Seq[i]))
    } else if (RTD$SS3[i]>0) {
      RTD$SpacerSeq[i]<-substr(RTD$Seq[i],1,nchar(RTD$Seq[i])-RTD$SS3[i])
      RTD$SSeq3[i]<-RTD$SS3_seq[i]
    } else {
      RTD$SpacerSeq[i]<-RTD$Seq[i]
    }
  } else if (RTD$Temp_type[i]=="A"){
    if (RTD$SS5[i]>0){
      RTD$SpacerSeq[i]<-rev.comp(substr(RTD$Seq[i],RTD$SS5[i]+1,nchar(RTD$Seq[i])))
      RTD$SSeq3[i]<-RTD$SS5_rev[i]
    } else if (RTD$SS3[i]>0) {
      RTD$SSeq5[i]<-RTD$SS3_rev[i]
      RTD$SpacerSeq[i]<-rev.comp(substr(RTD$Seq[i],1,nchar(RTD$Seq[i])-RTD$SS3[i]))
    } else {
      RTD$SpacerSeq[i]<-rev.comp(RTD$Seq[i])
    }
  }
}

### dinucleotides frq ana (SS>=4)
### SS5--sp--SS3
WT$SS5_1<-substr(WT$SSeq5,1,2)
WT$SS5_2<-substr(lapply(WT$SSeq5,FUN = function(x){intToUtf8(rev(utf8ToInt(x)))}),1,2)
WT$SS5_2<-lapply(WT$SS5_2,FUN = function(x){intToUtf8(rev(utf8ToInt(x)))})
WT$SP_1<-substr(WT$SpacerSeq,1,2)
WT$SP_2<-substr(lapply(WT$SpacerSeq,FUN = function(x){intToUtf8(rev(utf8ToInt(x)))}),1,2)
WT$SP_2<-lapply(WT$SP_2,FUN = function(x){intToUtf8(rev(utf8ToInt(x)))})
WT$SS3_1<-substr(WT$SSeq3,1,2)
WT$SS3_2<-substr(lapply(WT$SSeq3,FUN = function(x){intToUtf8(rev(utf8ToInt(x)))}),1,2)
WT$SS3_2<-lapply(WT$SS3_2,FUN = function(x){intToUtf8(rev(utf8ToInt(x)))})
  
RTD$SS5_1<-substr(RTD$SSeq5,1,2)
RTD$SS5_2<-substr(lapply(RTD$SSeq5,FUN = function(x){intToUtf8(rev(utf8ToInt(x)))}),1,2)
RTD$SS5_2<-lapply(RTD$SS5_2,FUN = function(x){intToUtf8(rev(utf8ToInt(x)))})
RTD$SP_1<-substr(RTD$SpacerSeq,1,2)
RTD$SP_2<-substr(lapply(RTD$SpacerSeq,FUN = function(x){intToUtf8(rev(utf8ToInt(x)))}),1,2)
RTD$SP_2<-lapply(RTD$SP_2,FUN = function(x){intToUtf8(rev(utf8ToInt(x)))})
RTD$SS3_1<-substr(RTD$SSeq3,1,2)
RTD$SS3_2<-substr(lapply(RTD$SSeq3,FUN = function(x){intToUtf8(rev(utf8ToInt(x)))}),1,2)
RTD$SS3_2<-lapply(RTD$SS3_2,FUN = function(x){intToUtf8(rev(utf8ToInt(x)))})

make_fq_df<-function(x,y){
  x<-data.frame(table(x[x!="" & nchar(x)>1]))
  y<-data.frame(table(y[y!="" & nchar(y)>1]))
  df<-merge(x[grep("N",x[,1],invert=T),],y[grep("N",y[,1],invert=T),],by=1,all=T)
  df[is.na(df)]<-0
  colnames(df)<-c("NT","WT","RTD")
  return(df)
}


pdf("MMSpacer_NTA8.pdf")
par(mfrow=c(3,2))
scol <- c(brewer.pal(8, "Set3"),brewer.pal(8, "Set1"))
df<-make_fq_df(WT$SS5_1[nchar(WT$SSeq5)>=4],RTD$SS5_1[nchar(RTD$SSeq5)>=4])
mp<-barplot(100*prop.table(as.matrix(df[,-1]),2),col=scol,axes=F,
            names.arg=colSums(df[,-1]),
            width=0.15,
            space=0.3,
            xlim=c(0,.8),
            legend.text=df$NT,
            adj=0,args.legend=list(x=0.4,xjust=0,y=50,yjust=0.5,bty="n",cex=0.7))
mtext("SS5-1",line=2,at=median(mp),cex=0.7)
mtext(c("WT","RTD"),line=0.25,at=mp,cex=0.7)
axis(1,labels=NA,at=c(0,0.4),lwd=1,lwd.ticks=0)
axis(1,labels=NA,at=mp,lwd=0,lwd.ticks=1)
axis(2,labels=c(0,25,50,75,100),las=1,at=c(0,25,50,75,100),pos=0,cex.axis=0.7)
par(xpd=T)
text(x=-0.15,y=50,"Spacers (%)",srt=90,adj=0.5,cex=0.7)


mp<-barplot(Last2,col=scol,axes=F,
            names.arg=c("WT","RTD"),
            width=0.15,
            space=0.3,
            xlim=c(0,.8),
            legend.text=rownames(Last2),
            adj=0,args.legend=list(x=0.4,xjust=0,y=50,yjust=0.5,bty="n",cex=0.7))
mtext("Spacers with non-encoded nucleotides
      (3' end of the spacer)",line=2,at=median(mp),cex=0.7)
mtext(Last2sum,line=0.25,at=mp,cex=0.7)
axis(1,labels=NA,at=c(0,0.4),lwd=1,lwd.ticks=0)
axis(1,labels=NA,at=mp,lwd=0,lwd.ticks=1)
axis(2,labels=c(0,25,50,75,100),las=1,at=c(0,25,50,75,100),pos=0,cex.axis=0.7)
par(xpd=T)
text(x=-0.15,y=50,"Spacers (%)",srt=90,adj=0.5,cex=0.7)

mp<-barplot(NTA2,col=scol,axes=F,
            names.arg=c("WT","RTD"),
            width=0.15,
            space=0.3,
            xlim=c(0,.8),
            legend.text=rownames(NTA2),
            adj=0,args.legend=list(x=0.4,xjust=0,y=50,yjust=0.5,bty="n",cex=0.7))
mtext("Spacers with non-encoded nucleotides
      (5' end of the NEU)",line=2,at=median(mp),cex=0.7)
mtext(NTAsum,line=0.25,at=mp,cex=0.7)
axis(1,labels=NA,at=c(0,0.4),lwd=1,lwd.ticks=0)
axis(1,labels=NA,at=mp,lwd=0,lwd.ticks=1)
axis(2,labels=c(0,25,50,75,100),las=1,at=c(0,25,50,75,100),pos=0,cex.axis=0.7)
par(xpd=T)
text(x=-0.15,y=50,"Spacers (%)",srt=90,adj=0.5,cex=0.7)

mp<-barplot(NTA2_3end,col=scol,axes=F,
            names.arg=c("WT","RTD"),
            width=0.15,
            space=0.3,
            xlim=c(0,.8),
            legend.text=rownames(NTA2_3end),
            adj=0,args.legend=list(x=0.4,xjust=0,y=50,yjust=0.5,bty="n",cex=0.7))
mtext("Spacers with non-encoded nucleotides
      (3' end of the NEU)",line=2,at=median(mp),cex=0.7)
mtext(NTA2_3endsum,line=0.25,at=mp,cex=0.7)
axis(1,labels=NA,at=c(0,0.4),lwd=1,lwd.ticks=0)
axis(1,labels=NA,at=mp,lwd=0,lwd.ticks=1)
axis(2,labels=c(0,25,50,75,100),las=1,at=c(0,25,50,75,100),pos=0,cex.axis=0.7)
par(xpd=T)
text(x=-0.15,y=50,"Spacers (%)",srt=90,adj=0.5,cex=0.7)
u <- par("usr")
mtext("*: p-value<0.001 (chi-squre test)", side=3, line=-20, outer=TRUE,adj=0.5,cex=0.7)
dev.off()

pdf("MMSpacer_NTA2.pdf")
par(mfrow=c(3,2))
scol <- c(brewer.pal(8, "Set3"),brewer.pal(8, "Set1"))
mp<-barplot(NTA2,col=scol,axes=F,
            names.arg=c("WT","RTD"),
            width=0.15,
            space=0.3,
            xlim=c(0,.8),
            legend.text=rownames(NTA2),
            adj=0,args.legend=list(x=0.4,xjust=0,y=50,yjust=0.5,bty="n",cex=0.7))
mtext("Spacers with non-encoded nucleotides
      (5' end of the NEU)",line=2,at=median(mp),cex=0.7)
mtext(NTAsum,line=0.25,at=mp,cex=0.7)
axis(1,labels=NA,at=c(0,0.4),lwd=1,lwd.ticks=0)
axis(1,labels=NA,at=mp,lwd=0,lwd.ticks=1)
axis(2,labels=c(0,25,50,75,100),las=1,at=c(0,25,50,75,100),pos=0,cex.axis=0.7)
par(xpd=T)
text(x=-0.1,y=50,"Spacers (%)",srt=90,adj=0.5)

mp<-barplot(First2,col=scol,axes=F,
            names.arg=c("WT","RTD"),
            width=0.15,
            space=0.3,
            xlim=c(0,.8),
            legend.text=rownames(First2),
            adj=0,args.legend=list(x=0.4,xjust=0,y=50,yjust=0.5,bty="n",cex=0.7))
mtext("Spacers with non-encoded nucleotides
      (5' end of the spacer)",line=2,at=median(mp),cex=0.7)
mtext(First2sum,line=0.25,at=mp,cex=0.7)
axis(1,labels=NA,at=c(0,0.4),lwd=1,lwd.ticks=0)
axis(1,labels=NA,at=mp,lwd=0,lwd.ticks=1)
axis(2,labels=c(0,25,50,75,100),las=1,at=c(0,25,50,75,100),pos=0,cex.axis=0.7)
par(xpd=T)
text(x=-0.1,y=50,"Spacers (%)",srt=90,adj=0.5)

df1<-WT[WT$Sp_type!="U",]
plotsp5NT(df1,"WT")
plotsp3NT(df1,"WT")

df1<-RTD[RTD$Sp_type!="U",]
plotsp5NT(df1,"RTD")
plotsp3NT(df1,"RTD")

dev.off()

### GC density
proteinGC<-read.table("../../MmRTCas/Ref_MMB1/protein.fa.tab",col.names=c("ID","Seq"))
proteinGC$GC<-(str_count(proteinGC$Seq,"G")+str_count(proteinGC$Seq,"C"))/
  nchar(proteinGC$Seq) * 100 
### arrange 1
### non-encoded nucleotides are always at 3' end, spacers without NEU is discarded
### sapcers without NEU are all flipped into sense strand
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
write.table(WT,"WT.rearr1",sep="\t",quote=F,row.names=F)
write.table(RTD,"RTD.rearr1",sep="\t",quote=F,row.names=F)


### arrange 2
### encoded spacers are always in sense orientation
### isolate spacers mapped to know gene (can define sense/antisense) 
### and with SS in at only one end

WT<-read.delim("MMB1_WT_sp_detail.info")
RTD<-read.delim("MMB1_RTdelta_sp_detail.info")

WT<-WT[WT$Temp_type!="U",]
WT<-WT[WT$SS5==0 | WT$SS3==0,]
RTD<-RTD[RTD$Temp_type!="U",]
RTD<-RTD[RTD$SS5==0 | RTD$SS3==0,]

### rearrange spacers into sense only, 2 type2 sense-SS, or SS-sense
WT$SSeq5<-""
WT$SpacerSeq<-""
WT$SSeq3<-""
for (i in 1:dim(WT)[1]){
  if (WT$Temp_type[i]=="S"){
    if (WT$SS5[i]>0){
      WT$SSeq5[i]<-WT$SS5_seq[i]
      WT$SpacerSeq[i]<-substr(WT$Seq[i],WT$SS5[i]+1,nchar(WT$Seq[i]))
    } else if (WT$SS3[i]>0) {
      WT$SpacerSeq[i]<-substr(WT$Seq[i],1,nchar(WT$Seq[i])-WT$SS3[i])
      WT$SSeq3[i]<-WT$SS3_seq[i]
    } else {
      WT$SpacerSeq[i]<-WT$Seq[i]
    }
  } else if (WT$Temp_type[i]=="A"){
    if (WT$SS5[i]>0){
      WT$SpacerSeq[i]<-rev.comp(substr(WT$Seq[i],WT$SS5[i]+1,nchar(WT$Seq[i])))
      WT$SSeq3[i]<-WT$SS5_rev[i]
    } else if (WT$SS3[i]>0) {
      WT$SSeq5[i]<-WT$SS3_rev[i]
      WT$SpacerSeq[i]<-rev.comp(substr(WT$Seq[i],1,nchar(WT$Seq[i])-WT$SS3[i]))
    } else {
      WT$SpacerSeq[i]<-rev.comp(WT$Seq[i])
    }
  }
}

RTD$SSeq5<-""
RTD$SpacerSeq<-""
RTD$SSeq3<-""
for (i in 1:dim(RTD)[1]){
  if (RTD$Temp_type[i]=="S"){
    if (RTD$SS5[i]>0){
      RTD$SSeq5[i]<-RTD$SS5_seq[i]
      RTD$SpacerSeq[i]<-substr(RTD$Seq[i],RTD$SS5[i]+1,nchar(RTD$Seq[i]))
    } else if (RTD$SS3[i]>0) {
      RTD$SpacerSeq[i]<-substr(RTD$Seq[i],1,nchar(RTD$Seq[i])-RTD$SS3[i])
      RTD$SSeq3[i]<-RTD$SS3_seq[i]
    } else {
      RTD$SpacerSeq[i]<-RTD$Seq[i]
    }
  } else if (RTD$Temp_type[i]=="A"){
    if (RTD$SS5[i]>0){
      RTD$SpacerSeq[i]<-rev.comp(substr(RTD$Seq[i],RTD$SS5[i]+1,nchar(RTD$Seq[i])))
      RTD$SSeq3[i]<-RTD$SS5_rev[i]
    } else if (RTD$SS3[i]>0) {
      RTD$SSeq5[i]<-RTD$SS3_rev[i]
      RTD$SpacerSeq[i]<-rev.comp(substr(RTD$Seq[i],1,nchar(RTD$Seq[i])-RTD$SS3[i]))
    } else {
      RTD$SpacerSeq[i]<-rev.comp(RTD$Seq[i])
    }
  }
}
write.table(WT,"WT.rearr2",sep="\t",quote=F,row.names=F)
write.table(RTD,"RTD.rearr2",sep="\t",quote=F,row.names=F)


pdf("spacer_GC.pdf")
par(mfrow=c(1,2))
WT<-read.delim("WT.rearr1")
RTD<-read.delim("RTD.rearr1")
WT<-WT[WT$Sp_type!="U",]
RTD<-RTD[RTD$Sp_type!="U",]
WT$GC<-(str_count(WT$SpacerSeq,"G")+str_count(WT$SpacerSeq,"C"))/
  nchar(WT$SpacerSeq) * 100 
RTD$GC<-(str_count(RTD$SpacerSeq,"G")+str_count(RTD$SpacerSeq,"C"))/
  nchar(RTD$SpacerSeq) * 100 

plot(density(WT$GC),bty="n",xlim=c(0,100),main="GC content (rearrange 1)",xlab="%",ylim=c(0,0.2))
lines(density(RTD$GC),col="red")
lines(density(proteinGC$GC),col="gray80")
legend("topleft",lty=1,col=c("gray80","black","red"),bty="n",cex=0.7,
       legend = c("Protein coding genes",
                  "WT spacers (with non-encoded nucleotides",
                  "RTD spacers (with non-encoded nucleotides"))

##raarr2
WT<-read.delim("WT.rearr2")
RTD<-read.delim("RTD.rearr2")
WT$GC<-(str_count(WT$SpacerSeq,"G")+str_count(WT$SpacerSeq,"C"))/
  nchar(WT$SpacerSeq) * 100 
RTD$GC<-(str_count(RTD$SpacerSeq,"G")+str_count(RTD$SpacerSeq,"C"))/
  nchar(RTD$SpacerSeq) * 100 

plot(density(WT$GC),bty="n",xlim=c(0,100),main="GC content (rearrange 2)",xlab="%",ylim=c(0,0.2))
lines(density(RTD$GC),col="red")
lines(density(proteinGC$GC),col="gray80")
legend("topleft",lty=1,col=c("gray80","black","red"),bty="n",cex=0.7,
       legend = c("Protein coding genes",
                  "WT spacers",
                  "RTD spacers"))
dev.off()

### tri and di frq
### read sim sp nu
WT_sim5_2<-read.table("sim_sp/sim_fa/WT_end5_2")
WT_sim5_3<-read.table("sim_sp/sim_fa/WT_end5_3")
WT_sim3_2<-read.table("sim_sp/sim_fa/WT_end3_2")
WT_sim3_3<-read.table("sim_sp/sim_fa/WT_end3_3")
RTD_sim5_2<-read.table("sim_sp/sim_fa/RTD_end5_2")
RTD_sim5_3<-read.table("sim_sp/sim_fa/RTD_end5_3")
RTD_sim3_2<-read.table("sim_sp/sim_fa/RTD_end3_2")
RTD_sim3_3<-read.table("sim_sp/sim_fa/RTD_end3_3")
##arr1
WT<-read.delim("WT.rearr1")
RTD<-read.delim("RTD.rearr1")
WT<-WT[WT$SpacerSeq!="",]
RTD<-RTD[RTD$SpacerSeq!="",]
### 
WT$end5_2<-substr(WT$SpacerSeq,1,2)
WT$end5_3<-substr(WT$SpacerSeq,1,3)
WT$end3_2<-substr(WT$SpacerSeq,nchar(WT$SpacerSeq)-1,nchar(WT$SpacerSeq))
WT$end3_3<-substr(WT$SpacerSeq,nchar(WT$SpacerSeq)-2,nchar(WT$SpacerSeq))

RTD$end5_2<-substr(RTD$SpacerSeq,1,2)
RTD$end5_3<-substr(RTD$SpacerSeq,1,3)
RTD$end3_2<-substr(RTD$SpacerSeq,nchar(RTD$SpacerSeq)-1,nchar(RTD$SpacerSeq))
RTD$end3_3<-substr(RTD$SpacerSeq,nchar(RTD$SpacerSeq)-2,nchar(RTD$SpacerSeq))

pdf("spacer-encoded1.pdf",width=8,height=12)
par(mfrow=c(2,2))
df<-make_fq_df(WT$end5_2[nchar(WT$END)>2],RTD$end5_2[nchar(RTD$END)>2])
df<-merge(df,WT_sim5_2,by=1)
df<-merge(df,RTD_sim5_2,by=1)
df<-df[,c(1,2,4,3,5)]
colnames(df)[c(3,5)]<-c("WT_sim","RTD_sim")
mp<-barplot(100*prop.table(as.matrix(df[,-1]),2),col=P16,axes=F,
            names.arg=colSums(df[,-1]),
            width=0.15,
            space=0.3,
            xlim=c(0,.8),
            legend.text=df$NT,
            adj=0,args.legend=list(x=0.8,xjust=0,y=50,yjust=0.5,bty="n",cex=0.7))
mtext("Spacers(encoded nucleotides, 5'end, NEU>2nt)",line=2,at=median(mp),cex=0.7)
mtext(c("WT","WT-sim","RTD-sim","RTD"),line=0.25,at=mp,cex=0.7)
axis(1,labels=NA,at=c(0,0.8),lwd=1,lwd.ticks=0)
axis(1,labels=NA,at=mp,lwd=0,lwd.ticks=1)
axis(2,labels=c(0,25,50,75,100),las=1,at=c(0,25,50,75,100),pos=0,cex.axis=0.7)
par(xpd=T)
text(x=-0.15,y=50,"Spacers (%)",srt=90,adj=0.5,cex=0.7)

df<-make_fq_df(WT$end5_2[nchar(WT$END)<3],RTD$end5_2[nchar(RTD$END)<3])
df<-merge(df,WT_sim5_2,by=1)
df<-merge(df,RTD_sim5_2,by=1)
df<-df[,c(1,2,4,3,5)]
colnames(df)[c(3,5)]<-c("WT_sim","RTD_sim")
mp<-barplot(100*prop.table(as.matrix(df[,-1]),2),col=P16,axes=F,
            names.arg=colSums(df[,-1]),
            width=0.15,
            space=0.3,
            xlim=c(0,.8),
            legend.text=df$NT,
            adj=0,args.legend=list(x=0.8,xjust=0,y=50,yjust=0.5,bty="n",cex=0.7))
mtext("Spacers(encoded nucleotides, 5'end, NEU 1-2nt)",line=2,at=median(mp),cex=0.7)
mtext(c("WT","WT-sim","RTD-sim","RTD"),line=0.25,at=mp,cex=0.7)
axis(1,labels=NA,at=c(0,0.8),lwd=1,lwd.ticks=0)
axis(1,labels=NA,at=mp,lwd=0,lwd.ticks=1)
axis(2,labels=c(0,25,50,75,100),las=1,at=c(0,25,50,75,100),pos=0,cex.axis=0.7)
par(xpd=T)
text(x=-0.15,y=50,"Spacers (%)",srt=90,adj=0.5,cex=0.7)

df<-make_fq_df(WT$end5_3[nchar(WT$END)>=3],RTD$end5_3[nchar(RTD$END)>=3])
df<-merge(df,WT_sim5_3,by=1)
df<-merge(df,RTD_sim5_3,by=1)
df<-df[,c(1,2,4,3,5)]
colnames(df)[c(3,5)]<-c("WT_sim","RTD_sim")
df<-df[df$NT%in%c("AAA","CCC","GGG","TTT"),]
mp<-barplot(100*prop.table(as.matrix(df[,-1]),2),col=P64,axes=F,
            names.arg=colSums(df[,-1]),
            width=0.15,
            space=0.3,
            xlim=c(0,.8),
            legend.text=df$NT,
            adj=0,args.legend=list(x=0.8,xjust=0,y=50,yjust=0.5,bty="n",cex=0.7))
mtext("Spacers(encoded nucleotides, 5'end, NEU>2nt)",line=2,at=median(mp),cex=0.7)
mtext(c("WT","WT-sim","RTD-sim","RTD"),line=0.25,at=mp,cex=0.7)
axis(1,labels=NA,at=c(0,0.8),lwd=1,lwd.ticks=0)
axis(1,labels=NA,at=mp,lwd=0,lwd.ticks=1)
axis(2,labels=c(0,25,50,75,100),las=1,at=c(0,25,50,75,100),pos=0,cex.axis=0.7)
par(xpd=T)
text(x=-0.15,y=50,"Spacers (%)",srt=90,adj=0.5,cex=0.7)

df<-make_fq_df(WT$end5_3[nchar(WT$END)<3],RTD$end5_3[nchar(RTD$END)<3])
df<-merge(df,WT_sim5_3,by=1)
df<-merge(df,RTD_sim5_3,by=1)
df<-df[,c(1,2,4,3,5)]
colnames(df)[c(3,5)]<-c("WT_sim","RTD_sim")
df<-df[df$NT%in%c("AAA","CCC","GGG","TTT"),]
mp<-barplot(100*prop.table(as.matrix(df[,-1]),2),col=P64,axes=F,
            names.arg=colSums(df[,-1]),
            width=0.15,
            space=0.3,
            xlim=c(0,.8),
            legend.text=df$NT,
            adj=0,args.legend=list(x=0.8,xjust=0,y=50,yjust=0.5,bty="n",cex=0.7))
mtext("Spacers(encoded nucleotides, 5'end, NEU 1-2nt)",line=2,at=median(mp),cex=0.7)
mtext(c("WT","WT-sim","RTD-sim","RTD"),line=0.25,at=mp,cex=0.7)
axis(1,labels=NA,at=c(0,0.8),lwd=1,lwd.ticks=0)
axis(1,labels=NA,at=mp,lwd=0,lwd.ticks=1)
axis(2,labels=c(0,25,50,75,100),las=1,at=c(0,25,50,75,100),pos=0,cex.axis=0.7)
par(xpd=T)
text(x=-0.15,y=50,"Spacers (%)",srt=90,adj=0.5,cex=0.7)
df<-make_fq_df(WT$end3_2[nchar(WT$END)>2],RTD$end3_2[nchar(RTD$END)>2])
df<-merge(df,WT_sim3_2,by=1)
df<-merge(df,RTD_sim3_2,by=1)
df<-df[,c(1,2,4,3,5)]
colnames(df)[c(3,5)]<-c("WT_sim","RTD_sim")
mp<-barplot(100*prop.table(as.matrix(df[,-1]),2),col=P16,axes=F,
            names.arg=colSums(df[,-1]),
            width=0.15,
            space=0.3,
            xlim=c(0,.8),
            legend.text=df$NT,
            adj=0,args.legend=list(x=0.8,xjust=0,y=50,yjust=0.5,bty="n",cex=0.7))
mtext("Spacers(encoded nucleotides, 3'end, NEU>2nt)",line=2,at=median(mp),cex=0.7)
mtext(c("WT","WT-sim","RTD-sim","RTD"),line=0.25,at=mp,cex=0.7)
axis(1,labels=NA,at=c(0,0.8),lwd=1,lwd.ticks=0)
axis(1,labels=NA,at=mp,lwd=0,lwd.ticks=1)
axis(2,labels=c(0,25,50,75,100),las=1,at=c(0,25,50,75,100),pos=0,cex.axis=0.7)
par(xpd=T)
text(x=-0.15,y=50,"Spacers (%)",srt=90,adj=0.5,cex=0.7)

df<-make_fq_df(WT$end3_2[nchar(WT$END)<3],RTD$end3_2[nchar(RTD$END)<3])
df<-merge(df,WT_sim3_2,by=1)
df<-merge(df,RTD_sim3_2,by=1)
df<-df[,c(1,2,4,3,5)]
colnames(df)[c(3,5)]<-c("WT_sim","RTD_sim")
mp<-barplot(100*prop.table(as.matrix(df[,-1]),2),col=P16,axes=F,
            names.arg=colSums(df[,-1]),
            width=0.15,
            space=0.3,
            xlim=c(0,.8),
            legend.text=df$NT,
            adj=0,args.legend=list(x=0.8,xjust=0,y=50,yjust=0.5,bty="n",cex=0.7))
mtext("Spacers(encoded nucleotides, 3'end, NEU 1-2nt)",line=2,at=median(mp),cex=0.7)
mtext(c("WT","WT-sim","RTD-sim","RTD"),line=0.25,at=mp,cex=0.7)
axis(1,labels=NA,at=c(0,0.8),lwd=1,lwd.ticks=0)
axis(1,labels=NA,at=mp,lwd=0,lwd.ticks=1)
axis(2,labels=c(0,25,50,75,100),las=1,at=c(0,25,50,75,100),pos=0,cex.axis=0.7)
par(xpd=T)
text(x=-0.15,y=50,"Spacers (%)",srt=90,adj=0.5,cex=0.7)

df<-make_fq_df(WT$end3_3[nchar(WT$END)>=3],RTD$end3_3[nchar(RTD$END)>=3])
df<-merge(df,WT_sim3_3,by=1)
df<-merge(df,RTD_sim3_3,by=1)
df<-df[,c(1,2,4,3,5)]
colnames(df)[c(3,5)]<-c("WT_sim","RTD_sim")
df<-df[df$NT%in%c("AAA","CCC","GGG","TTT"),]
mp<-barplot(100*prop.table(as.matrix(df[,-1]),2),col=P64,axes=F,
            names.arg=colSums(df[,-1]),
            width=0.15,
            space=0.3,
            xlim=c(0,.8),
            legend.text=df$NT,
            adj=0,args.legend=list(x=0.8,xjust=0,y=50,yjust=0.5,bty="n",cex=0.7))
mtext("Spacers(encoded nucleotides, 3'end, NEU>2nt)",line=2,at=median(mp),cex=0.7)
mtext(c("WT","WT-sim","RTD-sim","RTD"),line=0.25,at=mp,cex=0.7)
axis(1,labels=NA,at=c(0,0.8),lwd=1,lwd.ticks=0)
axis(1,labels=NA,at=mp,lwd=0,lwd.ticks=1)
axis(2,labels=c(0,25,50,75,100),las=1,at=c(0,25,50,75,100),pos=0,cex.axis=0.7)
par(xpd=T)
text(x=-0.15,y=50,"Spacers (%)",srt=90,adj=0.5,cex=0.7)

df<-make_fq_df(WT$end3_3[nchar(WT$END)<3],RTD$end3_3[nchar(RTD$END)<3])
df<-merge(df,WT_sim3_3,by=1)
df<-merge(df,RTD_sim3_3,by=1)
df<-df[,c(1,2,4,3,5)]
colnames(df)[c(3,5)]<-c("WT_sim","RTD_sim")
df<-df[df$NT%in%c("AAA","CCC","GGG","TTT"),]
mp<-barplot(100*prop.table(as.matrix(df[,-1]),2),col=P64,axes=F,
            names.arg=colSums(df[,-1]),
            width=0.15,
            space=0.3,
            xlim=c(0,.8),
            legend.text=df$NT,
            adj=0,args.legend=list(x=0.8,xjust=0,y=50,yjust=0.5,bty="n",cex=0.7))
mtext("Spacers(encoded nucleotides, 3'end, NEU 1-2nt)",line=2,at=median(mp),cex=0.7)
mtext(c("WT","WT-sim","RTD-sim","RTD"),line=0.25,at=mp,cex=0.7)
axis(1,labels=NA,at=c(0,0.8),lwd=1,lwd.ticks=0)
axis(1,labels=NA,at=mp,lwd=0,lwd.ticks=1)
axis(2,labels=c(0,25,50,75,100),las=1,at=c(0,25,50,75,100),pos=0,cex.axis=0.7)
par(xpd=T)
text(x=-0.15,y=50,"Spacers (%)",srt=90,adj=0.5,cex=0.7)
dev.off()





###########
##arr2
WT<-read.delim("WT.rearr2")
RTD<-read.delim("RTD.rearr2")
WT<-WT[WT$SpacerSeq!="",]
RTD<-RTD[RTD$SpacerSeq!="",]
### 
WT$end5_2<-substr(WT$SpacerSeq,1,2)
WT$end5_3<-substr(WT$SpacerSeq,1,3)
WT$end3_2<-substr(WT$SpacerSeq,nchar(WT$SpacerSeq)-1,nchar(WT$SpacerSeq))
WT$end3_3<-substr(WT$SpacerSeq,nchar(WT$SpacerSeq)-2,nchar(WT$SpacerSeq))

RTD$end5_2<-substr(RTD$SpacerSeq,1,2)
RTD$end5_3<-substr(RTD$SpacerSeq,1,3)
RTD$end3_2<-substr(RTD$SpacerSeq,nchar(RTD$SpacerSeq)-1,nchar(RTD$SpacerSeq))
RTD$end3_3<-substr(RTD$SpacerSeq,nchar(RTD$SpacerSeq)-2,nchar(RTD$SpacerSeq))

pdf("spacer-encoded2.pdf",width=8,height=12)
par(mfrow=c(2,2))
df<-make_fq_df(WT$end5_2[nchar(WT$SSeq3)>2],RTD$end5_2[nchar(RTD$SSeq3)>2])
df<-merge(df,WT_sim5_2,by=1)
df<-merge(df,RTD_sim5_2,by=1)
df<-df[,c(1,2,4,3,5)]
colnames(df)[c(3,5)]<-c("WT_sim","RTD_sim")
mp<-barplot(100*prop.table(as.matrix(df[,-1]),2),col=P16,axes=F,
            names.arg=colSums(df[,-1]),
            width=0.15,
            space=0.3,
            xlim=c(0,.8),
            legend.text=df$NT,
            adj=0,args.legend=list(x=0.8,xjust=0,y=50,yjust=0.5,bty="n",cex=0.7))
mtext("Spacers(encoded nucleotides, 5'end, NEU>2nt)",line=2,at=median(mp),cex=0.7)
mtext(c("WT","WT-sim","RTD-sim","RTD"),line=0.25,at=mp,cex=0.7)
axis(1,labels=NA,at=c(0,0.8),lwd=1,lwd.ticks=0)
axis(1,labels=NA,at=mp,lwd=0,lwd.ticks=1)
axis(2,labels=c(0,25,50,75,100),las=1,at=c(0,25,50,75,100),pos=0,cex.axis=0.7)
par(xpd=T)
text(x=-0.15,y=50,"Spacers (%)",srt=90,adj=0.5,cex=0.7)

df<-make_fq_df(WT$end5_2[nchar(WT$SSeq3)<3],RTD$end5_2[nchar(RTD$SSeq3)<3])
df<-merge(df,WT_sim5_2,by=1)
df<-merge(df,RTD_sim5_2,by=1)
df<-df[,c(1,2,4,3,5)]
colnames(df)[c(3,5)]<-c("WT_sim","RTD_sim")
mp<-barplot(100*prop.table(as.matrix(df[,-1]),2),col=P16,axes=F,
            names.arg=colSums(df[,-1]),
            width=0.15,
            space=0.3,
            xlim=c(0,.8),
            legend.text=df$NT,
            adj=0,args.legend=list(x=0.8,xjust=0,y=50,yjust=0.5,bty="n",cex=0.7))
mtext("Spacers(encoded nucleotides, 5'end, NEU 1-2nt)",line=2,at=median(mp),cex=0.7)
mtext(c("WT","WT-sim","RTD-sim","RTD"),line=0.25,at=mp,cex=0.7)
axis(1,labels=NA,at=c(0,0.8),lwd=1,lwd.ticks=0)
axis(1,labels=NA,at=mp,lwd=0,lwd.ticks=1)
axis(2,labels=c(0,25,50,75,100),las=1,at=c(0,25,50,75,100),pos=0,cex.axis=0.7)
par(xpd=T)
text(x=-0.15,y=50,"Spacers (%)",srt=90,adj=0.5,cex=0.7)

df<-make_fq_df(WT$end5_3[nchar(WT$SSeq3)>=3],RTD$end5_3[nchar(RTD$SSeq3)>=3])
df<-merge(df,WT_sim5_3,by=1)
df<-merge(df,RTD_sim5_3,by=1)
df<-df[,c(1,2,4,3,5)]
colnames(df)[c(3,5)]<-c("WT_sim","RTD_sim")
df<-df[df$NT%in%c("AAA","CCC","GGG","TTT"),]
mp<-barplot(100*prop.table(as.matrix(df[,-1]),2),col=P64,axes=F,
            names.arg=colSums(df[,-1]),
            width=0.15,
            space=0.3,
            xlim=c(0,.8),
            legend.text=df$NT,
            adj=0,args.legend=list(x=0.8,xjust=0,y=50,yjust=0.5,bty="n",cex=0.7))
mtext("Spacers(encoded nucleotides, 5'end, NEU>2nt)",line=2,at=median(mp),cex=0.7)
mtext(c("WT","WT-sim","RTD-sim","RTD"),line=0.25,at=mp,cex=0.7)
axis(1,labels=NA,at=c(0,0.8),lwd=1,lwd.ticks=0)
axis(1,labels=NA,at=mp,lwd=0,lwd.ticks=1)
axis(2,labels=c(0,25,50,75,100),las=1,at=c(0,25,50,75,100),pos=0,cex.axis=0.7)
par(xpd=T)
text(x=-0.15,y=50,"Spacers (%)",srt=90,adj=0.5,cex=0.7)

df<-make_fq_df(WT$end5_3[nchar(WT$SSeq3)<3],RTD$end5_3[nchar(RTD$SSeq3)<3])
df<-merge(df,WT_sim5_3,by=1)
df<-merge(df,RTD_sim5_3,by=1)
df<-df[,c(1,2,4,3,5)]
colnames(df)[c(3,5)]<-c("WT_sim","RTD_sim")
df<-df[df$NT%in%c("AAA","CCC","GGG","TTT"),]
mp<-barplot(100*prop.table(as.matrix(df[,-1]),2),col=P64,axes=F,
            names.arg=colSums(df[,-1]),
            width=0.15,
            space=0.3,
            xlim=c(0,.8),
            legend.text=df$NT,
            adj=0,args.legend=list(x=0.8,xjust=0,y=50,yjust=0.5,bty="n",cex=0.7))
mtext("Spacers(encoded nucleotides, 5'end, NEU 1-2nt)",line=2,at=median(mp),cex=0.7)
mtext(c("WT","WT-sim","RTD-sim","RTD"),line=0.25,at=mp,cex=0.7)
axis(1,labels=NA,at=c(0,0.8),lwd=1,lwd.ticks=0)
axis(1,labels=NA,at=mp,lwd=0,lwd.ticks=1)
axis(2,labels=c(0,25,50,75,100),las=1,at=c(0,25,50,75,100),pos=0,cex.axis=0.7)
par(xpd=T)
text(x=-0.15,y=50,"Spacers (%)",srt=90,adj=0.5,cex=0.7)
df<-make_fq_df(WT$end3_2[nchar(WT$SSeq3)>2],RTD$end3_2[nchar(RTD$SSeq3)>2])
df<-merge(df,WT_sim3_2,by=1)
df<-merge(df,RTD_sim3_2,by=1)
df<-df[,c(1,2,4,3,5)]
colnames(df)[c(3,5)]<-c("WT_sim","RTD_sim")
mp<-barplot(100*prop.table(as.matrix(df[,-1]),2),col=P16,axes=F,
            names.arg=colSums(df[,-1]),
            width=0.15,
            space=0.3,
            xlim=c(0,.8),
            legend.text=df$NT,
            adj=0,args.legend=list(x=0.8,xjust=0,y=50,yjust=0.5,bty="n",cex=0.7))
mtext("Spacers(encoded nucleotides, 3'end, NEU>2nt)",line=2,at=median(mp),cex=0.7)
mtext(c("WT","WT-sim","RTD-sim","RTD"),line=0.25,at=mp,cex=0.7)
axis(1,labels=NA,at=c(0,0.8),lwd=1,lwd.ticks=0)
axis(1,labels=NA,at=mp,lwd=0,lwd.ticks=1)
axis(2,labels=c(0,25,50,75,100),las=1,at=c(0,25,50,75,100),pos=0,cex.axis=0.7)
par(xpd=T)
text(x=-0.15,y=50,"Spacers (%)",srt=90,adj=0.5,cex=0.7)

df<-make_fq_df(WT$end3_2[nchar(WT$SSeq3)<3],RTD$end3_2[nchar(RTD$SSeq3)<3])
df<-merge(df,WT_sim3_2,by=1)
df<-merge(df,RTD_sim3_2,by=1)
df<-df[,c(1,2,4,3,5)]
colnames(df)[c(3,5)]<-c("WT_sim","RTD_sim")
mp<-barplot(100*prop.table(as.matrix(df[,-1]),2),col=P16,axes=F,
            names.arg=colSums(df[,-1]),
            width=0.15,
            space=0.3,
            xlim=c(0,.8),
            legend.text=df$NT,
            adj=0,args.legend=list(x=0.8,xjust=0,y=50,yjust=0.5,bty="n",cex=0.7))
mtext("Spacers(encoded nucleotides, 3'end, NEU 1-2nt)",line=2,at=median(mp),cex=0.7)
mtext(c("WT","WT-sim","RTD-sim","RTD"),line=0.25,at=mp,cex=0.7)
axis(1,labels=NA,at=c(0,0.8),lwd=1,lwd.ticks=0)
axis(1,labels=NA,at=mp,lwd=0,lwd.ticks=1)
axis(2,labels=c(0,25,50,75,100),las=1,at=c(0,25,50,75,100),pos=0,cex.axis=0.7)
par(xpd=T)
text(x=-0.15,y=50,"Spacers (%)",srt=90,adj=0.5,cex=0.7)

df<-make_fq_df(WT$end3_3[nchar(WT$SSeq3)>=3],RTD$end3_3[nchar(RTD$SSeq3)>=3])
df<-merge(df,WT_sim3_3,by=1)
df<-merge(df,RTD_sim3_3,by=1)
df<-df[,c(1,2,4,3,5)]
colnames(df)[c(3,5)]<-c("WT_sim","RTD_sim")
df<-df[df$NT%in%c("AAA","CCC","GGG","TTT"),]
mp<-barplot(100*prop.table(as.matrix(df[,-1]),2),col=P64,axes=F,
            names.arg=colSums(df[,-1]),
            width=0.15,
            space=0.3,
            xlim=c(0,.8),
            legend.text=df$NT,
            adj=0,args.legend=list(x=0.8,xjust=0,y=50,yjust=0.5,bty="n",cex=0.7))
mtext("Spacers(encoded nucleotides, 3'end, NEU>2nt)",line=2,at=median(mp),cex=0.7)
mtext(c("WT","WT-sim","RTD-sim","RTD"),line=0.25,at=mp,cex=0.7)
axis(1,labels=NA,at=c(0,0.8),lwd=1,lwd.ticks=0)
axis(1,labels=NA,at=mp,lwd=0,lwd.ticks=1)
axis(2,labels=c(0,25,50,75,100),las=1,at=c(0,25,50,75,100),pos=0,cex.axis=0.7)
par(xpd=T)
text(x=-0.15,y=50,"Spacers (%)",srt=90,adj=0.5,cex=0.7)

df<-make_fq_df(WT$end3_3[nchar(WT$SSeq3)<3],RTD$end3_3[nchar(RTD$SSeq3)<3])
df<-merge(df,WT_sim3_3,by=1)
df<-merge(df,RTD_sim3_3,by=1)
df<-df[,c(1,2,4,3,5)]
colnames(df)[c(3,5)]<-c("WT_sim","RTD_sim")
df<-df[df$NT%in%c("AAA","CCC","GGG","TTT"),]
mp<-barplot(100*prop.table(as.matrix(df[,-1]),2),col=P64,axes=F,
            names.arg=colSums(df[,-1]),
            width=0.15,
            space=0.3,
            xlim=c(0,.8),
            legend.text=df$NT,
            adj=0,args.legend=list(x=0.8,xjust=0,y=50,yjust=0.5,bty="n",cex=0.7))
mtext("Spacers(encoded nucleotides, 3'end, NEU 1-2nt)",line=2,at=median(mp),cex=0.7)
mtext(c("WT","WT-sim","RTD-sim","RTD"),line=0.25,at=mp,cex=0.7)
axis(1,labels=NA,at=c(0,0.8),lwd=1,lwd.ticks=0)
axis(1,labels=NA,at=mp,lwd=0,lwd.ticks=1)
axis(2,labels=c(0,25,50,75,100),las=1,at=c(0,25,50,75,100),pos=0,cex.axis=0.7)
par(xpd=T)
text(x=-0.15,y=50,"Spacers (%)",srt=90,adj=0.5,cex=0.7)
dev.off()


### AT rich at ends?----------------------------------------------------------------------------------------------------------
##arr1
WT<-read.delim("WT.rearr1")
RTD<-read.delim("RTD.rearr1")
### 
WT$SpacerSeq[WT$tail==0]<-WT$Seq[WT$tail==0]
WT$end5_3<-substr(WT$SpacerSeq,1,3)
WT$end3_3<-substr(WT$SpacerSeq,nchar(WT$SpacerSeq)-2,nchar(WT$SpacerSeq))
RTD$SpacerSeq[RTD$tail==0]<-RTD$Seq[RTD$tail==0]
RTD$end5_3<-substr(RTD$SpacerSeq,1,3)
RTD$end3_3<-substr(RTD$SpacerSeq,nchar(RTD$SpacerSeq)-2,nchar(RTD$SpacerSeq))

pdf("spacer-encoded1.pdf",width=8,height=12)
par(mfrow=c(2,1),mar=c(2,4,5,5))
df<-make_fq_df(WT$end5_3[WT$tail==0],RTD$end5_3[RTD$tail==0])
df<-merge(df,make_fq_df(WT$end5_3[WT$tail==1],RTD$end5_3[RTD$tail==1]),by=1)
df<-merge(df,make_fq_df(WT$end5_3[WT$tail==2],RTD$end5_3[RTD$tail==2]),by=1)
df<-merge(df,make_fq_df(WT$end5_3[WT$tail==3],RTD$end5_3[RTD$tail==3]),by=1)
df<-merge(df,make_fq_df(WT$end5_3[WT$tail==4],RTD$end5_3[RTD$tail==4]),by=1)
df<-merge(df,make_fq_df(WT$end5_3[WT$tail>4],RTD$end5_3[RTD$tail>4]),by=1)
df$AT<-str_count(df$NT,"A")+str_count(df$NT,"T")
agg<-aggregate(.~AT,data=df[,-1],sum)
mp<-barplot(100*prop.table(as.matrix(agg[,-1]),2),col=P16,axes=F,
            names.arg=colSums(agg[,-1]),cex.names=0.7,
            width=0.15,
            space=0.3,
            xlim=c(0,2.3),
            legend.text=agg$AT,
            adj=0,args.legend=list(x=2.35,xjust=0,y=50,yjust=0.5,bty="n",cex=0.7,title = 'Number of A/T'))
mtext("5' end trinucleotides",line=2,at=median(mp))
mtext(rep(c("WT","RTD"),each=6),line=0.25,at=mp,cex=0.7)
mtext(c("0","1","2","3","4",">=5"),line=1,
      at=c(mean(mp[1:2]),mean(mp[3:4]),mean(mp[5:6]),mean(mp[7:8]),mean(mp[9:10]),mean(mp[11:12])),cex=0.7)
mtext("Non-encoded nucleotides",line=1,at=-0.1,cex=0.7)
axis(1,labels=NA,at=c(0,2.4),lwd=1,lwd.ticks=0)
axis(1,labels=NA,at=mp,lwd=0,lwd.ticks=1)
axis(2,labels=c(0,25,50,75,100),las=1,at=c(0,25,50,75,100),pos=0,cex.axis=0.7)
par(xpd=T)
segments(mp[1]-0.075,104.5,mp[2]+0.075,104.5)
segments(mp[3]-0.075,104.5,mp[4]+0.075,104.5)
segments(mp[5]-0.075,104.5,mp[6]+0.075,104.5)
segments(mp[7]-0.075,104.5,mp[8]+0.075,104.5)
segments(mp[9]-0.075,104.5,mp[10]+0.075,104.5)
segments(mp[11]-0.075,104.5,mp[12]+0.075,104.5)
text(x=-0.15,y=50,"Spacers (%)",srt=90,adj=0.5,cex=0.7)

df<-make_fq_df(WT$end3_3[WT$tail==0],RTD$end3_3[RTD$tail==0])
df<-merge(df,make_fq_df(WT$end3_3[WT$tail==1],RTD$end3_3[RTD$tail==1]),by=1)
df<-merge(df,make_fq_df(WT$end3_3[WT$tail==2],RTD$end3_3[RTD$tail==2]),by=1)
df<-merge(df,make_fq_df(WT$end3_3[WT$tail==3],RTD$end3_3[RTD$tail==3]),by=1)
df<-merge(df,make_fq_df(WT$end3_3[WT$tail==4],RTD$end3_3[RTD$tail==4]),by=1)
df<-merge(df,make_fq_df(WT$end3_3[WT$tail>4],RTD$end3_3[RTD$tail>4]),by=1)
df$AT<-str_count(df$NT,"A")+str_count(df$NT,"T")
agg<-aggregate(.~AT,data=df[,-1],sum)
mp<-barplot(100*prop.table(as.matrix(agg[,-1]),2),col=P16,axes=F,
            names.arg=colSums(agg[,-1]),cex.names=0.7,
            width=0.15,
            space=0.3,
            xlim=c(0,2.3),
            legend.text=agg$AT,
            adj=0,args.legend=list(x=2.35,xjust=0,y=50,yjust=0.5,bty="n",cex=0.7,title = 'Number of A/T'))
mtext("3' end trinucleotides",line=2,at=median(mp))
mtext(rep(c("WT","RTD"),6),line=0.25,at=mp,cex=0.7)
mtext(c("0","1","2","3","4",">=5"),line=1,
      at=c(mean(mp[1:2]),mean(mp[3:4]),mean(mp[5:6]),mean(mp[7:8]),mean(mp[9:10]),mean(mp[11:12])),cex=0.7)
mtext("Non-encoded nucleotides",line=1,at=-0.1,cex=0.7)
axis(1,labels=NA,at=c(0,2.4),lwd=1,lwd.ticks=0)
axis(1,labels=NA,at=mp,lwd=0,lwd.ticks=1)
axis(2,labels=c(0,25,50,75,100),las=1,at=c(0,25,50,75,100),pos=0,cex.axis=0.7)
par(xpd=T)
segments(mp[1]-0.075,104.5,mp[2]+0.075,104.5)
segments(mp[3]-0.075,104.5,mp[4]+0.075,104.5)
segments(mp[5]-0.075,104.5,mp[6]+0.075,104.5)
segments(mp[7]-0.075,104.5,mp[8]+0.075,104.5)
segments(mp[9]-0.075,104.5,mp[10]+0.075,104.5)
segments(mp[11]-0.075,104.5,mp[12]+0.075,104.5)
text(x=-0.15,y=50,"Spacers (%)",srt=90,adj=0.5,cex=0.7)
dev.off()