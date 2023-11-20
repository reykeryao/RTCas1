rm(list=ls())
setwd("/stor/work/Lambowitz/yaojun/Work/RTCas1/Documents/Script/")
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
    RTD[i,18]<-RTD[i,8]
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
### chisq test of ss nt composition
data_df<-WT[WT$Sp_type!="U",]
tmp<-data.frame("Nucl"=c("A","G","C","T"))
for (i in 1:5){
  tmp<-merge(tmp,data.frame(table(substr(data_df$END,i,i))),by=1,all.x=T)
}
colnames(tmp)[2:6]<-c("1","2","3","4","5")
WT_SS<-tmp
data_df<-RTD[RTD$Sp_type!="U",]
tmp<-data.frame("Nucl"=c("A","G","C","T"))
for (i in 1:5){
  tmp<-merge(tmp,data.frame(table(substr(data_df$END,i,i))),by=1,all.x=T)
}
colnames(tmp)[2:6]<-c("1","2","3","4","5")
RTD_SS<-tmp
NT<-c("A","C","G","T")
for (i in 1:4){
  sample<-NT[i]
  for (j in 2:6){
    a<-matrix(c(WT_SS[i,j],RTD_SS[i,j],
                   sum(WT_SS[-i,j]),sum(RTD_SS[-i,j])),
                 nrow=2,ncol=2)
    if ((a[1,1]/sum(a[1,]))>(a[2,1]/sum(a[2,]))){
      p<-fisher.test(a,simulate.p.value = T,alternative = "greater")$p.value
      if (p<0.001){
        print(paste0("p-value of WT has more ",sample, " at position ",j-1, " is <0.001"))
      } else {
        p<-round(p,3)
        print(paste0("p-value of WT has more ",sample, " at position ",j-1, " is ",p))
      }
      
    } else {
      p<-fisher.test(a,simulate.p.value = T,alternative = "less")$p.value
      if (p<0.001){
        print(paste0("p-value of WT has less ",sample, " at position ",j-1, " is <0.001"))
      } else {
        p<-round(p,3)
        print(paste0("p-value of WT has less ",sample, " at position ",j-1, " is ",p))
      }
    }
  }
}


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


## Fig1D-ver2 based on rearrangement 2
pdf("Fig1D-2.pdf",width=8,height=6)
par(mfrow=c(2,2),cex=0.7)
#from the start of the soft-clipped seq (+1 to +5 from the end of the sapcer)
plotspcar_stack5<-function(data_df,title){
  tmp<-data.frame("Nucl"=c("A","G","C","T"))
  for (i in 1:5){
    tmp<-merge(tmp,data.frame(table(substr(data_df$SSeq5,nchar(data_df$SSeq5)+1-i,nchar(data_df$SSeq5)+1-i))),
               by=1,all.x=T)
  }
  y_max<-round(max(colSums(tmp[,2:6]))/500)*500
  tmp<-tmp[,c(1,6,5,4,3,2)]
  colnames(tmp)[2:6]<-c("-5","-4","-3","-2","-1")
  df.bar <-barplot(as.matrix(tmp[,2:6]),ylab="Spacers", 
                   main=paste0(title," (5' end non-encoded nucleotides)"),
                   width = 0.1, ylim=c(0,y_max),
                   beside = F,col=c("green","blue","chocolate","red"),
                   args.legend = list(x="right",bty="n"),legend.text = tmp$Nucl,
                   space=0.1,xlim=c(0,0.7),xlab="Position from the begining of spacer")
}

plotspcar_stack3<-function(data_df,title){
  tmp<-data.frame("Nucl"=c("A","G","C","T"))
  for (i in 1:5){
    tmp<-merge(tmp,data.frame(table(substr(data_df$SSeq3,i,i))),by=1,all.x=T)
  }
  y_max<-round(max(colSums(tmp[,2:6]))/500)*500
  colnames(tmp)[2:6]<-c("1","2","3","4","5")
  df.bar <-barplot(as.matrix(tmp[,2:6]),ylab="Spacers", 
                   main=paste0(title," (3' end non-encoded nucleotides)"),
                   width = 0.1, ylim=c(0,y_max),
                   beside = F,col=c("green","blue","chocolate","red"),
                   args.legend = list(x="right",bty="n"),legend.text = tmp$Nucl,
                   space=0.1,xlim=c(0,0.7),xlab="Position from the end of spacer")
}
df1<-WT[WT$Sp_type!="U",]
plotspcar_stack5(df1,"WT")
plotspcar_stack3(df1,"WT")

df1<-RTD[RTD$Sp_type!="U",]
plotspcar_stack5(df1,"RTdelat")
plotspcar_stack3(df1,"RTdelat")
dev.off()

### chisq test of ss nt composition
data_df<-WT[WT$Sp_type!="U",]
tmp<-data.frame("Nucl"=c("A","G","C","T"))
for (i in 1:5){
  tmp<-merge(tmp,data.frame(table(substr(data_df$SSeq5,nchar(data_df$SSeq5)+1-i,nchar(data_df$SSeq5)+1-i))),
             by=1,all.x=T)
}
colnames(tmp)[2:6]<-c("1","2","3","4","5")
WT_SS5<-tmp
tmp<-data.frame("Nucl"=c("A","G","C","T"))
for (i in 1:5){
  tmp<-merge(tmp,data.frame(table(substr(data_df$SSeq3,i,i))),by=1,all.x=T)
}
colnames(tmp)[2:6]<-c("1","2","3","4","5")
tmp$Nucl<-sapply(tmp$Nucl,rev.comp)
tmp<-tmp[order(tmp$Nucl),]
WT_SS3<-tmp

## 
NT<-c("A","C","G","T")
for (i in 1:4){
  sample<-NT[i]
  for (j in 2:6){
    a<-matrix(c(WT_SS5[i,j],WT_SS3[i,j],
                sum(WT_SS5[-i,j]),sum(WT_SS3[-i,j])),
              nrow=2,ncol=2)
    if ((a[1,1]/sum(a[1,]))>(a[2,1]/sum(a[2,]))){
      p<-fisher.test(a,simulate.p.value = T,alternative = "greater")$p.value
      if (p<0.001){
        print(paste0("p-value of WT has more ",sample, " at position ",j-1, " is <0.001"))
      } else {
        p<-round(p,3)
        print(paste0("p-value of WT has more ",sample, " at position ",j-1, " is ",p))
      }
      
    } else {
      p<-fisher.test(a,simulate.p.value = T,alternative = "less")$p.value
      if (p<0.001){
        print(paste0("p-value of WT has less ",sample, " at position ",j-1, " is <0.001"))
      } else {
        p<-round(p,3)
        print(paste0("p-value of WT has less ",sample, " at position ",j-1, " is ",p))
      }
    }
  }
}

data_df<-RTD[RTD$Sp_type!="U",]
tmp<-data.frame("Nucl"=c("A","G","C","T"))
for (i in 1:5){
  tmp<-merge(tmp,data.frame(table(substr(data_df$SSeq5,nchar(data_df$SSeq5)+1-i,nchar(data_df$SSeq5)+1-i))),
             by=1,all.x=T)
}
colnames(tmp)[2:6]<-c("1","2","3","4","5")
RTD_SS5<-tmp
tmp<-data.frame("Nucl"=c("A","G","C","T"))
for (i in 1:5){
  tmp<-merge(tmp,data.frame(table(substr(data_df$SSeq3,i,i))),by=1,all.x=T)
}
colnames(tmp)[2:6]<-c("1","2","3","4","5")
tmp$Nucl<-sapply(tmp$Nucl,rev.comp)
tmp<-tmp[order(tmp$Nucl),]
RTD_SS3<-tmp
## 
NT<-c("A","C","G","T")
for (i in 1:4){
  sample<-NT[i]
  for (j in 2:6){
    a<-matrix(c(RTD_SS5[i,j],RTD_SS3[i,j],
                sum(RTD_SS5[-i,j]),sum(RTD_SS3[-i,j])),
              nrow=2,ncol=2)
    if ((a[1,1]/sum(a[1,]))>(a[2,1]/sum(a[2,]))){
      p<-fisher.test(a,simulate.p.value = T,alternative = "greater")$p.value
      if (p<0.001){
        print(paste0("p-value of WT has more ",sample, " at position ",j-1, " is <0.001"))
      } else {
        p<-round(p,3)
        print(paste0("p-value of WT has more ",sample, " at position ",j-1, " is ",p))
      }
      
    } else {
      p<-fisher.test(a,simulate.p.value = T,alternative = "less")$p.value
      if (p<0.001){
        print(paste0("p-value of WT has less ",sample, " at position ",j-1, " is <0.001"))
      } else {
        p<-round(p,3)
        print(paste0("p-value of WT has less ",sample, " at position ",j-1, " is ",p))
      }
    }
  }
}



NT<-c("A","C","G","T")
for (i in 1:4){
  sample<-NT[i]
  for (j in 2:6){
    a<-matrix(c(WT_SS5[i,j],RTD_SS5[i,j],
                sum(WT_SS5[-i,j]),sum(RTD_SS5[-i,j])),
              nrow=2,ncol=2)
    if ((a[1,1]/sum(a[1,]))>(a[2,1]/sum(a[2,]))){
      p<-fisher.test(a,simulate.p.value = T,alternative = "greater")$p.value
      if (p<0.001){
        print(paste0("p-value of WT has more ",sample, " at position ",j-1, " is <0.001"))
      } else {
        p<-round(p,3)
        print(paste0("p-value of WT has more ",sample, " at position ",j-1, " is ",p))
      }
      
    } else {
      p<-fisher.test(a,simulate.p.value = T,alternative = "less")$p.value
      if (p<0.001){
        print(paste0("p-value of WT has less ",sample, " at position ",j-1, " is <0.001"))
      } else {
        p<-round(p,3)
        print(paste0("p-value of WT has less ",sample, " at position ",j-1, " is ",p))
      }
    }
  }
}

NT<-c("A","C","G","T")
for (i in 1:4){
  sample<-NT[i]
  for (j in 2:6){
    a<-matrix(c(WT_SS3[i,j],RTD_SS3[i,j],
                sum(WT_SS3[-i,j]),sum(RTD_SS3[-i,j])),
              nrow=2,ncol=2)
    if ((a[1,1]/sum(a[1,]))>(a[2,1]/sum(a[2,]))){
      p<-fisher.test(a,simulate.p.value = T,alternative = "greater")$p.value
      if (p<0.001){
        print(paste0("p-value of WT has more ",sample, " at position ",j-1, " is <0.001"))
      } else {
        p<-round(p,3)
        print(paste0("p-value of WT has more ",sample, " at position ",j-1, " is ",p))
      }
      
    } else {
      p<-fisher.test(a,simulate.p.value = T,alternative = "less")$p.value
      if (p<0.001){
        print(paste0("p-value of WT has less ",sample, " at position ",j-1, " is <0.001"))
      } else {
        p<-round(p,3)
        print(paste0("p-value of WT has less ",sample, " at position ",j-1, " is ",p))
      }
    }
  }
}



### arrange 3
### encoded spacers are always in sense orientation
### isolate spacers mapped to know gene (can define sense/antisense) 
### and with SS in at only one end
### but RTD is not rearranged, it's as
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
  if (RTD$SS5[i]>0){
    RTD$SSeq5[i]<-RTD$SS5_seq[i]
    RTD$SpacerSeq[i]<-substr(RTD$Seq[i],RTD$SS5[i]+1,nchar(RTD$Seq[i]))
  } else if (RTD$SS3[i]>0) {
    RTD$SpacerSeq[i]<-substr(RTD$Seq[i],1,nchar(RTD$Seq[i])-RTD$SS3[i])
    RTD$SSeq3[i]<-RTD$SS3_seq[i]
  } else {
    RTD$SpacerSeq[i]<-RTD$Seq[i]
  }
}


## Fig1D-ver3 based on rearrangement 3
pdf("Fig1D-3.pdf",width=8,height=6)
par(mfrow=c(2,2),cex=0.7)
#from the start of the soft-clipped seq (+1 to +5 from the end of the sapcer)
plotspcar_stack5<-function(data_df,title){
  tmp<-data.frame("Nucl"=c("A","G","C","T"))
  for (i in 1:5){
    tmp<-merge(tmp,data.frame(table(substr(data_df$SSeq5,nchar(data_df$SSeq5)+1-i,nchar(data_df$SSeq5)+1-i))),
               by=1,all.x=T)
  }
  y_max<-round(max(colSums(tmp[,2:6]))/500)*500
  tmp<-tmp[,c(1,6,5,4,3,2)]
  colnames(tmp)[2:6]<-c("-5","-4","-3","-2","-1")
  df.bar <-barplot(as.matrix(tmp[,2:6]),ylab="Spacers", 
                   main=paste0(title," (5' end non-encoded nucleotides)"),
                   width = 0.1, ylim=c(0,y_max),
                   beside = F,col=c("green","blue","chocolate","red"),
                   args.legend = list(x="right",bty="n"),legend.text = tmp$Nucl,
                   space=0.1,xlim=c(0,0.7),xlab="Position from the begining of spacer")
}

plotspcar_stack3<-function(data_df,title){
  tmp<-data.frame("Nucl"=c("A","G","C","T"))
  for (i in 1:5){
    tmp<-merge(tmp,data.frame(table(substr(data_df$SSeq3,i,i))),by=1,all.x=T)
  }
  y_max<-round(max(colSums(tmp[,2:6]))/500)*500
  colnames(tmp)[2:6]<-c("1","2","3","4","5")
  df.bar <-barplot(as.matrix(tmp[,2:6]),ylab="Spacers", 
                   main=paste0(title," (3' end non-encoded nucleotides)"),
                   width = 0.1, ylim=c(0,y_max),
                   beside = F,col=c("green","blue","chocolate","red"),
                   args.legend = list(x="right",bty="n"),legend.text = tmp$Nucl,
                   space=0.1,xlim=c(0,0.7),xlab="Position from the end of spacer")
}
df1<-WT[WT$Sp_type!="U",]
plotspcar_stack5(df1,"WT")
plotspcar_stack3(df1,"WT")

df1<-RTD[RTD$Sp_type!="U",]
plotspcar_stack5(df1,"RTdelat")
plotspcar_stack3(df1,"RTdelat")
dev.off()


###
data_df<-RTD[RTD$Sp_type!="U",]
tmp<-data.frame("Nucl"=c("A","G","C","T"))
for (i in 1:5){
  tmp<-merge(tmp,data.frame(table(substr(data_df$SSeq5,nchar(data_df$SSeq5)+1-i,nchar(data_df$SSeq5)+1-i))),
             by=1,all.x=T)
}
colnames(tmp)[2:6]<-c("1","2","3","4","5")
RTD_SS5_2<-tmp
tmp<-data.frame("Nucl"=c("A","G","C","T"))
for (i in 1:5){
  tmp<-merge(tmp,data.frame(table(substr(data_df$SSeq3,i,i))),by=1,all.x=T)
}
colnames(tmp)[2:6]<-c("1","2","3","4","5")
tmp$Nucl<-sapply(tmp$Nucl,rev.comp)
tmp<-tmp[order(tmp$Nucl),]
RTD_SS3_2<-tmp

### RTD_SS5_2 vs SS3_2
NT<-c("A","C","G","T")
for (i in 1:4){
  sample<-NT[i]
  for (j in 2:6){
    a<-matrix(c(RTD_SS5_2[i,j],RTD_SS3_2[i,j],
                sum(RTD_SS5_2[-i,j]),sum(RTD_SS3_2[-i,j])),
              nrow=2,ncol=2)
    if ((a[1,1]/sum(a[1,]))>(a[2,1]/sum(a[2,]))){
      p<-fisher.test(a,simulate.p.value = T,alternative = "greater")$p.value
      if (p<0.001){
        print(paste0("p-value of WT has more ",sample, " at position ",j-1, " is <0.001"))
      } else {
        p<-round(p,3)
        print(paste0("p-value of WT has more ",sample, " at position ",j-1, " is ",p))
      }
      
    } else {
      p<-fisher.test(a,simulate.p.value = T,alternative = "less")$p.value
      if (p<0.001){
        print(paste0("p-value of WT has less ",sample, " at position ",j-1, " is <0.001"))
      } else {
        p<-round(p,3)
        print(paste0("p-value of WT has less ",sample, " at position ",j-1, " is ",p))
      }
    }
  }
}

for (i in 1:4){
  sample<-NT[i]
  for (j in 2:6){
    a<-matrix(c(RTD_SS3_2[i,j],RTD_SS3[i,j],
                sum(RTD_SS3_2[-i,j]),sum(RTD_SS3[-i,j])),
              nrow=2,ncol=2)
    if ((a[1,1]/sum(a[1,]))>(a[2,1]/sum(a[2,]))){
      p<-fisher.test(a,simulate.p.value = T,alternative = "greater")$p.value
      if (p<0.001){
        print(paste0("p-value of WT has more ",sample, " at position ",j-1, " is <0.001"))
      } else {
        p<-round(p,3)
        print(paste0("p-value of WT has more ",sample, " at position ",j-1, " is ",p))
      }
      
    } else {
      p<-fisher.test(a,simulate.p.value = T,alternative = "less")$p.value
      if (p<0.001){
        print(paste0("p-value of WT has less ",sample, " at position ",j-1, " is <0.001"))
      } else {
        p<-round(p,3)
        print(paste0("p-value of WT has less ",sample, " at position ",j-1, " is ",p))
      }
    }
  }
}





### Fig1E not affected by rearrangement

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
