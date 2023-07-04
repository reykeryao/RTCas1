rm(list=ls())
setwd("/stor/work/Lambowitz/yaojun/Work/RTCas1/Documents/Script/")
library(RColorBrewer)
library(stringi)
col<-brewer.pal(6,"Paired")[c(1,3,5,2,4,6)]

### simulate WT/RTD spacers with random nucleotides with similar length and ACGT
### distribution, 10k sequences each
WT<-read.delim("MMB1_WT_sp_detail.info")
RTD<-read.delim("MMB1_RTdelta_sp_detail.info")

### WT
### 6701 antisense spacer, 54 undecided, 8743 sense spacer
### WT ignore "U", flip the antisense spacer to sense spacer.
WT$SpacerSeq<-""
for (i in 1:dim(WT)[1]){
  #  if (WT[i,16]!="A"){
  WT[i,17]<-substr(WT[i,3],WT[i,4]+1,nchar(WT[i,3])-WT[i,5])
  #  } else if (WT[i,16]=="A"){
  #   WT[i,17]<-rev.comp(substr(WT[i,3],WT[i,4]+1,nchar(WT[i,3])-WT[i,5]))
  #  }
}

RTD$SpacerSeq<-""
for (i in 1:dim(RTD)[1]){
  RTD[i,17]<-substr(RTD[i,3],RTD[i,4]+1,nchar(RTD[i,3])-RTD[i,5])
}

WT$Sp_len<-nchar(WT$SpacerSeq)
RTD$Sp_len<-nchar(RTD$SpacerSeq)

### simulate sequences
len_prob<-data.frame(table(WT$Sp_len))
colnames(len_prob)<-c("Len","WT")
                     
len_prob<-merge(len_prob,data.frame(table(RTD$Sp_len)),by=1,all=T)
len_prob$Len<-as.integer(as.character(len_prob$Len))
len_prob<-merge(data.frame("Len"=c(min(len_prob$Len):max(len_prob$Len))),len_prob,by=1,all=T)
len_prob[is.na(len_prob)]<-0
colnames(len_prob)[3]<-"RTD"
len_prob[,2:3]<-prop.table(as.matrix(len_prob[,2:3]),2)

tmp<-paste0(WT$SpacerSeq,collapse = "")
nt_prob<-data.frame(c(stri_count(tmp,fixed="A"),
                      stri_count(tmp,fixed="C"),
                      stri_count(tmp,fixed="G"),
                      stri_count(tmp,fixed="T")))
tmp<-paste0(RTD$SpacerSeq,collapse = "")
nt_prob<-cbind(nt_prob,data.frame(c(stri_count(tmp,fixed="A"),
                      stri_count(tmp,fixed="C"),
                      stri_count(tmp,fixed="G"),
                      stri_count(tmp,fixed="T"))))
rm(tmp)
colnames(nt_prob)<-c("WT","RTD")                    
rownames(nt_prob)<-c("A","C","G","T")
nt_prob<-data.frame(prop.table(as.matrix(nt_prob),2))
nt_prob$NT<-rownames(nt_prob)
sim_num<-10000
for (rep in 1:100){
  sim_sapcer<-data.frame("WT_len"=sample(len_prob$Len,sim_num,replace = T,prob = len_prob$WT),
                         "RTD_len"=sample(len_prob$Len,sim_num,replace = T,prob = len_prob$RTD))
  sim_sapcer$WT_sp<-""
  sim_sapcer$RTD_sp<-""
  sim_sapcer$WT_sp_withRTDlen<-""
  sim_sapcer$RTD_sp_withWTlen<-""
  for (i in 1:sim_num){
    sp_len<-sim_sapcer$WT_len[i]
    sp_seq<-sample(nt_prob$NT,sp_len,replace = T,prob = nt_prob$WT)
    sim_sapcer$WT_sp[i]<-paste(sp_seq,collapse = "")
    sp_seq<-sample(nt_prob$NT,sp_len,replace = T,prob = nt_prob$RTD)
    sim_sapcer$RTD_sp_withWTlen[i]<-paste(sp_seq,collapse = "")
    
    sp_len<-sim_sapcer$RTD_len[i]
    sp_seq<-sample(nt_prob$NT,sp_len,replace = T,prob = nt_prob$RTD)
    sim_sapcer$RTD_sp[i]<-paste(sp_seq,collapse = "")
    sp_seq<-sample(nt_prob$NT,sp_len,replace = T,prob = nt_prob$WT)
    sim_sapcer$WT_sp_withRTDlen[i]<-paste(sp_seq,collapse = "")
  }
  sim_sapcer$ID<-paste0("SimSP",1:sim_num)
  write.table(sim_sapcer,paste0("sim_sp/simulated_spacer",rep),quote=F,sep="\t",row.names=F)
  if (rep%%10==0){
    print(paste0("Finished:",rep,"%"))
    }
}


### convert to fa 
### 
### and run RNAfold
### RNAfold --noPS -i WT_sim10k.fa |tr \\n \\t|tr ">" \\n| \
###  sed 's/ (/ /g'|cut -f 1,2,4|sed 's/[ |)]//g' > WT_sim10k.mfe

WTMFE<-read.table("MMB1_WT_spacer_noSS_forMFE.mfe",sep="\t",col.names = c("ID","Seq","MFE"))
RTDMFE<-read.table("MMB1_RTD_spacer_noSS_forMFE.mfe",sep="\t",col.names = c("ID","Seq","MFE"))

RTD_col="#FF000010"
pdf("FigS7B.pdf",width=10,height=5)
par(mfrow=c(1,2))
plot(density(WTMFE$MFE),bty="n",xlab="MFE (kcal/mol)",main="Spacer w/o SS",ylim=c(0,0.16))
lines(density(RTDMFE$MFE),col="red")
for (i in 1:100){
  WTSIM<-read.table(paste0("sim_sp/sim_mfe/WT_sim",i,".mfe"),
                    sep="\t",col.names = c("ID","Seq","MFE"))
  RTDSIM<-read.table(paste0("sim_sp/sim_mfe/RTD_sim",i,".mfe"),
                     sep="\t",col.names = c("ID","Seq","MFE"))
  lines(density(WTSIM$MFE),lty=2,col="gray75")
  lines(density(RTDSIM$MFE),col=RTD_col,lty=2)
}
lines(density(WTMFE$MFE),lwd=2)
lines(density(RTDMFE$MFE),col="red",lwd=3)
legend(-18,0.12,col=rep(c("black","red"),each=2),
       lty=c(1,2,1,2),bty="n",
       legend = c("WT","simulated WT","RTD","simulated RTD"))
#dev.off()


plot(ecdf(WTMFE$MFE),xlab="MFE (kcal/mol)",ylab="ECDF",
     do.points=F,col.01line = NULL,main=NA,bty="n",
     xlim=c(-10.5,-0.4))
lines(ecdf(RTDMFE$MFE),col="red",do.points=F,col.01line = NULL)
for (i in 1:100){
  WTSIM<-read.table(paste0("sim_sp/sim_mfe/WT_sim",i,".mfe"),
                    sep="\t",col.names = c("ID","Seq","MFE"))
  RTDSIM<-read.table(paste0("sim_sp/sim_mfe/RTD_sim",i,".mfe"),
                     sep="\t",col.names = c("ID","Seq","MFE"))
  lines(ecdf(WTSIM$MFE),col="gray75",do.points=F,col.01line = NULL,lty=2)
  lines(ecdf(RTDSIM$MFE),col=RTD_col,do.points=F,col.01line = NULL,lty=2)
}
lines(ecdf(WTMFE$MFE),do.points=F,col.01line = NULL,lwd=3)
lines(ecdf(RTDMFE$MFE),col="red",do.points=F,col.01line = NULL,lwd=3)
abline(h=0.5,lty=3)

dev.off()