rm(list=ls())
gc()
library(tidyr)
library(RColorBrewer)
setwd("/stor/work/Lambowitz/yaojun/Work/RTCas1/JA23020/")
args = commandArgs(trailingOnly=TRUE)
name<-args[1]
ID<-read.delim("meta.info")
ID$Mn<-paste0("Mn",ID$Mn)
ID$Primer[ID$Primer=="-"]<-"no primer"
ID<-ID[ID$ID==name,]
#MAFFT ana, new should be used, better but take longer time
## only need to run once.
for (sample in 1:length(ID$ID)){
  print(paste0("Start processing Dataset ",ID$ID[sample]))
  dat<-read.table(paste0(ID$ID[sample],".output"),
                  col.names=c("ID","Align"))
  dat<-separate(dat,ID,into=c("ID","Reads","Seq"),sep="@",convert=T)
  dat<-separate(dat,Align,into=c("Empty",paste0("Pos",1:ID$Len[sample])),sep="")
  dat<-dat[-4]
  #length of the reads
  dat$Read_len<-0
  #length of MAFFT alignment
  dat$Align_len<-0
  #start position of alignment, the first matched nucleotides
  dat$Start<-6
  #end position of alignment, the last matched nucleotides
  dat$End<-ID$Len[sample]-5
  #5'end NTA
  dat$SS5<-""
  #3'end NTA
  dat$SS3<-""
  #number of matched nucleotides to template
  dat$Match<-0
  #err rate (NTA and removed insertion in reads not included, this is only for mismatch in the middle)
  dat$Err_rate<-0
  #Is 5'end of the read match the first aligned nucleotide in read, this is check if there is any truncation at 5'end
  dat$Match5<-""
  #Is 3'end of the read match the first aligned nucleotide in read, this is check if there is any truncation at 3'end
  dat$Match3<-""
  #insertion that removed by MAFFT, only meaningful when Match5 and Match3 are true
  dat$Insertion<-0
  #Is it a good read that should be included in analysis
  dat$Good<-""
  FLAG=0
  for (i in 2:dim(dat)[1]){
    if (i/dim(dat)[1]>=0.1 & FLAG==0){
      print(paste0(ID$ID[sample],": 10% finished"))
      FLAG<-FLAG+1
      }
    if (i/dim(dat)[1]>=0.2 & FLAG==1){
      print(paste0(ID$ID[sample],": 20% finished"))
      FLAG<-FLAG+1
    }
    if (i/dim(dat)[1]>=0.3 & FLAG==2){
      print(paste0(ID$ID[sample],": 30% finished"))
      FLAG<-FLAG+1
    }
    if (i/dim(dat)[1]>=0.4 & FLAG==3){
      print(paste0(ID$ID[sample],": 40% finished"))
      FLAG<-FLAG+1
    }
    if (i/dim(dat)[1]>=0.5 & FLAG==4){
      print(paste0(ID$ID[sample],": 50% finished"))
      FLAG<-FLAG+1
    }
    if (i/dim(dat)[1]>=0.6 & FLAG==5){
      print(paste0(ID$ID[sample],": 60% finished"))
      FLAG<-FLAG+1
    }
    if (i/dim(dat)[1]>=0.7 & FLAG==6){
      print(paste0(ID$ID[sample],": 70% finished"))
      FLAG<-FLAG+1
    }
    if (i/dim(dat)[1]>=0.8 & FLAG==7){
      print(paste0(ID$ID[sample],": 80% finished"))
      FLAG<-FLAG+1
    }
    if (i/dim(dat)[1]>=0.9 & FLAG==8){
      print(paste0(ID$ID[sample],": 90% finished"))
      FLAG<-FLAG+1
    }
    #get aligned sequences (gap, space as "-" included)
    Seq<-dat[i,4:(ID$Len[sample]+3)]
    Seq[which(Seq=="N")]<-"-"
    #read length
    dat$Read_len[i]<-nchar(dat$Seq[i])
    #read aligned length (including NTA, mismatch, deletion, but not insertion)
    dat$Align_len[i]<-sum(Seq!="-")
    #number of matched nucleotide to template
    Seq_temp_match<-(Seq==dat[1,4:(ID$Len[sample]+3)])
    dat$Match[i]<-sum(Seq_temp_match)
    #start positon of the alignment, the first nucleotide that matches the template
    dat$Start[i]<-min(which(Seq_temp_match))
    #end positon of the alignment, the last nucleotide that matches the template
    dat$End[i]<-max(which(Seq_temp_match))
    #get 5'end NTA
    SS5<-Seq[1:dat$Start[i]-1]
    #remove sapce "-"
    SS5<-SS5[SS5!="-"]
    #get SS5
    if (length(SS5)>0){
      dat$SS5[i]<-paste0(SS5,collapse = "")
    }
    #get 3'end NTA'
    SS3<-Seq[(dat$End[i]+1):length(Seq)]
    #remove space
    SS3<-SS3[SS3!="-"]
    if (length(SS3)>0){
      dat$SS3[i]<-paste0(SS3,collapse = "")
    }
    #calculate err rate (mismatch rate exclude NTA), revised to comparing the matched part of template and reads
    Mismatch<-sum(dat[1,(dat$Start[i]+3):(dat$End[i]+3)]!=dat[i,(dat$Start[i]+3):(dat$End[i]+3)])
    Hyphen<-sum(dat[i,(dat$Start[i]+3):(dat$End[i]+3)]=="-")
    # insertion/deletion is not counted as error
    Error<-Mismatch-Hyphen
    dat$Err_rate[i]<-100*Error/(dat$Align_len[i]-nchar(dat$SS3[i])-nchar(dat$SS5[i]))
    #5'end of read matches 5' end of the aligned results?
    dat$Match5[i]<-substring(dat$Seq[i],1,1)==Seq[Seq!="-"][1]
    dat$Match3[i]<-substring(dat$Seq[i],dat$Read_len[i],dat$Read_len[i])==tail(Seq[Seq!="-"],n=1)
    #calculate insertion, only meaningful when match5 and Match3 are true
    dat$Insertion[i]<-dat$Read_len[i]-dat$Align_len[i]
    #check if the read should be used for analysis
    # err rate<=10%, 5' and 3' NTA <=5 nt, insertion<=10 nt, and match5 and match3 both true
    dat$Good[i]<-((dat$Err_rate[i]<=10) & (nchar(dat$SS3[i])<=5) & 
      (nchar(dat$SS5[i])<=5) & (dat$Match3[i]=="TRUE") & 
        (dat$Match5[i]=="TRUE") & (dat$Insertion[i]<=10))
  }
  write.table(dat,gzfile(paste0(ID$ID[sample],".results.gz")),
              quote=F,sep="\t",row.names=F)
  print(paste0(ID$ID[sample],": 100% finished"))
}
