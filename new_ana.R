setwd("/stor/work/Lambowitz/yaojun/Work/RTCas1/Documents/Script/")
library(tidyr)
library(ggplot2)
library(ggseqlogo)
library(cowplot)
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
    if(xx[bbb]=="U") y[bbb]<-"A"
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

template<-list(T1="NNNNNGCAAUAAUCUAUACAAUACAACACAUACAAACAAAUUCUUAAGGUAAAAACddNNNNN",
               T2="NNNNNGCAAUAAUCUAUACAAUACAACACAUACAAACAAAUUCUUAAGGUCCCAACddNNNNN",
               T3="NNNNNGCAAUAAUCUAUACAAUACAACACAUACAAACAAAUUCUUAAGGUGGGAACddNNNNN",
               T4="NNNNNGCAAUAAUCUAUACAAUACAACACAUACAAACAAAUUCUUAAGGUUUUAACddNNNNN",
               T5="NNNNNGCAAUAAUCUAUACAAUACAACACAUACAAACAAAUUCUUAAGGUCGCAACddNNNNN",
               T6="NNNNNAACACAUACAAACAAAUUCUUAAGGUCCCAAAAAACddNNNNN",
               T7="NNNNNCAACACAUACAAACAAAUUCUUAAGGUCCCAAAAAACddNNNNN",
               T8="NNNNNACAACACAUACAAACAAAUUCUUAAGGUCCCAAAAAACddNNNNN",
               T9="NNNNNTACAACACAUACAAACAAAUUCUUAAGGUCCCAAAAAACddNNNNN",
               T10="NNNNNACAACACAUACAAACAAAUUCUUAAGGUCCCAAAAACddNNNNN",
               T11="NNNNNTACAACACAUACAAACAAAUUCUUAAGGUCCCAAAACddNNNNN")

##JA23020
meta<-read.delim("../../JA23020/meta.info")
name<-meta$ID
meta$Temp<-paste0("T",meta$Temp)
file_name<-paste0("../../JA23020/",name,".results.gz")
cs1 = make_col_scheme(chars=c('A', 'C', 'G', 'T','N','X'), groups=c('A', 'C', 'G', 'T','',''), 
                      cols=c('darkgreen', 'blue', 'orange', 'red','white','white'))

theme_logo1 <- function(base_size=12, base_family=''){
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) %+replace% 
    theme(panel.grid = element_blank(), legend.position = 'bottom', 
          axis.text.x=element_blank(),
          axis.text.y=element_blank())
}
ggseqlogo1 <- function(data, title,facet='wrap', scales='free_x', 
                       ncol=NULL, nrow=NULL, ...){
  
  # Generate the plot with default theme
  p = ggplot() + geom_logo(data = data, ...) + theme_logo1() +
    ggtitle(title)
  
  # If it's an inidivdual sequence logo, return plot
  if(!'list' %in% class(data)) return(p)
  
  # If we have more than one plot, facet
  facet_opts = c('grid', 'wrap')
  pind = pmatch(facet, facet_opts)
  facet = facet_opts[pind]
  if(is.na(facet)) stop("facet option must be set to 'wrap' or 'grid'")
  
  if(facet == 'grid'){
    p = p + facet_grid(~seq_group, scales = scales)
  }else if(facet == 'wrap'){
    p = p + facet_wrap(~seq_group, scales = scales, nrow = nrow, ncol = ncol)
  }
  
  # Return plot
  return(p)
}


for (i in 1:length(name)){
  dat<-read.delim(gzfile(file_name[i]))
  edge=dim(dat)[2]-12
  dat[,4:edge]<-sapply(dat[,4:edge],FUN=function(i){substr(as.character(i),1,1)})
  dat$Align<-apply(dat[,4:edge],MARGIN =1, FUN=function(x){paste0(x,collapse="")})
  dat<-dat[dat$Read_len>0 & dat$Good,]
  plot_title<-paste(paste(meta[i,c(2:5)],collapse =":"),sum(dat$Reads),sep=":")
  #Total_reads<-sum(dat$Reads)
  #check the strat-end distribution
  #dat$SE<-paste(dat$Start,dat$End,sep="-")
  #FL_rate<-(max(table(rep(dat$End,dat$Reads)))-dat$Reads[1])/dat$Reads[1]
  #Max_end<-as.integer(names(table(rep(dat$End,dat$Reads)))[table(rep(dat$End,dat$Reads))==
  #                                                           max(table(rep(dat$End,dat$Reads)))])
  #if (FL_rate<0.9 | (dat$End[1]-Max_end)>2){
  #  print(paste0("Warning:",name[i],"-Low % of reads reach to the end",":",round(FL_rate,2),":",Max_end))
  #} else {
  #  print(paste0(name[i],"-Pass the check"))
  #}
  # Only need to look at the start position
  
  # get the top10 start position
  top10<-data.frame(sort(table(rep(dat$Start,dat$Reads)),decreasing = T))
  top10$Var1<-as.integer(as.character(top10$Var1))
  top10$Cum<-cumsum(top10$Freq)/sum(top10$Freq)
  top10$Freq<-top10$Freq/sum(top10$Freq)
  top10<-top10[order(top10$Var1),]
  #pdf(paste0(name[i],"_start_pos.pdf"))
  #plot(top10$Freq~top10$Var1,type="l",bty="n",xlim=rev(range(top10$Var1)))
  #dev.off()
  top10<-top10[top10$Freq>=0.01,]
  print(paste(name[i],sum(top10$Freq>=0.01),tail(top10$Cum[top10$Freq>=0.01],1)),sep="-")
  #number of positons and the corresponds cumulative % of reads for positions >=1% reads
  # dataset|number of start position to show|cumulative % of reads
  # "1 15 0.906704572355091"
  # "2 4 0.950312524473178"
  # "3 9 0.916294753306187"
  # "4 13 0.913184305140112"
  # "5 7 0.967051070840198"
  # "6 17 0.801579405059424"
  # "7 7 0.944310182632108"
  # "8 10 0.919773770676621"
  # "9 16 0.895699228307924"
  # "10 5 0.946229936705384"
  # "11 7 0.96822049684144"
  # "12 7 0.972872380359807"
  # "13 4 0.948952979808744"
  # "14 5 0.969754350368311"
  # "15 5 0.979612386502393"
  # "16 8 0.879002996034353"
  # "17 7 0.902180564998907"
  # "18 7 0.944660927021912"
  # "19 9 0.93851032475047"
  # "20 11 0.9424260893981"
  # "21 8 0.975270550785254"
  # "22 6 0.899534236759468"
  # "23 8 0.765010911318662"
  # "24 12 0.699949045794332"
  # "25 12 0.843990767745306"
  # "26 7 0.939074049484354"
  # "27 10 0.96925362898617"
  # "28 10 0.890733548492971"
  # "29 12 0.920589703069117"
  # "30 11 0.959040506902014"
  # "31 14 0.95170861469204"
  # "32 11 0.82781677133358"
  positions<-top10$Var1
  for (j in 1:length(positions)){
    Pos<-positions[j]
    tmp<-dat[dat$Start==Pos,]
    #create frequency matrix
    for (k in edge:4){
      if (k==edge){
        frq<-data.frame(table(rep(tmp[,k],tmp$Reads)))
      } else {
        frq<-merge(frq,data.frame(table(rep(tmp[,k],tmp$Reads))),by=1,all=T)
      }
    }
    frq[is.na(frq)]<-0
    rownames(frq)<-frq$Var1
    frq<-prop.table(as.matrix(frq[,-1]),2)
    rownames(frq)[1]<-"X"
    colnames(frq)<-paste0("Pos",1:dim(frq)[2])
    if (j==1){
      logo_list<-list(frq)
      names(logo_list)[j]<-paste(name[i],paste0("Pos",top10$Var1[j]),
                                 paste0(round(100*top10$Freq[j],1),"%"),sep=":")
    } else {
      logo_list<-c(logo_list,list(frq))
      names(logo_list)[j]<-paste(name[i],paste0("Pos",top10$Var1[j]),
                                 paste0(round(100*top10$Freq[j],1),"%"),sep=":")
    }
    
  }
  gglogo<-ggseqlogo1(logo_list,plot_title,ncol=1,method="p",
                     col_scheme=cs1,namespace="ACGTNX")
  height_N<-length(logo_list)
  ggsave2(plot=gglogo, file=paste0("../../JA23020/",name[i],"_seqlogo.pdf"), width=11, height=0.75*height_N)
}  


cs2 = make_col_scheme(chars=c('A', 'C', 'G', 'T','U','d','N'), 
                      groups=c('A', 'C', 'G', 'T','U','d','N'), 
                      cols=c('darkgreen', 'blue', 'orange', 'red','red','black','white'))

gglogo<-ggseqlogo1(template,"template",ncol=1,method="p",col_scheme=cs2,namespace="ACGTUdN")
ggsave2(plot=gglogo, file="../../JA23020/template_seqlogo.pdf", width=11, height=0.75*length(template))


###-------full length snap-back ana-------------
### dataset 5 11 16 using CGC
for (i in 17:20){
  dat<-read.delim(paste0("../../JA23020/snap_back/",i,".output"),header=F)
  dat<-separate(dat,"V1",into=c("Reads","Seq"),sep="@")
  dat<-separate(dat,"Reads",into=c("ID","Reads"),sep="_",remove=F)
  dat$Seq[1]<-dat$V2[1]
  dat$Reads[1]<-0
  dat$ID[-1]<-paste(dat$ID[-1],1:(dim(dat)[1]-1),sep="_")
  if (i==17){
    sheet<-list(dat)
  } else {
    sheet<-c(sheet,list(dat))
  }
}
names(sheet)<-paste0("dataset_",17:20)
write.xlsx(sheet, file = "../../JA23020/snap_back.xlsx")

### 1-10 cDNA start plot
### 1% product cutoff
### 1$ break point is 24
### 5$ break point is 35
cutoff<-1
break_point<-24
bcol<-c("green","blue","darkorange","red")
pdf("Fig5AB.pdf",width=12,height=8)
par(pch=16,cex=0.7,mfcol=c(5,2),xpd=FALSE)
par(mar=c(1,4,2,1))
for (sample in 1:15){
  dat<-read.delim(gzfile(file_name[sample]))
  edge=dim(dat)[2]-12
  dat[,4:edge]<-sapply(dat[,4:edge],FUN=function(i){substr(as.character(i),1,1)})
  dat$Align<-apply(dat[,4:edge],MARGIN =1, FUN=function(x){paste0(x,collapse="")})
  dat<-dat[dat$Read_len>0 & dat$Good,]
  dat$PureAlign<-gsub("-","",dat$Align)
  plot_title<-paste0(meta$Template[sample],", Mn",meta$Mn[sample],", ",meta$Primer[sample]," primer")
  if (dim(dat)[1]==0){next}
  temp<-template[meta$Temp[sample]]
  #cDNA start
  tmp<-aggregate(dat$Reads~dat$Start,dat,sum)
  colnames(tmp)<-c("Len","Freq")
  tmp$Len<-56-tmp$Len
  tmp$Nuc3<-sapply(tmp$Len,FUN=function(x){substr(temp,x+3,x+3)})
  tmp$Nuc2<-sapply(tmp$Len,FUN=function(x){substr(temp,x+4,x+4)})
  tmp$Nuc1<-sapply(tmp$Len,FUN=function(x){substr(temp,x+5,x+5)})
  Nuc3<-aggregate(tmp$Freq~tmp$Nuc3,dat,sum)
  Nuc2<-aggregate(tmp$Freq~tmp$Nuc2,dat,sum)
  Nuc1<-aggregate(tmp$Freq~tmp$Nuc1,dat,sum)
  Len<-data.frame("Len"=c(1:(meta$Length[sample]-10)))
  tmp<-merge(Len,tmp[,1:2],by=1,all=T)
  tmp[is.na(tmp)]<-0
  tmp$Freq<-100*tmp$Freq/sum(tmp$Freq)
  y_lim<-ceiling(max(tmp$Freq)/10)*10
  tmp$Freq[tmp$Freq<cutoff]<-NA
  tmp<-tmp[(break_point-5):dim(tmp)[1],]
  plot(tmp$Freq~tmp$Len,type="h",xlim=range(tmp$Len)+c(0,9),ylim=c(-y_lim*0.02,y_lim),
       lwd=2,bty="n",xaxt="n",xlab=NA,yaxt="n",
       ylab="Reads (%)",main=plot_title)
  axis(2,at=seq(0,y_lim,10),labels =seq(0,y_lim,10),las=2)
  par(xpd=TRUE)
  text(18,-y_lim*0.06,"5'-",cex=1)
  for (j in 1:3){
    text(j+18,-y_lim*0.06,substr(temp,j+5,j+5),col="red",cex=1)
  }
  segments(21.5,-y_lim*0.06,22,-y_lim*0.06)
  segments(23,-y_lim*0.06,23.5,-y_lim*0.06)
  segments(22.4,-y_lim*0.08,23.2,-y_lim*0.04)
  segments(21.9,-y_lim*0.08,22.7,-y_lim*0.04)
  for (j in 24:50){
    text(j,-y_lim*0.06,substr(temp,j+5,j+5),col="red",cex=1)
  }
  text(52,-y_lim*0.06,"ddC-3'",cex=1)
  ### add 1,2,3 bars
  Nuc<-Nuc3
  colnames(Nuc)<-c("Nuc","Freq")
  Nuc<-Nuc[Nuc$Nuc!="N",]
  Nuc$Nuc<-c("T","G","C","A")
  Nuc<-Nuc[order(Nuc$Nuc),]
  Nuc$Freq<-100*Nuc$Freq/sum(Nuc$Freq)
  Nuc$Freq<-Nuc$Freq*y_lim/100
  Nuc$Freq<-cumsum(Nuc$Freq)
  for (j in 1:4){
    if (j==1){
      y_base=0
    } else {
      y_base=Nuc$Freq[j-1]
    }
    rect(54,y_base,55.5,sum(Nuc$Freq[j]),
         col=bcol[j])
 
  }
  ## pos2
  Nuc<-Nuc2
  colnames(Nuc)<-c("Nuc","Freq")
  Nuc<-Nuc[Nuc$Nuc!="N",]
  Nuc$Nuc<-c("T","G","C","A")
  Nuc<-Nuc[order(Nuc$Nuc),]
  Nuc$Freq<-100*Nuc$Freq/sum(Nuc$Freq)
  Nuc$Freq<-Nuc$Freq*y_lim/100
  Nuc$Freq<-cumsum(Nuc$Freq)
  for (j in 1:4){
    if (j==1){
      y_base=0
    } else {
      y_base=Nuc$Freq[j-1]
    }
    rect(56,y_base,57.5,sum(Nuc$Freq[j]),
         col=bcol[j])
    
  }
  ## pos1
  Nuc<-Nuc1
  colnames(Nuc)<-c("Nuc","Freq")
  Nuc<-Nuc[Nuc$Nuc!="N",]
  Nuc$Nuc<-c("T","G","C","A")
  Nuc<-Nuc[order(Nuc$Nuc),]
  Nuc$Freq<-100*Nuc$Freq/sum(Nuc$Freq)
  Nuc$Freq<-Nuc$Freq*y_lim/100
  Nuc$Freq<-cumsum(Nuc$Freq)
  for (j in 1:4){
    if (j==1){
      y_base=0
    } else {
      y_base=Nuc$Freq[j-1]
    }
    rect(58,y_base,59.5,sum(Nuc$Freq[j]),
         col=bcol[j])
    
  }
  text(seq(54.75,58.75,2),-y_lim*0.06,labels = paste0("N",3:1),cex=1)
  #axis(4,las=2,at=seq(0,y_lim,y_lim/4),labels = seq(0,100,25))
  #text(66,y_lim/2,"Nucleotide (%)",srt=90)
}
dev.off()


### 21-32
cutoff<-1
bcol<-c("green","blue","darkorange","red")
pdf("Fig5C1.pdf",width=11,height=9.6)
par(pch=16,cex=0.7,mfcol=c(6,2),xpd=FALSE)
par(mar=c(1,4,2,1))
for (sample in 21:32){
  break_point<-meta$Length[sample]-36
  dat<-read.delim(gzfile(file_name[sample]))
  edge=dim(dat)[2]-12
  dat[,4:edge]<-sapply(dat[,4:edge],FUN=function(i){substr(as.character(i),1,1)})
  dat$Align<-apply(dat[,4:edge],MARGIN =1, FUN=function(x){paste0(x,collapse="")})
  dat<-dat[dat$Read_len>0 & dat$Good,]
  dat$PureAlign<-gsub("-","",dat$Align)
  plot_title<-paste0(meta$Template[sample],", Mn",meta$Mn[sample],", ",meta$Primer[sample]," primer")
  if (dim(dat)[1]==0){next}
  temp<-template[meta$Temp[sample]]
  #cDNA start
  tmp<-aggregate(dat$Reads~dat$Start,dat,sum)
  colnames(tmp)<-c("Len","Freq")
  tmp$Len<-meta$Length[sample]-4-tmp$Len
  tmp$Nuc3<-sapply(tmp$Len,FUN=function(x){substr(temp,x+3,x+3)})
  tmp$Nuc2<-sapply(tmp$Len,FUN=function(x){substr(temp,x+4,x+4)})
  tmp$Nuc1<-sapply(tmp$Len,FUN=function(x){substr(temp,x+5,x+5)})
  Nuc3<-aggregate(tmp$Freq~tmp$Nuc3,dat,sum)
  Nuc2<-aggregate(tmp$Freq~tmp$Nuc2,dat,sum)
  Nuc1<-aggregate(tmp$Freq~tmp$Nuc1,dat,sum)
  Len<-data.frame("Len"=c(1:(meta$Length[sample]-10)))
  tmp<-merge(Len,tmp[,1:2],by=1,all=T)
  tmp[is.na(tmp)]<-0
  tmp$Freq<-100*tmp$Freq/sum(tmp$Freq)
  y_lim<-ceiling(max(tmp$Freq)/10)*10
  tmp$Freq[tmp$Freq<cutoff]<-NA
  tmp<-tmp[(break_point-5):dim(tmp)[1],]
  plot(tmp$Freq~tmp$Len,type="h",xlim=range(tmp$Len)+c(0,9),
       ylim=c(-y_lim*0.02,y_lim),
       lwd=2,bty="n",xaxt="n",xlab=NA,yaxt="n",
       ylab="Reads (%)",main=plot_title)
  axis(2,at=seq(0,y_lim,10),labels =seq(0,y_lim,10),las=2)
  par(xpd=TRUE)
  text(break_point-6,-y_lim*0.06,"5'-",cex=1)
  for (j in 1:3){
    text(j+break_point-6,-y_lim*0.06,substr(temp,j+5,j+5),col="red",cex=1)
  }
  segments(break_point-2.5,-y_lim*0.06,break_point-2,-y_lim*0.06)
  segments(break_point-1,-y_lim*0.06,break_point-0.5,-y_lim*0.06)
  segments(break_point-1.6,-y_lim*0.08,break_point-0.8,-y_lim*0.04)
  segments(break_point-2.1,-y_lim*0.08,break_point-1.3,-y_lim*0.04)
  for (j in (break_point):(meta$Length[sample]-10)){
    text(j,-y_lim*0.06,substr(temp,j+5,j+5),col="red",cex=1)
  }
  text(meta$Length[sample]-8,-y_lim*0.06,"ddC-3'",cex=1)
  ### add 1,2,3 bars
  Nuc<-Nuc3
  colnames(Nuc)<-c("Nuc","Freq")
  Nuc<-Nuc[Nuc$Nuc!="N",]
  Nuc$Nuc<-c("T","G","C","A")
  Nuc<-Nuc[order(Nuc$Nuc),]
  Nuc$Freq<-100*Nuc$Freq/sum(Nuc$Freq)
  Nuc$Freq<-Nuc$Freq*y_lim/100
  Nuc$Freq<-cumsum(Nuc$Freq)
  for (j in 1:4){
    if (j==1){
      y_base=0
    } else {
      y_base=Nuc$Freq[j-1]
    }
    rect(meta$Length[sample]-6,y_base,meta$Length[sample]-4.5,sum(Nuc$Freq[j]),
         col=bcol[j])
    
  }
  ## pos2
  Nuc<-Nuc2
  colnames(Nuc)<-c("Nuc","Freq")
  Nuc<-Nuc[Nuc$Nuc!="N",]
  Nuc$Nuc<-c("T","G","C","A")
  Nuc<-Nuc[order(Nuc$Nuc),]
  Nuc$Freq<-100*Nuc$Freq/sum(Nuc$Freq)
  Nuc$Freq<-Nuc$Freq*y_lim/100
  Nuc$Freq<-cumsum(Nuc$Freq)
  for (j in 1:4){
    if (j==1){
      y_base=0
    } else {
      y_base=Nuc$Freq[j-1]
    }
    rect(meta$Length[sample]-4,y_base,meta$Length[sample]-2.5,sum(Nuc$Freq[j]),
         col=bcol[j])
    
  }
  ## pos1
  Nuc<-Nuc1
  colnames(Nuc)<-c("Nuc","Freq")
  Nuc<-Nuc[Nuc$Nuc!="N",]
  Nuc$Nuc<-c("T","G","C","A")
  Nuc<-Nuc[order(Nuc$Nuc),]
  Nuc$Freq<-100*Nuc$Freq/sum(Nuc$Freq)
  Nuc$Freq<-Nuc$Freq*y_lim/100
  Nuc$Freq<-cumsum(Nuc$Freq)
  for (j in 1:4){
    if (j==1){
      y_base=0
    } else {
      y_base=Nuc$Freq[j-1]
    }
    rect(meta$Length[sample]-2,y_base,meta$Length[sample]-0.5,sum(Nuc$Freq[j]),
         col=bcol[j])
    
  }
  text(seq(meta$Length[sample]-5.25,meta$Length[sample]-1.25,2),-y_lim*0.06,labels = paste0("N",3:1),cex=1)
  #axis(4,las=2,at=seq(0,y_lim,y_lim/4),labels = seq(0,100,25))
  #text(66,y_lim/2,"Nucleotide (%)",srt=90)
}
dev.off()


###JA22482

### ---make -bar-----
meta<-read.delim("../../JA22482/meta.info")
meta<-meta[c(12:14,18:20),]
meta$Temp<-rep(1:3,2)
meta$xMax<-rep(c(8,9,11),2)
file_name<-paste0("../../JA22482/",meta$ID,".results.gz")
template<-list(R29_3dA="UUUCUCGAGUCAUCUUUUAGGGCUCCAAGAAAAA",
               R29_4dA="UUUCUCGAGUCAUCUUUUAGGGCUCCAAGAAAAAA",
               R29_6dA="UUUCUCGAGUCAUCUUUUAGGGCUCCAAGAAAAAAAA")
### T14-16 T18-20
cutoff<-1
bcol<-c("green","blue","darkorange","red")
pdf("Fig5C2.pdf",width=11,height=9.6)
par(pch=16,cex=0.7,mfcol=c(6,2),xpd=FALSE)
par(mar=c(1,4,2,1))
for (sample in 1:6){
  break_point<-nchar(template[[meta$Temp[sample]]])-25
  dat<-read.delim(gzfile(file_name[sample]))
  edge=dim(dat)[2]-12
  dat[,4:edge]<-sapply(dat[,4:edge],FUN=function(i){substr(as.character(i),1,1)})
  dat$Align<-apply(dat[,4:edge],MARGIN =1, FUN=function(x){paste0(x,collapse="")})
  dat<-dat[dat$Read_len>0 & dat$Good,]
  dat$PureAlign<-gsub("-","",dat$Align)
  plot_title<-paste0(meta$Template[sample],", Mn",meta$Mn[sample],", ",meta$Primer[sample]," primer")
  if (dim(dat)[1]==0){next}
  temp_for_x<-template[[meta$Temp[sample]]]
  temp<-template[[2]]
  #cDNA start
  tmp<-aggregate(dat$Reads~dat$Start,dat,sum)
  colnames(tmp)<-c("Len","Freq")
  tmp$Len<-meta$Length[sample]-4-tmp$Len
  tmp$Nuc3<-sapply(tmp$Len,FUN=function(x){substr(temp,x-2,x-2)})
  tmp$Nuc2<-sapply(tmp$Len,FUN=function(x){substr(temp,x-1,x-1)})
  tmp$Nuc1<-sapply(tmp$Len,FUN=function(x){substr(temp,x,x)})
  Nuc3<-aggregate(tmp$Freq~tmp$Nuc3,dat,sum)
  Nuc2<-aggregate(tmp$Freq~tmp$Nuc2,dat,sum)
  Nuc1<-aggregate(tmp$Freq~tmp$Nuc1,dat,sum)
  
  Len<-data.frame("Len"=c(1:(meta$Length[sample]-10)))
  tmp<-merge(Len,tmp[,1:2],by=1,all=T)
  tmp[is.na(tmp)]<-0
  tmp$Freq<-100*tmp$Freq/sum(tmp$Freq)
  y_lim<-ceiling(max(tmp$Freq)/10)*10
  tmp$Freq[tmp$Freq<cutoff]<-NA
  tmp<-tmp[(break_point-5):dim(tmp)[1],]
  plot(tmp$Freq~tmp$Len,type="h",xlim=range(tmp$Len)+c(0,meta$xMax[sample]),
       ylim=c(-y_lim*0.02,y_lim),
       lwd=2,bty="n",xaxt="n",xlab=NA,yaxt="n",
       ylab="Reads (%)",main=plot_title)
  axis(2,at=seq(0,y_lim,10),labels =seq(0,y_lim,10),las=2)
  par(xpd=TRUE)
  text(break_point-6,-y_lim*0.06,"5'-",cex=1)
  for (j in 1:3){
    text(j+break_point-6,-y_lim*0.06,substr(temp,j,j),col="red",cex=1)
  }
  segments(break_point-2.5,-y_lim*0.06,break_point-2,-y_lim*0.06)
  segments(break_point-1,-y_lim*0.06,break_point-0.5,-y_lim*0.06)
  segments(break_point-1.6,-y_lim*0.08,break_point-0.8,-y_lim*0.04)
  segments(break_point-2.1,-y_lim*0.08,break_point-1.3,-y_lim*0.04)
  for (j in (break_point):(nchar(temp_for_x))){
    text(j,-y_lim*0.06,substr(temp_for_x,j,j),col="red",cex=1)
  }
  text(meta$xMax[sample]+28,-y_lim*0.06,"ddC-3'",cex=1)
  ### add 1,2,3 bars
  Nuc<-merge(c("A","C","G","U"),Nuc3,by=1,all=T)
  Nuc[is.na(Nuc)]<-0
  colnames(Nuc)<-c("Nuc","Freq")
  Nuc<-Nuc[Nuc$Nuc!="N",]
  Nuc<-Nuc[Nuc$Nuc!="",]
  Nuc$Nuc<-c("T","G","C","A")
  Nuc<-Nuc[order(Nuc$Nuc),]
  Nuc$Freq<-100*Nuc$Freq/sum(Nuc$Freq)
  Nuc$Freq<-Nuc$Freq*y_lim/100
  Nuc$Freq<-cumsum(Nuc$Freq)
  for (j in 1:4){
    if (j==1){
      y_base=0
    } else {
      y_base=Nuc$Freq[j-1]
    }
    rect(meta$xMax[sample]+37-6,y_base,meta$xMax[sample]+37-4.5,sum(Nuc$Freq[j]),
         col=bcol[j])
    
  }
  ## pos2
  Nuc<-merge(c("A","C","G","U"),Nuc2,by=1,all=T)
  Nuc[is.na(Nuc)]<-0
  colnames(Nuc)<-c("Nuc","Freq")
  Nuc<-Nuc[Nuc$Nuc!="N",]
  Nuc$Nuc<-c("T","G","C","A")
  Nuc<-Nuc[order(Nuc$Nuc),]
  Nuc$Freq<-100*Nuc$Freq/sum(Nuc$Freq)
  Nuc$Freq<-Nuc$Freq*y_lim/100
  Nuc$Freq<-cumsum(Nuc$Freq)
  for (j in 1:4){
    if (j==1){
      y_base=0
    } else {
      y_base=Nuc$Freq[j-1]
    }
    rect(meta$xMax[sample]+37-4,y_base,meta$xMax[sample]+37-2.5,sum(Nuc$Freq[j]),
         col=bcol[j])
    
  }
  ## pos1
  Nuc<-merge(c("A","C","G","U"),Nuc1,by=1,all=T)
  Nuc[is.na(Nuc)]<-0
  colnames(Nuc)<-c("Nuc","Freq")
  Nuc<-Nuc[Nuc$Nuc!="N",]
  Nuc$Nuc<-c("T","G","C","A")
  Nuc<-Nuc[order(Nuc$Nuc),]
  Nuc$Freq<-100*Nuc$Freq/sum(Nuc$Freq)
  Nuc$Freq<-Nuc$Freq*y_lim/100
  Nuc$Freq<-cumsum(Nuc$Freq)
  for (j in 1:4){
    if (j==1){
      y_base=0
    } else {
      y_base=Nuc$Freq[j-1]
    }
    rect(meta$xMax[sample]+37-2,y_base,meta$xMax[sample]+37-0.5,sum(Nuc$Freq[j]),
         col=bcol[j])
    
  }
  text(seq(meta$xMax[sample]+37-5.25,meta$xMax[sample]+37-1.25,2),-y_lim*0.06,labels = paste0("N",3:1),cex=1)
  #axis(4,las=2,at=seq(0,y_lim,y_lim/4),labels = seq(0,100,25))
  #text(66,y_lim/2,"Nucleotide (%)",srt=90)
}
dev.off()

