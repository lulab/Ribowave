#!/usr/bin/env Rscript

###################################################################
#    This file is part of RiboWave.
#    RiboWave is powerful Ribo-seq analysis tool that is able to 
#    denoise the Ribo-seq data and serve for multiple functions.
#       
#    RiboWave can be used for multiple purposes: 
#               1. denoise the raw Ribo-seq data
#               2. define translated ORFs
#               3. estimate the abundance of actively elongating ribosomes
#               4. estimate translation efficiency(TE)
#               5. identify potential frameshift candidates
#    
#    Author: Zhiyu Xu, Long Hu
#
#    Copyright (C) 2017  Zhiyu Xu
#
#    RiboWave is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#    Contact: xanthexu18@gmail.com
#######################################################################


Args    =       commandArgs(TRUE);
inputF1 =       Args[1];#the denoised signal 
inputF2	=	Args[2];#the raw signal
inputF3 =       Args[3];#the ORF file;
outputF1 =      Args[4];#the output;
cores	=	as.numeric(Args[5]);
options(scipen=1000);
options(warn=-1);
library(parallel);

source(paste(Args[6],"function.R",sep = "/"));
########################################        sub functions
########        pre_stored ORFs:
ORFs    =       read.table(inputF3,head=F,sep="\t",stringsAsFactors=F);
colnames(ORFs)= c("ORFID","TRtype","label","codon","getorf","orflonger","TRID","Start","Stop");

mx      =       read.table(inputF1,head=F,sep="\t",stringsAsFactors=F);
raw      =       read.table(inputF2,head=F,sep="\t",stringsAsFactors=F);

my_feat  =               function(i){
        tmp             =       mx[i,];
        TRID            =       as.character(tmp[1]);
        denoise		=	as.numeric(unlist(strsplit(as.character(tmp[2]),split=",",fixed=T)));
	raw		=	as.numeric(unlist(strsplit(as.character(raw[which(raw[,1]==TRID),2]),split=",",fixed=T)));
	TR_ORF          =       ORFs[which(ORFs[,7]==TRID),c(8,9)];
        if(max(denoise)>0){
                feats   =       feat_calc(TR_ORF,denoise,raw);
        }else{
                feats   =       matrix(rep(c(1,0,0),nrow(TR_ORF)),ncol=3,byrow=T);
		colnames(feats)=c("pvalue","OUT_3nt","startPsite");
        }
	cov		=	cov_check(TR_ORF,raw);
	colnames(cov)	=	c('CovORF');
        out             =       cbind(ORFs[which(ORFs[,7]==TRID),c(1,7,8,9)],cov,feats);
}

output	 =		mclapply(X=c(1:dim(mx)[1]),FUN=my_feat,mc.cores=cores,mc.preschedule=T);

all	 =		do.call(rbind.data.frame,output);

write.table(all,outputF1,sep='\t',col.names=T,row.names=F,quote=F,append=F);
