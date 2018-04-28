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
inputF1 =       Args[1];#the transcript psite
inputF2 =       Args[2];#the ORF file;
outputF1=       Args[3];#the output;
options(scipen=1000);
options(warn=-1);
library(parallel);
cores	=	as.numeric(Args[4]);

source(paste(Args[5],"function.R",sep = "/"));
########################################        sub functions
########        pre_stored ORFs:
ORFs    =       read.table(inputF2,head=F,sep="\t",stringsAsFactors=F);
colnames(ORFs)= c("ORFID","TRtype","label","codon","getorf","orflonger","TRID","Start","Stop");

mx	=       read.table(inputF1,head=F,sep="\t",stringsAsFactors=F);

my_cov	=		function(i){
	tmp             =       mx[i,];
        TRID            =       as.character(tmp[1]);
        TR_ORF          =       ORFs[which(ORFs[,7]==TRID),c(8,9)];
        signal          =       as.numeric(unlist(strsplit(as.character(tmp[2]),split=",",fixed=T)));
	density		=	den_check(TR_ORF,signal);
	colnames(density)	=	c("trans","ORF","ORF_DEN");	
	out             =       cbind(ORFs[which(ORFs[,7]==TRID),c(1,7,8,9)],density);
}

output   =              mclapply(X=c(1:dim(mx)[1]),FUN=my_cov,mc.cores=cores,mc.preschedule=T);

all      =              do.call(rbind.data.frame,output);

write.table(all,outputF1,sep='\t',col.names=T,row.names=F,quote=F,append=F);

