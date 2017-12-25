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
#    Author: Zhiyu Xu
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
outputF1 =      Args[2];#the output;
cores	=	as.numeric(Args[3]);
options(scipen=1000);
options(warn=-1);
library(parallel);

source(paste(Args[4],"function.R",sep = "/"));
########################################        sub functions
########        pre_stored ORFs:
mx      =       read.table(inputF1,head=F,sep="\t",stringsAsFactors=F);

denoise  =               function(i){
        tmp             =       mx[i,];
	TRID            =       as.character(tmp[1]);
        Psite           =       as.numeric(unlist(strsplit(as.character(tmp[2]),split=",",fixed=T)));
        PF_Psite	=	my_DWPT(Psite);
	out		=	data.frame(TRID,paste(PF_Psite,sep=',',collapse=','));	
	colnames(out)	=	c('TRID','PF_Psite');
	out;
}

output	 =		mclapply(X=c(1:dim(mx)[1]),FUN=denoise,mc.cores=cores,mc.preschedule=T);

all	 =		do.call(rbind.data.frame,output);

write.table(all,outputF1,sep='\t',col.names=F,row.names=F,quote=F,append=F);
