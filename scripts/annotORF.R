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
#    Author: Long Hu, Zhiyu Xu
#
#    Copyright (C) 2017  Zhiyu Xu
#
#    RiboWave is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#    Contact: xanthexu18@gmail.com
#######################################################################

Args <-commandArgs(TRUE);
inputF1= Args[1];#input
output1= Args[2];#output 
options(scipen=1000);

mx		=	read.table(inputF1,head=F,sep="\t")
mx2		=	as.matrix(mx[,1:4]);
for(i in 1:nrow(mx)){
	TR	=	as.vector(mx[i,1]);
	STR	=	as.vector(mx[i,2]);
	exon	=	sort(unique(eval(parse(text=as.character(mx[i,3])))));
	Start	=	sort(unique(eval(parse(text=as.character(mx[i,4])))));
	Stop	=	sort(unique(eval(parse(text=as.character(mx[i,5])))));
	if(STR=="-"){
	exon	=	rev(exon);
	L	=	which(exon==max(Start));
	R	=	which(exon==max(Stop))-1;
	}else{
	L	=	which(exon==min(Start));
	R	=	which(exon==min(Stop))-1;
	}
	tmp	=	c(TR,STR,L,R)
	mx2[i,]	=	tmp;
}

write.table(mx2,file=output1,quote=F,sep="\t",col.name=F,row.name=F);
rm(list=ls());

