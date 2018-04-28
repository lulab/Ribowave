#!/usr/bin/env Rscript

###################################################################
#    This file is part of RiboWave.
#    RiboWave is powerful Ribo-seq analysis tool that is able to 
#    denoise the Ribo-seq data and serve for multiple functions.
#	
#    RiboWave can be used for multiple purposes: 
#    		1. denoise the raw Ribo-seq data
#		2. define translated ORFs
#		3. estimate the abundance of actively elongating ribosomes
#		4. estimate translation efficiency(TE)
#		5. identify potential frameshift candidates
#    
#    Author: Long Hu
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
inputF1= Args[1];#
output1= Args[2];#
options(scipen=1000);

my_max_psite	=	function(X,L){
initial_psite_set	=	8:16;
out			=	X[initial_psite_set];
best_psite		=	initial_psite_set[which.max(out)];
}

mx	=	read.table(inputF1,sep="\t",head=F);
SUMall	=	sum(as.vector(mx[,-1]))
out2	=	matrix(0,nrow=nrow(mx),ncol=(ncol(mx)+2));
for(i in 1:nrow(mx)){
X	=	as.vector(as.matrix(mx[i,-1]));
L	=	as.vector(as.matrix(mx[i,1]));
X	=	X[1:L];
SUM1	=	sum(X);
SUM	=	round(sum(X)/SUMall,4);
tmp	=	my_max_psite(X,L);
out2[i,c(1,2,3,tmp+3)]=	c(L,SUM1,SUM,1)
}
write.table(out2,file=output1,col.name=F,row.name=F,quote=F,sep="\t");

rm(list=ls());

