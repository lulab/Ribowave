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
inputF1= Args[1];#input the CP h5 file of this chromsome.str;
inputF2= Args[2];#input coor file
inputF3= Args[3];#input strand
output1= Args[4];#output psite of each transcripts

library(methods);
library(rhdf5);
####	assign function:
my_fuc=function(X,STR,CP){
x0	=	eval(parse(text=X))
x0	=	sort(unique(x0));
tmp	=	CP[x0];
if(STR=="-"){
tmp	=	rev(tmp);
}
y	=	paste(tmp,collapse=",");
}
####	
CP	=	h5read(inputF1,"HDF5file");
TR	=	read.table(inputF2,sep="\t",head=F);
TRID	=	as.character(as.vector(TR[,1]));
coor	=	as.character(as.vector(TR[,3]));
tmp	=	sapply(coor,my_fuc,STR=inputF3,CP=CP,USE.NAMES=FALSE);
out	=	cbind(TRID,tmp);
write.table(out,file=output1,quote=F,sep="\t",col.name=F,row.name=F);
rm(list=ls());
