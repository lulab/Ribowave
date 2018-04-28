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
inputF1= Args[1];#
output1= Args[2];#
options(scipen=1000);

mx	=	read.table(inputF1,sep="\t",head=F);
mx2	=	mx[,-1];
SUM1	=	rowSums(mx2);
mx3	=	mx2/SUM1;
SUM2	=	SUM1/sum(SUM1)
out	=	cbind(mx[,1],SUM1,SUM2,mx3)


write.table(out,file=output1,col.name=F,row.name=F,quote=F,sep="\t");

rm(list=ls());
