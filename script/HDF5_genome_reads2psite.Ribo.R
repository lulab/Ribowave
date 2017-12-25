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
inputF1= Args[1];#input the reads coor file
inputF2= Args[2];#read_length, psite position;
inputF3= Args[3];#chromosome length;
inputF4= Args[4];#strand;
output1= Args[5];#output the chromosomes psite HDF5 file
library(methods);
library(rhdf5)
options(warn=-1);

########	sub function
my_parse	=	function(x,y){
	tmp	=	eval(parse(text=x));
	tmp[y];
}

reads2psite	=	function(reads_mx,readsL_psite){
	out	=	c();
	for(i in 1:nrow(readsL_psite)){
		readsL	=	readsL_psite[i,1];
		psite_pos=	readsL_psite[i,2];
		ID	=	which(reads_mx[,1]==readsL);
		coor_mx	=	as.vector(reads_mx[ID,2]);
		tmp	=	sapply(coor_mx,my_parse,y=psite_pos,USE.NAMES = FALSE);
		out	=	c(out,tmp);
	}
	out;
}

######## main function:
chrVecL		=	as.numeric(inputF3);
chrVec		<-	rep(0,length=chrVecL);
readsL_psite	=	read.table(inputF2);
if(inputF4=="-"){
readsL_psite[,2]=	readsL_psite[,1]-readsL_psite[,2]+1	
}

Nlines		=       10000;
INid		=       file(inputF1,"r");                      
IN		=	readLines(INid,n=Nlines);                       
while(length(IN) != 0){
	reads_mx	=	read.table(text=IN,sep="\t",head=F,quote="",check.names=F);
	tmp		=	reads2psite(reads_mx,readsL_psite);
	if(length(tmp)>0){
	counts		=	table(tmp);
	coors		=	as.numeric(names(counts));
	chrVec[coors]	<-	chrVec[coors]+counts;
	}
	IN		=	readLines(INid,n=Nlines);                       
}
close(INid);

h5createFile(output1);
h5write(chrVec,output1,"HDF5file");


rm(list=ls());

