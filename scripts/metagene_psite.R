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
inputF1= Args[1];#
output1= Args[2];#
output2= Args[3];#
options(warn=-1);
options(scipen=1000);
library(reshape);
library(ggplot2);

####	assign function:
my_fuc=function(X,STR,mark){
	x0	=	eval(parse(text=X))
	if(STR=="-"){
		x0	=	rev(x0);
	}
	y	=	which(x0==mark);
	out	=	c(length(x0),y);
}
####	
TR_all	=	read.table(inputF1,sep="\t",head=F);
out	=	matrix(0,nrow=nrow(TR_all),ncol=2);
for(i in 1:nrow(TR_all)){
	X	=	as.vector(TR_all[i,3]);
	STR	=	as.vector(TR_all[i,2]);
	mark	=	as.numeric(as.vector(TR_all[i,4]));
	out[i,]	=	my_fuc(X,STR,mark);
}

L	=	sort(unique(out[,1]));
L	=	L[L<37];
out2	=	matrix(0,nrow=length(L),ncol=max(L));
for(i in 1:length(L)){
	tmp	=	out[which(out[,1]==L[i]),2];
	for(j in 1:max(L)){
		out2[i,j]	=	length(which(tmp==j));
	}
}
psite_coor	=	1:max(L);
colnames(out2)	=	psite_coor;
L		=	paste(L, " nt reads",sep="");
rownames(out2)	=	L;
write.table(out2,file=output1,col.name=F,row.name=T,quote=F,sep="\t");



out4		=	melt(out2,id=1)
colnames(out4)	=	c("readL","psite_coor","counts");
#out4$psite_coor	=	factor(out4$psite_coor,levels=psite_coor)

hah=ggplot(out4,aes(x=psite_coor,y=counts))+
geom_bar(position=position_dodge(),stat="identity",colour="white",fill="grey")+
facet_grid(readL ~. ,scales="free_y")+
xlab("coordinates")+
ylab("counts")+
theme_bw() +
theme(axis.title.x = element_text(size=0),axis.text.x  = element_text(size=15),axis.title.y = element_text(size=0),axis.text.y  = element_text(size=15))+
pdf(output2,width=3.5,height=12)
print(hah);
dev.off()


print("ALLDONE");
rm(list=ls());
