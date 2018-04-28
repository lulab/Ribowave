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
outputF1=       Args[2];#the output;
options(scipen=1000);
options(warn=-1);
library(parallel);

source(paste(Args[3],"function.R",sep = "/"));
source(paste(Args[3],"CRF_function.R",sep = "/"));
########################################        sub functions
########        pre_stored ORFs:
INid    =       file(inputF1,"r");
IN      =       readLines(INid,n=1);

write.table(paste('ORF','gap_start','gap_stop','pattern_after','position_after','CRF',sep='\t',collapse='\t'),outputF1,sep='\t',col.names=F,row.names=F,quote=F,append=F);

while(length(IN)        !=      0){
        tmp             =       unlist(strsplit(IN,split="\t",fixed=T));
        TRID            =       tmp[1];
        frame           =       as.numeric(tmp[2]);
        start           =       as.numeric(tmp[3]);
        stop            =       as.numeric(tmp[4]);
        high3nt         =       as.numeric(unlist(strsplit(tmp[5],split=",",fixed=T)));
        chk_ORF         =       as.numeric(tmp[2:4]);
        ID3nt           =       which(high3nt>0);
        threent_start   =       min(ID3nt);
        threent_stop    =       max(ID3nt);
        frms            =       (ID3nt+2)%%3;
        brkPots         =       my_brkPots(frms);
        if(brkPots>0){
                gaps    	=       my_gap(ID3nt,frms);
                gap_start       =       gaps$gap_start;
		gap_stop        =       gaps$gap_stop;
		ID_tmp		=	intersect(which(gap_start > start), which(gap_start < stop));
		gap_start	=	gap_start[ID_tmp];
		gap_stop	=	gap_stop[ID_tmp];
		if(length(gap_start)>0){
                	for (k in 1:length(gap_start)){
                        	gap_start_tmp   =       gap_start[k];
                        	gap_stop_tmp    =       gap_stop[k];
                        	ID3nt_after     =       ID3nt[ID3nt>gap_stop_tmp];
                        	frms_after      =       (ID3nt_after+2)%%3;
                        	pattern_after   =       my_pattern(frms_after);
                        	position_after  =       my_position(ID3nt_after,frms_after);
                        	gap_out 	=       my_calculation(high3nt,frame,gap_start_tmp,gap_stop_tmp,start,threent_stop);
                        	if(!is.na(gap_out)){
                                	line=paste(paste(TRID,frame,start,stop,sep="_",collapse="_"),gap_start_tmp,gap_stop_tmp,pattern_after,position_after,gap_out,sep='\t',collapse='\t');
                                	write.table(line,outputF1,sep='\t',col.names=F,row.names=F,quote=F,append=T);
                        	}
                	}
        	}
	}
        IN              =       readLines(INid,n=1);
}

####    
close(INid);
rm(list=ls());


