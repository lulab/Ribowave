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


library(methods,quietly=T)
library(wmtsa,quietly=T)

########        sub functions
####    signal is the raw signal, a numeric vector
##      signal is the raw signal, if shorter than 64 nt, zeros will be padded 
##      higher3nt is the raw signal with lower 3nt positions' intensity flatted to zeros.
##      mx is the DWPT parameter matrix, for each frequency as rows, each signal position's coefficients as columns.
my_DWPT         =       function(signal){
library(wmtsa);
TR_lth          =       length(signal);
if(TR_lth       <       64){
        signal  =       c(signal,rep(0,64-length(signal)));
}
W1              =       wavMODWPT(signal, wavelet="s4",n.levels=6);
W2              =       wavShift(W1);
bands           =       32:51;# use the 0.2~0.5 Hz components only.
mx              =       matrix(0,nrow=length(bands),ncol=length(signal));
for(i in 1:length(bands)){
        tmp     =       paste("w6.",bands[i],sep="");
        mx[i,]  =       W2$data[[tmp]];
}
ID_signal       =       which(signal>0);        #the positions with signal;
mx[,-ID_signal] =       0;                      #remove noise;
minus3nt        =       mx[-11,];
only3nt         =       mx[11,];
BKgrnd          =       apply(minus3nt,2,max);
ID1             =       which(only3nt > BKgrnd);#the positions with 3nt energy higher than other frequency;
higher3nt       =       signal;
higher3nt[-ID1] =       0;
out             =       higher3nt[1:TR_lth];
}

########        sub functions
####    TR_ORF is a two column dataframe, Start nt, Stop nt -1;
####    TR_lth is the length of whole transcript.
##      TR_ORF2 is a five column dataframe, Upstream_Stop, Start, Stop-1, Downstream_Start-1,ORF's first one third or first 30nt, whicever comes first;
ORF_with_flank  =       function(TR_ORF,TR_lth){
stops           =       as.vector(as.matrix(TR_ORF[,2]))+1;
Up_stops        =       unique(c(0,stops));
starts          =       as.vector(as.matrix(TR_ORF[,1]))-1;
Dn_starts       =       unique(c(starts,TR_lth));
TR_ORF2         =       cbind(TR_ORF[,1],TR_ORF,TR_ORF[,2],TR_ORF[,1]);
colnames(TR_ORF2)=      c("UpStop","Start","Stop","DnStart","first30nt");
for(i in 1:nrow(TR_ORF)){
ID1             =       which(Up_stops<TR_ORF[i,1]);
TR_ORF2[i,1]    =       max(Up_stops[ID1]);
ID2             =       which(Dn_starts>TR_ORF[i,2]);
TR_ORF2[i,4]    =       min(Dn_starts[ID2]);
tmp1            =       round((TR_ORF[i,2]-TR_ORF[i,1]+1)/3);
TR_ORF2[i,5]    =       TR_ORF[i,1]+tmp1-1;
}
TR_ORF2;
}

########        sub functions
####    vec is a numeric vector.
##      out is the scaled vector.
my_relative     =       function(vec){
a               =       min(vec);
b               =       max(vec);
L               =       length(vec);
if(a==b){
        if(a==0){
        out     =       rep(0,L);
        }else{
        out     =       rep(1/L,L);
        }
}else{
out             =       (vec-a)/(b-a);
}
out;
}

########        sub functions
####    in_exp is the expected in-orf length
####    OT_exp is the expected out-orf length
####    in_3nt is the 3nt signal's in-orf length
####    OT_3nt is the 3nt signal's out-orf length
my_chisq        =       function(in_exp,OT_exp,in_3nt,OT_3nt){
if(in_3nt*OT_exp > OT_3nt*in_exp){
       tmp      =       c(in_3nt,in_exp,OT_3nt,OT_exp);
       dim(tmp)=        c(2,2);
       if(in_3nt<=5||OT_exp<=5||OT_3nt<=5||in_exp<=5){
        tmp2    =       fisher.test(tmp)
        P       =       tmp2$p.value;
       }else{
        tmp2    =       chisq.test(tmp)
        P       =       tmp2$p.value;
       }
}else{
        P               =       1;
}
P;
}


########        sub functions
####    TR_ORF2 is a two column dataframe, Start, Stop-1;
####    high3nt is the raw signal with lower3nt positions' intensity flatted to zero.
####    signal is the raw signal.
####    main function;
OVLP_ORF        =       function(high3nt,TR_ORF2){
TR_lth          =       length(high3nt);
IDTR            =       1:TR_lth;
ID3nt           =       which(high3nt>0);
N               =       nrow(TR_ORF2);
SUM3nt          =       sum(high3nt);
LEN3nt          =       length(ID3nt);
####    initial values
pvalue		=       rep(0,N);
OUT_3nt    =       rep(0,N);

for(i in 1:N){
IDfrm           =       seq(TR_ORF2[i,2],TR_ORF2[i,3],by=3);
IDORF           =       seq(TR_ORF2[i,2],TR_ORF2[i,3],by=1);
IDOUT           =       IDTR[-IDORF];

ID_in_3nt       =       intersect(IDfrm,ID3nt);
ID_OT_3nt       =       intersect(IDOUT,ID3nt);

in_exp          =       length(IDfrm);
OT_exp          =       length(IDOUT);

in_3nt_pos      =       length(ID_in_3nt);
OT_3nt_pos      =       length(ID_OT_3nt);
in_3nt_int      =       sum(high3nt[ID_in_3nt]);
OT_3nt_int      =       sum(high3nt[ID_OT_3nt]);

# in vs OUT chisq test
pvalue[i]       =       my_chisq(3*in_exp,OT_exp,in_3nt_pos,OT_3nt_pos);
if(LEN3nt>0){
OUT_3nt[i] =       (OT_3nt_pos/LEN3nt) * (OT_3nt_int/SUM3nt);
}
}

out             =       cbind(pvalue,OUT_3nt);
out;
}

########        sub functions
####    TR_ORF is a two column dataframe, Start, Stop-1;
####    signal is the raw signal.
####    main function;
cov_check       =       function(TR_ORF,signal){
TR_lth          =       length(signal);
N               =       nrow(TR_ORF);
CovORF          =       rep(0,N);
for(i in 1:N){
IDORF           =       seq(TR_ORF[i,1],TR_ORF[i,2],by=1);
CovORF[i]       =       length(which(signal[IDORF]>0))/length(IDORF);
}
out             =       cbind(CovORF);
out;
}


########        sub functions
####    TR_ORF is a two column dataframe, Start, Stop-1;
####    signal is the raw signal.
####    main function;
den_check       =       function(TR_ORF,signal){
TR_lth          =       length(signal);
N               =       nrow(TR_ORF);
trans		=	rep(0,N);
ORF		=	rep(0,N);
DenORF          =       rep(0,N);
for(i in 1:N){
trans[i]	=	sum(signal);
IDORF           =       seq(TR_ORF[i,1],TR_ORF[i,2],by=1);
ORF[i]		=	sum(signal[IDORF]);
DenORF[i]       =       sum(signal[IDORF])/length(IDORF);
}
out             =       cbind(trans,ORF,DenORF);
out;
}



########        main function
feat_calc       =       function(TR_ORF,high3nt,signal){
TR_lth          =       length(signal);
TR_ORF2         =       ORF_with_flank(TR_ORF,TR_lth);
out3nt          =       OVLP_ORF(high3nt,TR_ORF2);
N               =       nrow(TR_ORF);
start_psites    =       rep(0,N);
for(i in 1:N){
ID1             =       (TR_ORF[i,1]-1):(TR_ORF[i,1]+1);
ID2             =       ID1[ID1>=1&ID1<=TR_lth];
start_psites[i] =       max(signal[ID2]);
}
startPsite      =       my_relative(start_psites);
OUT             =       cbind(out3nt,startPsite);
OUT;
}



