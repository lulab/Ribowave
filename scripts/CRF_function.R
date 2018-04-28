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


library(methods,quietly=T)
library(wmtsa,quietly=T)

########        sub functions
## to define how many breakpoints in the sequence
my_brkPots      =       function(frms){
count           =       0;
if(length(frms)>1){
for(j in 2:length(frms)){
if(frms[j]!=frms[j-1]){
count           =       count+1;
}
}
}
count;
}


## to define the distribution of all the signals by frm
my_pattern      =       function(frms){
L               =       length(frms);
tmp             =       c();
foo             =       frms[1];
count_foo       =       1;
if(L>1){
	for(i in 2:L){
		if(frms[i]==foo){
			count_foo       =       count_foo+1;
		}else{
			tmp2            =       paste(foo,count_foo,sep="_");
			tmp             =       c(tmp,tmp2);
			foo             =       frms[i];
			count_foo       =       1;
		}
	}
}
tmp2            =       paste(foo,count_foo,sep="_");
tmp             =       c(tmp,tmp2);
out             =       paste(tmp,collapse=",");
out;
}

## to define the distribution of all position by its frm
my_position     =       function(ID3nt,frms){
L               =       length(ID3nt);
tmp             =       c();
foo             =       frms[1];
tmp2            =       ID3nt[1];
if(L    > 1){
	for (i in 2:L){
		if(frms[i]==foo){
			tmp2    =       paste(tmp2,ID3nt[i],sep=',');
		}else{
			tmp2    =       paste(tmp2,ID3nt[i],sep=';');
			foo     =       frms[i];
		}
	}
}
tmp2;
}


## to define the gap region for each breakpoint
my_gap     =       function(ID3nt,frms){
L               =       length(ID3nt);
gap_start       =       c();
gap_stop        =       c();
foo             =       frms[1];
tmp2            =       ID3nt[1];
for (i in 2:L){
        if(frms[i]==foo){
                tmp2    =       ID3nt[i];
        }else{
                gap_start    =       c(gap_start,tmp2+1);
                gap_stop     =       c(gap_stop,ID3nt[i]-1);
                foo     =      frms[i];
                tmp2    =      ID3nt[i];
        }
}
out     =       list(gap_start=gap_start,gap_stop=gap_stop);
out;
}


Percentage      =       function(psite,start,stop){
        F0              =               seq(1,length(psite),by=3);
        F1              =               seq(2,length(psite),by=3);
        F2              =               seq(3,length(psite),by=3);
        psite_tot       =               psite[start:stop];
        psite_F0        =               psite[intersect(start:stop,F0)];
        psite_perF0     =               sum(psite_F0)/sum(psite_tot);
        psite_F1        =               psite[intersect(start:stop,F1)];
        psite_perF1     =               sum(psite_F1)/sum(psite_tot);
        psite_F2        =               psite[intersect(start:stop,F2)];
        psite_perF2     =               sum(psite_F2)/sum(psite_tot);
        out             =               c(psite_perF0,psite_perF1,psite_perF2);
        out
}


density =       function(psite,start,stop){
        F0              =               seq(1,length(psite),by=3);
        F1              =               seq(2,length(psite),by=3);
        F2              =               seq(3,length(psite),by=3);
        psite_F0        =               psite[intersect(start:stop,F0)];
        sum_F0          =               sum(psite_F0);
        psite_F1        =               psite[intersect(start:stop,F1)];
        sum_F1          =               sum(psite_F1);
        psite_F2        =               psite[intersect(start:stop,F2)];
        sum_F2          =               sum(psite_F2);
        tot             =               sum(sum_F0,sum_F1,sum_F2);
        out             =               c(sum_F0/tot,sum_F1/tot,sum_F2/tot);
        out;
}


F_score =       function(a,b){
        score   =       2*a*b/(a+b);
        score;
}


## to calculate the frame difference
#  normalize to 0-1 
F1score =       function(F_before,F_after){ 
        F_tmp  =       F_score(F_before,F_after);
        F      =       1-F_tmp;
        F;
}


## to calculate the frame dominancy
#  normalize to 0-1
F2score =       function(F_before,F_after){
        F      =       F_score(F_before,F_after);
        F;
}


## main_function
# here the signal meant the 3nt signal specially for the detemination of frameshift events potentially
my_calculation  =       function(signal,frame,gap_start,gap_stop,start,stop){
        if(sum(signal)>0 && gap_start-1 > start && gap_stop+1 < stop && sum(signal[start:(gap_start-1)])>0 && sum(signal[(gap_stop+1):stop]>0)>0 ){    
                signal_per_before          =       Percentage(signal,start,gap_start-1)[frame+1];
                signal_per_after           =       Percentage(signal,gap_stop+1,stop)[frame+1];
                signal_den_before          =       max(density(signal,start,gap_start-1));
                signal_den_after           =       max(density(signal,gap_stop+1,stop));
                if(signal_den_before >= 0.5 && signal_den_after >= 0.5 && signal_per_before >0){ # before coverage should be larger than 0 and density from both before and after(frame regardless) should be greater than 0.5)
                        signal_den_before_tmp  =       2*signal_den_before-1;
                        signal_den_after_tmp   =       2*signal_den_after-1;
                        score1              =       F1score(signal_per_before,signal_per_after);
                        score2              =       F2score(signal_den_before_tmp,signal_den_after_tmp);
                        final_score         =       score1*score2;
                        out                 =       final_score;
                }else{
                        out                 =       NA;
                }
        }else{
                out     =       NA;
        }
        out;
}



