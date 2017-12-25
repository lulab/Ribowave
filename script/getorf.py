#!/usr/bin/env python

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


import sys

f=open(sys.argv[1],'r')
p=open(sys.argv[2],'w')

def ORF_detection(seq,n):
	seq=seq.upper()
	codon=[seq[i:i+n] for i in range(0, len(seq), n)]
	stop_indice=[i for i, x in enumerate(codon) if (x=='TAA' or x=='TGA' or x=='TAG')]	
	ORF=[]
	if stop_indice:
		string=codon[:stop_indice[0]]
		start_indice=[j for j,x in enumerate(string) if x in ['ATG','CTG','GTG','TTG']]
		ORF=ORF+[[j*3+1,stop_indice[0]*3,''.join(string[j:])]for j in start_indice]
		for i in range(1,len(stop_indice)):
			string=codon[stop_indice[i-1]+1:stop_indice[i]]
			start_indice=[j for j,x in enumerate(string) if x in ['ATG','CTG','GTG','TTG']]
			ORF=ORF+[[(stop_indice[i-1]+j+1)*3+1,stop_indice[i]*3,''.join(string[j:])] for j in start_indice]
	return ORF 
		
				
for i in f.readlines():
	line=i.strip().split('\t')
	id=line[0]
	fa=line[1]
	for i in range(0,3):
		ORF=ORF_detection(fa[i:],3)
		for j in ORF:
			if not 'N' in j[2]:
				p.write(id+'_'+str(i)+'_'+str(int(j[0])+i)+'_'+str(int(j[1])+i)+'\n'+j[2]+'\n')

p.close()
f.close()
