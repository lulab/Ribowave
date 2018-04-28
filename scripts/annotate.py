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

for i in f:
	line=i.strip().split('\t')
	gene='NA'
	trans='NA'
	name='NA'
	gene_biotype='NA'
	transcript_biotype='NA'
	tmp=line[8].split(';')
	for k in tmp:
		if 'gene_id' in k:
			gene=k.split('\"')[1]
		if 'transcript_id' in k:
			trans=k.split('\"')[1]
		if 'gene_name' in k or 'gene_symbol' in k:
			name=k.split('\"')[1]
		if 'gene_type' in k or 'gene_biotype' in k:
			gene_biotype=k.split('\"')[1]
		if 'transcript_type' in k or 'transcript_biotype' in k:
			transcript_biotype=k.split('\"')[1]
	p.write(gene+'\t'+trans+'\t'+gene_biotype+'\t'+name+'\t'+transcript_biotype+'\n')
