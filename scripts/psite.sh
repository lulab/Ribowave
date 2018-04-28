#!/bin/bash

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


# Print help message if no parameter given
if [ $# -ne 5 ];then
echo "Usage: ./psite.sh input_bam  input_bed	bin_dir		output_dir	header"
exit;
fi
########################	input_parameters
input_bam="`readlink -e $1`";
input_bed="`readlink -e	$2`";
bin_dir="`readlink -e $3`";
output_dir="`readlink -e $4`";
header=$5;

########################	processing
mkdir -p $output_dir
cd $output_dir;

bamToBed -i 	$input_bam  -split 				>	$output_dir/$header.bed;
intersectBed -a $output_dir/$header.bed -b $input_bed -s -wa -wb	>	$output_dir/$header.overlap;
awk 	'{print $4"\t"$9}'	$output_dir/$header.overlap		>	$output_dir/$header.foo1;
sort -k1,1 -k2,2g	$output_dir/$header.foo1	-u		>	$output_dir/$header.foo2;
awk 'BEGIN{foo="";info=$2;}{if(foo==$1){info=info","$2;}else{print foo"\t"info;foo=$1;info=$2;}}END{
print foo"\t"info;}'	$output_dir/$header.foo2	|sed '1d'	>	$output_dir/$header.foo3;
awk '{print $4"\t"$6"\t"$2"\t"$3}'	$output_dir/$header.bed 	>	$output_dir/$header.foo4;
perl	$bin_dir/is_A_in_B.pl		$output_dir/$header.foo4		$output_dir/$header.foo3		$output_dir/$header.foo5;
sort -k1,1 -k2,2 -k3,3g -k4,4g		$output_dir/$header.foo5	>	$output_dir/$header.foo6;
awk 'BEGIN{foo="";info=$3+1":"$4;}{if(foo==($1"\t"$2)){info=info","$3+1":"$4;}else{print foo"\tc("info")";foo=$1"\t"$2;info=$3+1":"$4;}}END{
print foo"\tc("info")";}' $output_dir/$header.foo6|sed '1d'	>	$output_dir/$header.foo7;
paste	$output_dir/$header.foo7		$output_dir/$header.foo3	|awk '$1==$4'	|cut -f 1-3,5	>	$output_dir/$header.foo8;
echo "meta_gene";
Rscript	$bin_dir/metagene_psite.R	$output_dir/$header.foo8		$output_dir/$header.psite.txt	$output_dir/$header.psite.pdf


####	clean
rm -rf		$output_dir/$header.bed	$output_dir/$header.overlap	$output_dir/$header.foo*;

