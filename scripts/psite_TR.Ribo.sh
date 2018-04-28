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


# Print help message if no parameter given
input_bam=$1
input_gtf=$2;
input_offset=$3;
input_genome=$4;
binDir=$5;
outDir=$6;

########################	pre-processing
mkdir -p $outDir
cd $outDir;
file=${input_bam##*/};


awk -F '[\t;]' '{for (x=1;x<=NF;x++) if ($x~"transcript_id") print $x"\t"$1"."$7"\t"$4"\t"$5}'          $input_gtf              |sed 's/transcript_id//g;s/"//g' | awk 'BEGIN{OFS="\t"}{NF=NF;print $0}'|sort -k1,1 -k2,2 -k3,3n -k4,4n       >       $outDir/foo;
awk -F '\t' 'BEGIN{ID="";chr="";coor=""}{if(ID!=$1){print ID"\t"chr"\tc("coor")";ID=$1;chr=$2;coor=$3":"$4;}else{coor=coor","$3":"$4;}}END{print ID"\t"chr"\tc("coor")";}'	$outDir/foo	|sed '1d'								>       $outDir/exons.coor;
bedtools bamtobed  -i	$input_bam -bed12	>	$outDir/$file.bed12;

######## cut by chr and str
rm -rf	$outDir/exons.coor.*	$outDir/$file.bed12.*;
awk -F '\t' -v prefix="$outDir/exons.coor"  '{print $0 >> prefix"."$2}'	$outDir/exons.coor;
awk -F '\t' -v prefix="$outDir/$file.bed12" '{
if($10==1){
	print $11"\tc("$7+1":"$8")"	>>	prefix"."$1"."$6;
}else{
	split($11,a,",");split($12,b,",");len=a[1];coor=$7+1+b[1]":"$7+b[1]+a[1];
	for(i=2;i<=$10;i++){len=len+a[i];coor=coor","$7+1+b[i]":"$7+b[i]+a[i];}
	print len"\tc("coor")"		>>	prefix"."$1"."$6;
}}'	$outDir/$file.bed12;
######## reads2psite
ls $outDir/exons.coor.*		|awk -F '.' '{print $(NF-1)"."$NF}'	>	$outDir/foo;
ls $outDir/$file.bed12.*	|awk -F '.' '{print $(NF-1)"."$NF}'	>>	$outDir/foo;
sort $outDir/foo|uniq -c|awk '$1==2{print $2}'				>	$outDir/chrstr;

for i in `cat $outDir/chrstr`;do 
chr=`echo $i|awk -F '.' '{print $(NF-1)}'`;
str=`echo $i|awk -F '.' '{print $NF}'`; 
chrLen=`awk -v chr="$chr" '$1==chr{print $2}'	$input_genome`;
echo -e "$chr\t$str\t$chrLen";
rm -rf	$outDir/$file.HDF5.$chr.$str		$outDir/exons.psite.$chr.$str;
Rscript	$binDir/HDF5_genome_reads2psite.Ribo.R	$outDir/$file.bed12.$chr.$str	$input_offset	$chrLen		$str	$outDir/$file.HDF5.$chr.$str;
Rscript	$binDir/HDF5_TR_psite_track.R		$outDir/$file.HDF5.$chr.$str	$outDir/exons.coor.$chr.$str	$str	$outDir/exons.psite.$chr.$str;
done

######## final.psite
cat	$outDir/exons.psite.*		>	$outDir/final.psite;
rm -rf	$outDir/foo	$outDir/chrstr	$outDir/exons.*	$outDir/$file.*;
echo "ALLDONE";
