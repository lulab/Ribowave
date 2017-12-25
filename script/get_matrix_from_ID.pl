#!/usr/bin/perl
#
#####################################################################
###    This file is part of RiboWave.
###    RiboWave is powerful Ribo-seq analysis tool that is able to 
###    denoise the Ribo-seq data and serve for multiple functions.
###       
###    RiboWave can be used for multiple purposes: 
###               1. denoise the raw Ribo-seq data
###               2. define translated ORFs
###               3. estimate the abundance of actively elongating ribosomes
###               4. estimate translation efficiency(TE)
###               5. identify potential frameshift candidates
###    
###    Author: Long Hu
###
###    Copyright (C) 2017  Zhiyu Xu
###
###    RiboWave is distributed in the hope that it will be useful,
###    but WITHOUT ANY WARRANTY; without even the implied warranty of
###    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
###
###    Contact: xanthexu18@gmail.com
#########################################################################
##
#
use strict;
use warnings;


### header not considered here, the output file will inherit the order of row of the input file. 
open FileTarget	,	"$ARGV[0]" or die; ##first column is the ID, second and latter is one value of prob, all line are outputted
open FileLinenum,	"$ARGV[1]" or die; ##first column is the IDs.
open OUT	,	">$ARGV[2]" or die;##output matched lines 
my $line	=	"";
my %hashTarget;

#process the target file into a hash. prepare for the index search. first is the binID, second 
while($line=<FileTarget>){
	chomp($line);
	my @arr	=	split /\t/,	$line;
	$hashTarget{$arr[0]}=$line;
}
#read the transcripts.binID file, binID as index and ouput corresponding line. 
while($line=<FileLinenum>){
        chomp($line);
	my @arr =       split /\t/,     $line;
        if(defined($hashTarget{$arr[0]})){
                print OUT "$hashTarget{$arr[0]}\n";
        }
}


