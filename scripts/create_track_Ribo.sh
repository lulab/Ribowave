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

#------
#Set fonts for Help.
#------
NORM=`tput sgr0`
BOLD=`tput bold`
REV=`tput smso`
#get script name
#------
SCRIPT=`basename ${BASH_SOURCE[0]}`
#------
#HELP function
#------
function HELP
{
        printf -v line '%*s' "100"; echo -e "${line// /-}";
        echo -e "${BOLD}RiboWave : version 1.0 ${NORM}";
        printf '%-90s\n%-90s\n' "${BOLD}This step is the last step for data processing where P-site track is generated for each transcript." "${NORM}";
        printf -v line '%*s' "100"; echo -e "${line// /-}\n";
        echo -e "${BOLD}Usage:${NORM}\n\t ${SCRIPT} [-h] -i Ribo.bam -G exons.gtf -g genome_size -P psite.position -o out_dir -s scripts_dir [-n study name]"\\n;
        echo -e "${BOLD}Options:${NORM}";
        printf "\t%2s\t%-10s\t%1s%-50s%1s\n" "-G" "<filename>" "(" "Ribo-seq bam" ")"
        printf "\t%2s\t%-10s\t%1s%-50s%1s\n" "-g" "<filename>" "(" "Genome size annotation file" ")"
	printf "\t%2s\t%-10s\t%1s%-50s%1s\n" "-P" "<filename>" "(" "P-site position (offset+1)" ")"
        printf "\t%2s\t%-10s\t%1s%-50s%1s\n" "-o" "<directory>" "(" "Output annotation directory" ")"
        printf "\t%2s\t%-10s\t%1s%-50s%1s\n" "-s" "<directory>" "(" "Script directory " ")"
        printf "\t%2s\t%-10s\t%1s%-50s%1s\n" "-n" "<string>" "(" "study name, defult:test " ")"
        printf "\t%2s\t%-10s\t%1s%-50s%1s\n" "-h" "" "(" "Help" ")"
        printf -v line '%*s' "100";echo ${line// /-}    
        exit 0
}
if [ $# == 0 ] ;then
        HELP;
        exit 1
fi

## Start getopts code ###
while getopts :i:G:g:P:o:s:n:h FLAG; do
        case $FLAG in
                i) #set Ribo-seq bam "i"
                bam=$OPTARG
                #echo "-i used:$OPTARG"
                ;;
                G) #set exon gtf "G"
                gtf=$OPTARG
                #echo "-G used:$OPTARG"
                ;;
                g) #set genome size  "f"
                genome=$OPTARG
                #echo "-g used:$OPTARG"
                ;;
		P) #set P-site position  "f"
                P_site=$OPTARG
                #echo "-P used:$OPTARG"
                ;;
		o) #set output directory "o"
                out_dir=$OPTARG
                #echo "-o used:$OPTARG"
                ;;
                s) #set script directory "s"
                script_dir=$OPTARG
                #echo "-s used:$OPTARG"
                ;;
                n) #set script directory "s"
                name=$OPTARG
                #echo "-n used:$OPTARG"
                ;;
                h) #set hellp
                HELP
                ;;
                \?) #unrecognized option - show help;
                HELP
                ;;
        esac
done
shift $((OPTIND-1))


### End getopts code ###
### Main loop to process files ###
while [ $# -ne 0 ]; do
  FILE=$1
  #TEMPFILE=`basename $FILE`
  TEMPFILE="${FILE##*/}"  #
  FILE_BASE=`echo "${TEMPFILE%.*}"`  #file without extension
  FILE_EXT="${TEMPFILE##*.}"
  shift  #Move on to next input file.
done
### End main loop ###

if ! [[ -f $bam ]];then
        echo "Error!!!!!   Ribo-seq bam not found!."
        exit 1
fi
if ! [[ -f $gtf ]];then
        echo "Error!!!!!   Transcript exon annotation file not found!."
        exit 1
fi
if ! [[ -f $genome ]];then
        echo "Error!!!!!   Genome size not found!."
        exit 1
fi
if ! [[ -f $P_site ]];then
        echo "Error!!!!!   P-site position not specified!."
        exit 1
fi
if ! [[ -d $out_dir ]];then
        echo "Error!!!!!   Output directory not found!."
        exit 1
fi
if ! [[ -d $script_dir ]];then
        echo "Error!!!!!   Script_dir not found!."
        exit 1
fi
if ! [[ $name ]];then
        name="test";
fi


bam_full="`readlink -e $bam`";

gtf_full="`readlink -e $gtf`";

genome_full="`readlink -e $genome`";

P_site_full="`readlink -e $P_site`";

out_dir_full="`readlink -e $out_dir`";

script_dir_full=`readlink -e $script_dir`


echo "--- generating signal track ---";
date;
echo "Parameters used:"
echo ""


echo "<bam> 				$bam_full"
echo "<exon gtf> 			$gtf_full"
echo "<genome size> 			$genome_full"
echo "<P_site>					$P_site_full"
echo "<output folder> 		$out_dir_full"
echo "<study name>			$name"
echo "<scripts_dir> 			$script_dir_full"
echo ""
echo "---------------"
echo ""

echo "creating directory..."
mkdir -p $out_dir_full;

cd $out_dir_full;

######################## generate bedgraph
echo "bedgraph tracking";
mkdir -p $out_dir_full/bedgraph;
rm -rf $out_dir_full/bedgraph/$name;
mkdir -p $out_dir_full/bedgraph/$name;
bash    $script_dir_full/psite_TR.Ribo.sh $bam_full     $gtf_full        $P_site_full     $genome_full     $script_dir_full  $out_dir_full/bedgraph/$name;

echo "--- generating signal track, Done! ---";
date;
