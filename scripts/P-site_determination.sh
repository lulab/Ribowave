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
        printf '%-90s\n%-90s\n' "${BOLD}This step is for determining the P-site position for each read length." "${NORM}";
        printf -v line '%*s' "100"; echo -e "${line// /-}\n";
        echo -e "${BOLD}Usage:${NORM}\n\t ${SCRIPT} [-h] -i Ribo.bam -S start_codon.bed -o out_dir -s scripts_dir [-n study name]"\\n;
        echo -e "${BOLD}Options:${NORM}";
        printf "\t%2s\t%-10s\t%1s%-50s%1s\n" "-i" "<filename>" "(" "Ribo-seq bam" ")"
        printf "\t%2s\t%-10s\t%1s%-50s%1s\n" "-S" "<filename>" "(" "start codon annotation file" ")"
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
while getopts :i:S:o:s:n:h FLAG; do
        case $FLAG in
                i) #set Ribo-seq bam "i"
                bam=$OPTARG
                #echo "-i used:$OPTARG"
                ;;
                S) #set start codon  "f"
                start_codon=$OPTARG
                #echo "-S used:$OPTARG"
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
if ! [[ -f $start_codon ]];then
        echo "Error!!!!!   Start codon annotation file not found!."
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


####
bam_full="`readlink -e $bam`";

start_codon_full="`readlink -e $start_codon`";

out_dir_full="`readlink -e $out_dir`";

script_dir_full=`readlink -e $script_dir`


echo "--- determining P-sites ---";
date;
echo "Parameters used:"
echo ""


echo "<bam> 				$bam_full"
echo "<start codon> 			$start_codon_full"
echo "<output folder> 		$out_dir_full"
echo "<study name>			$name"
echo "<scripts directory> 		$script_dir_full"
echo ""
echo "---------------"
echo ""

echo "creating directory..."
mkdir -p $out_dir_full;

cd $out_dir_full;

####################### determine P-site position 
echo "overlapping reads with the annotated start codon";
mkdir -p $out_dir_full/P-site;
bash $script_dir_full/psite.sh    $bam_full    $start_codon_full         $script_dir_full  $out_dir_full/P-site $name;
sed 's/ nt reads//g'    $out_dir_full/P-site/$name.psite.txt   >       $out_dir_full/P-site/foo
Rscript $script_dir_full/psite_1nt_wholeReads.R   $out_dir_full/P-site/foo              $out_dir_full/P-site/foo2;
awk '$2>5000&&$3>0.01'  $out_dir_full/P-site/foo2     |cut -f 1,4-    >       $out_dir_full/P-site/foo3;
awk '{for(i=2;i<=NF;i++){if($i==1){print $1"\t"(i-1)}}}'        $out_dir_full/P-site/foo3     >       $out_dir_full/P-site/$name.psite1nt.txt;
echo "P-site determination output ..."
echo ""
echo "---------------"
echo ""
head -n 20                              $out_dir_full/P-site/$name.psite1nt.txt;
rm -rf  $out_dir_full/P-site/*foo*;

echo ""
echo "---------------"
echo ""
echo "--- determining P-sites, Done! ---";
date
