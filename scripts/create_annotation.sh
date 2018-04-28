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
        printf '%-90s\n%-90s\n' "${BOLD}This step is set for the purpose of genome annotation." "Several of the output is necessary for the following steps.${NORM}";
        printf -v line '%*s' "100"; echo -e "${line// /-}\n";
        echo -e "${BOLD}Usage:${NORM}\n\t ${SCRIPT} [-h] -G annotation.gtf -f fasta -o annotation_dir -s scripts_dir"\\n;
        echo -e "${BOLD}Options:${NORM}";
        printf "\t%2s\t%-10s\t%1s%-50s%1s\n" "-G" "<filename>" "(" "annotation.gtf" ")"
        printf "\t%2s\t%-10s\t%1s%-50s%1s\n" "-f" "<filename>" "(" "genome fasta" ")"
        printf "\t%2s\t%-10s\t%1s%-50s%1s\n" "-o" "<directory>" "(" "Output annotation directory" ")"
        printf "\t%2s\t%-10s\t%1s%-50s%1s\n" "-s" "<directory>" "(" "Script directory " ")"
        printf "\t%2s\t%-10s\t%1s%-50s%1s\n" "-h" "" "(" "Help" ")"
        printf -v line '%*s' "100";echo ${line// /-}    
	exit 0
}
if [ $# == 0 ] ;then
        HELP;
        exit 1
fi

## Start getopts code ###
while getopts :G:f:o:s:h FLAG; do
	case $FLAG in
		G) #set option gtf "G"
		gtf=$OPTARG
		#echo "-G used:$OPTARG"
		;;
		f) #set fasta sequence "f"
		fa=$OPTARG
		#echo "-f used:$OPTARG"
		;;
		o) #set annoation output directory "o"
		annotation=$OPTARG
		#echo "-o used:$OPTARG"
		;;
		s) #set script directory "s"
		script_dir=$OPTARG
		#echo "-s used:$OPTARG"
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


if ! [[ -f $gtf ]];then
        echo "Error!!!!!   Genome annotation file not found!."
        exit 1
fi
if ! [[ -f $fa ]];then
        echo "Error!!!!!   Genome fasta not found!."
        exit 1
fi
if ! [[ -d $annotation ]];then
        echo "Error!!!!!   Annotation output directory not found!."
        exit 1
fi
if ! [[ -d $script_dir ]];then
        echo "Error!!!!!   script_dir not found!."
        exit 1
fi

#### 
annotation_full=`readlink -e $annotation`

gtf_full=`readlink -e $gtf`

fa_full=`readlink -e $fa`

script_dir_full=`readlink -e $script_dir`

echo "Define ORF ...";
echo "Parameters used:"
echo ""

echo "<input_gtf> 			$gtf_full"
echo "<fasta>				$fa_full";
echo "<genome annotation reference> 	$annotation_full"
echo "<scripts_dir> 			$script_dir_full"
echo ""
echo "---------------"
echo ""

cd $annotation_full;

### find the start codon position(bed format)
echo "------ start codon annotation ------";
awk -F '[\t;]' '$3=="start_codon"'      $gtf_full      >       $annotation_full/transcript.exon.CDS.foo;
awk -F '[\t;]'  '($5-$4)==2{for (x=1;x<=NF;x++) if ($x~"transcript_id")print $1"\t"$4"\t"$5"\t"$x"\t"$7}' $annotation_full/transcript.exon.CDS.foo|sed 's/transcript_id//g;s/"//g'     |awk 'BEGIN{OFS="\t"}{NF=NF;print $0}'       >       $annotation_full/transcript.exon.CDS.foo1;
awk -F '\t'     '{if($5=="+"){print $1"\t"$2-1"\t"$2"\t"$4"\t0\t"$5;}else{print $1"\t"$3-1"\t"$3"\t"$4"\t0\t"$5}}'      $annotation_full/transcript.exon.CDS.foo1 >       $annotation_full/start_codon.bed;

### scan for the ORFs
echo "------ scan ORFs ------";
awk -F '\t' '$3=="exon"{print $0}'	$gtf_full	>	$annotation_full/exons.gtf;
perl    $script_dir_full/gtf2Bed.pl       $annotation_full/exons.gtf                 >       $annotation_full/exons.bed;
###    exon fasta and peptide;
bedtools getfasta -fi 	$fa_full    -name -s -split -bed $annotation_full/exons.bed   -fo $annotation_full/exons.fa;
####    getorf
awk -F '\t' '{if($1~/>/){if(NR>1){print foo;}foo=$0"\t";}else{foo=foo""$1;}}
END{print foo}'     $annotation_full/exons.fa                      >       $annotation_full/foo;
python  $script_dir_full/getorf.py    $annotation_full/foo                       $annotation_full/exons.ORFs.fa;
awk -F '\t' '$1~/>/{gsub(">","",$1);b=split($1,a,"_");foo=a[1];for(i=2;i<b-2;i++){foo=foo"_"a[i];}
print foo","a[3]","a[4]}'       $annotation_full/exons.ORFs.fa >       $annotation_full/exons.ORFs;
rm -rf $annotation_full/*foo*;

####   annotation the genome file
python $script_dir_full/annotate.py	$annotation_full/exons.gtf	$annotation_full/annotate.foo;
cat $annotation_full/annotate.foo|sort|uniq 	>	$annotation_full/annotation.list;


####    annotate ORF
####    1. all          <       unannotated     vs      annotated
echo "======  annotate ORF ======";
awk -F '[\t;]' '$3=="exon"||$3=="start_codon"||$3=="stop_codon"'        $gtf_full       >       $annotation_full/transcript.exon.CDS.foo;
awk -F '[\t;]' '{for (x=1;x<=NF;x++) if ($x~"transcript_id")print $x"\t"$3"\t"$4"\t"$5"\t"$7}'       $annotation_full/transcript.exon.CDS.foo  |sed 's/transcript_id//g;s/"//g'     |awk 'BEGIN{OFS="\t"}{NF=NF;print $0}'	>     $annotation_full/transcript.exon.CDS.foo1;
awk -F '[\t]'  '$2=="exon"'     $annotation_full/transcript.exon.CDS.foo1|sort -k1,1 -k3,3g -k4,4g|awk '{tmp=$1"\t"$5;
if(NR==1){foo=tmp;info=$3":"$4;}else{if(tmp==foo){info=info","$3":"$4;}else{print foo"\tc("info")";foo=tmp;info=$3":"$4;}
}}END{print foo"\tc("info")"}'                  >          $annotation_full/transcript.exon.CDS.foo2.exon;
awk -F '[\t]'  '$2=="start_codon"'                         $annotation_full/transcript.exon.CDS.foo1|sort -k1,1 -k3,3g -k4,4g|awk '{tmp=$1"\t"$5;
if(NR==1){foo=tmp;info=$3":"$4;}else{if(tmp==foo){info=info","$3":"$4;}else{print foo"\tc("info")";foo=tmp;info=$3":"$4;}
}}END{print foo"\tc("info")"}'                  >          $annotation_full/transcript.exon.CDS.foo2.start;
awk -F '[\t]'  '$2=="stop_codon"'                          $annotation_full/transcript.exon.CDS.foo1|sort -k1,1 -k3,3g -k4,4g|awk '{tmp=$1"\t"$5;
if(NR==1){foo=tmp;info=$3":"$4;}else{if(tmp==foo){info=info","$3":"$4;}else{print foo"\tc("info")";foo=tmp;info=$3":"$4;}
}}END{print foo"\tc("info")"}'                  >          $annotation_full/transcript.exon.CDS.foo2.stop;
#awk -F ',' '{print NF}' $annotation_full/transcript.exon.CDS.foo2.start   |sort|uniq -c
#awk -F ',' '{print NF}' $annotation_full/transcript.exon.CDS.foo2.stop    |sort|uniq -c
cat     $annotation_full/transcript.exon.CDS.foo2.*       |cut -f1|sort|uniq -c|awk '$1==3{print $2}'     >       $annotation_full/transcript.exon.CDS.foo3;
perl    $script_dir_full/get_matrix_from_ID.pl    $annotation_full/transcript.exon.CDS.foo2.exon    $annotation_full/transcript.exon.CDS.foo3 $annotation_full/transcript.exon.CDS.foo4.exon;
perl    $script_dir_full/get_matrix_from_ID.pl    $annotation_full/transcript.exon.CDS.foo2.start   $annotation_full/transcript.exon.CDS.foo3 $annotation_full/transcript.exon.CDS.foo4.start;
perl    $script_dir_full/get_matrix_from_ID.pl    $annotation_full/transcript.exon.CDS.foo2.stop    $annotation_full/transcript.exon.CDS.foo3 $annotation_full/transcript.exon.CDS.foo4.stop;


### select annotated ORF
paste   $annotation_full/transcript.exon.CDS.foo4.exon    $annotation_full/transcript.exon.CDS.foo4.start   $annotation_full/transcript.exon.CDS.foo4.stop    |\
awk -F '\t' '$1==$4 && $1==$7{print $1"\t"$2"\t"$3"\t"$6"\t"$9}'                          >       $annotation_full/transcript.exon.CDS.foo4;
Rscript         $script_dir_full/annotORF.R                       $annotation_full/transcript.exon.CDS.foo4         $annotation_full/transcript.exon.CDS.foo5;
awk -F '\t' '(($4-$3+1)%3)==0{print $1","$3","$4}'      $annotation_full/transcript.exon.CDS.foo5 >       $annotation_full/transcript.annotated.ORF;
rm -rf  $annotation_full/*foo*;
wc -l $annotation_full/transcript.annotated.ORF;

### annotate ORF by its start codon (ATG,NTG);
awk -F '\t' '{if($1~/>/){if(NR>1){print foo;}foo=$0"\t";}else{foo=foo""$1;}}
END{print foo}'     $annotation_full/exons.fa                      >       $annotation_full/foo;
awk -F '>' '{print $2}' $annotation_full/foo      >       $annotation_full/foo2;
awk -F ',' '{print $1"\t"$1"_"($2+2)%3"_"$2"_"$3}' $annotation_full/exons.ORFs    >       $annotation_full/exons.ORFs.foo;
perl $script_dir_full/get_matrix_from_ID.pl       $annotation_full/foo2     $annotation_full/exons.ORFs.foo   $annotation_full/foo3;
paste $annotation_full/exons.ORFs.foo     $annotation_full/foo3     |awk -F '\t' '$1==$3{print $2"\t"$NF}'  >       $annotation_full/exons.ORFs.fa.2;
awk -F '\t' 'split($1,a,"_"){print $1"\t"substr($2,a[3],3)}'    $annotation_full/exons.ORFs.fa.2  >       $annotation_full/exons.ORFs.NTG;


### annotate ORF by NTG , coding type, ORF type
awk -F '\t' 'split($1,a,"_") && a[4]-a[3]+1>=30{print $0}'      $annotation_full/exons.ORFs.NTG   >       $annotation_full/exons.ORFs.2;
awk -F ',' '{print $1"_"($2+2)%3"_"$2"_"$3}' $annotation_full/transcript.annotated.ORF    >       $annotation_full/transcript.annotated.ORF.2;
perl $script_dir_full/is_A_not_in_B.pl    $annotation_full/exons.ORFs.2     $annotation_full/transcript.annotated.ORF.2       $annotation_full/unanno.foo;
awk -F '\t' '{print $1"\tunanno\t"$2}'  $annotation_full/unanno.foo       >       $annotation_full/unanno.foo2;

#1. anno
perl $script_dir_full/is_A_in_B.pl        $annotation_full/exons.ORFs.2     $annotation_full/transcript.annotated.ORF.2       $annotation_full/anno.foo;
awk -F '\t' '{print $1"\tanno\t"$2}'    $annotation_full/anno.foo         >       $annotation_full/anno.foo2;

#2. sharestop
awk -F '\t' 'split($1,a,"_"){print a[1]"_"a[4]}'        $annotation_full/transcript.annotated.ORF.2       >       $annotation_full/anno_stop.foo;
awk -F '\t' 'split($1,a,"_"){print a[1]"_"a[4]"\t"$0}'  $annotation_full/unanno.foo2                      >       $annotation_full/unanno.foo3;
perl $script_dir_full/is_A_in_B.pl        $annotation_full/unanno.foo3      $annotation_full/anno_stop.foo    $annotation_full/unanno.foo4;
awk -F '\t' '{print $2"\tsharestop\t"$4}'       $annotation_full/unanno.foo4      >       $annotation_full/sharestop.foo2;
perl $script_dir_full/is_A_not_in_B.pl    $annotation_full/unanno.foo3      $annotation_full/anno_stop.foo    $annotation_full/unanno.foo5;
cut -f 2-       $annotation_full/unanno.foo5      >       $annotation_full/unanno.foo2;
cat $annotation_full/anno.foo2    $annotation_full/sharestop.foo2   $annotation_full/unanno.foo2      >       $annotation_full/final.annot.1;

#3. transcript type
awk -F' \t' 'split($1,a,"_"){print a[1]"\t"$0}' $annotation_full/final.annot.1    >       $annotation_full/final.annot.foo;
awk -F '\t' '{print $2"\t"$0}'  $annotation_full/annotation.list 		  >       $annotation_full/foo;
perl $script_dir_full/is_A_in_B.pl	$annotation_full/final.annot.foo	$annotation_full/foo	$annotation_full/final.annot.foo2;
perl $script_dir_full/get_matrix_from_ID.pl       $annotation_full/foo      $annotation_full/final.annot.foo2  $annotation_full/foo2;
paste   $annotation_full/final.annot.foo2  $annotation_full/foo2     |awk -F '\t' '$1==$5{print $2"\t"$3"\t"$4"\t"$(NF-2)"\t"$NF}'   >       $annotation_full/final.annot.2;
awk -F '\t' 'split($1,a,"_"){if($4=="protein_coding" && $5=="protein_coding"){print $1"\tcoding\t"$2"\t"$3"\t"a[3]"\t"a[4]}else{print $1"\tNC\t"$2"\t"$3"\t"a[3]"\t"a[4]}}'      $annotation_full/final.annot.2    >       $annotation_full/final.ORFs.annot3;

##### considering ATG only
awk -F '\t' '$3=="anno"{print $1}'	$annotation_full/final.ORFs.annot3	>	$annotation_full/aORF;
awk -F '\t' '$4=="ATG"{print $0}'	$annotation_full/final.ORFs.annot3	>	$annotation_full/final.ORFs;
rm -rf $annotation_full/*final.annot.*;
rm -rf $annotation_full/*foo*;
rm -rf $annotation_full/exons.ORFs.fa*;
