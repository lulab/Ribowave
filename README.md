# RiboWave 

RiboWave is a funtional Ribo-seq analysis tool to identify translated ORF based on Ribo-seq data.


The RiboWave workflow consists of:


- Pre-processing :
  - Create the annotation file for the subsequent analysis. [`create_annotation.sh`]
  - Determine the P-site position of Ribo-seq. [`P-site_determination.sh`]
  - Generate P-site tracks from Ribo-seq, dependent on P-sites calculation cutoffs.[`create_track_Ribo.sh`]

- Main function :
  - Predict translated ORFs [`main_function.sh  pvalue`]
  - Estimate reads density for each given ORF [`main_function.sh  density`]
  - Estimate frameshift potential for each given ORF [`main_function.sh  CRF`]


## Requirements
### software
* R 
* bedtools v2.25.0 

### R packages
* reshape
* ggplot2
* rhdf5
* methods
* wmtsa
* parallel

## Before running 
It is **recommanded** to make a new directory and move the `Ribo-seq bam file` into that directory;

## Pre-processing

### 0. Create annotation

This step scans for and annotates all putative ORFs 

```
Usage: ./create_annotation.sh <genome.gtf> <fasta> <annotation_dir> <scripts_dir>

Example: scripts/create_annotation.sh annotation_yeast  annotation_yeast/Saccharomyces_cerevisiae.R64-1-1.90.gtf  annotation_yeast/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa scripts;
```

#### Input files:
- <annotation.gtf> : the annotation gtf should contain **start_codon** and **stop_codon** information,eg: `Saccharomyces_cerevisiae.R64-1-1.90.gtf` 
- <genome.fasta> : eg: `Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa` 

- <annotation_dir> :  the directory for all the annotation output
- <scripts_dir> 	: the directory of all the scripts in the package

#### Output files:
**`annotation`** directory, including :

* start_codon.bed 	: annotated start codon 

* final.ORFs 	: all ORFs. eg:`YAL001C_0_1_3480` where YAL001C refers to the transcript, 0 refers to the reading frame, 1 refers to the start site, 3480 refers to the stop codon.  

### 1. P-site determination

This step determine the P-site position for each read length by overlapping with the annotated start codon 

```
Usage: ./P-site_determination.sh  <Ribo_bam>  <start_codon.bed> <out_dir> <output_identifier> <scripts_dir>

Example: scripts/P-site_determination.sh  GSE52968/SRR1042853.sort.bam  annotation_yeast/start_codon.bed  GSE52968  SRR1042853  scripts;
```

#### Input files:
- <Ribo_bam> :  **secondary alignment removed** and sorted

- annotation_dir  : 
  - <start_codon.bed> : annotated start site `start_codon.bed` 

- out_dir 	: the directory of the output result, eg: `GSE52968`

- output_name 	: the name of all the output file, default: test. eg: `SRR1042853` 

- scripts_dir 	: the directory of all the scripts in the package


#### Output files:
**`P-site`** directory, including :

* _name_.psite1nt.txt 	: the P-sites position (= offset + 1) for each length. It may look this this : 
  
  ```
  25	10
  26	11
  27	12
  28	13
  29	11
  30	12
  ```
  
* _name_.psite.pdf 	: the pdf displaying the histogram of aggregated reads


### 2. Generating P-site track 

This step creats the P-site track for transcripts of interests

```
Usage: ./create_track_Ribo.sh <Ribo_bam>  <transcripts.exon.gtf>  <genome>  <out_dir> <output_name> <scripts_dir>

Example1: scripts/create_track_Ribo.sh  GSE52968/SRR1042853.sort.bam  annotation_yeast/exons.gtf  annotation_yeast/genome GSE52968  SRR1042853  scripts;

Example2: scripts/create_track_Ribo.sh  GSE52799/SRR1039770.sort.bam  annotation_fly/chrX.exons.gtf annotation_fly/genome GSE52799  SRR1039770  scripts;
```

#### Input files:

- <Ribo_bam> 

- <exons.gtf> :  a gtf file for only the exons from transcripts of intersect, eg: `chrX.exons.gtf`, `exons.gtf`

- <genome\> :  the file including all the chromosomes and its length, `genome` may look like this:
    
    ```
    I 230218
    II  813184
    III 316620
    IV  1531933
    ```
    
- out_dir 	: the directory of the output result, eg: `GSE52968`

- output_name 	: the name of all the output file, default: test. eg: `SRR1042853` 

- scripts_dir 	: the directory of all the scripts in the package

#### Output files:

**`bedgraph`** directory, including :

* final.psite 	: P-site track for each interested transcript. It may look like this : 
  
  ```
  YAL044C 3,0,0,3,10,0,0,0,0,2,0,0,36,0,0,4,0,0,0,6,0,12,0,7,32,0,0,6,7,9,19,2,5,28,0,0,0,0,0,0,0,0,0,0,4,0,0,0,24,0,1,34,0,1,9,2,0,8,0,0,0,0,0,38,0,4,33,0,10,24,0,8,2,0,6,16,0,0,2,0,0,4,0,0,0,2,0,
  YAL045C 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
  YAL046C 3,0,0,2,0,0,4,0,0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,1,0,0,16,0,0,6,0,0,0,0,0,12,0,0,19,0,0,0,0,0,5,0,1,5,0,0,0,0,1,6,0,0,0,0,0,14,0,0,2,0,0,7,0
  ```
  
## Main function

### 3. RiboWave

This step can achieve multiple functions : 

  - providing predicted p value for each given ORF to identify its translation status [`pvalue`]
  
  - providing reads density (P-site/PF P-site) for each given ORF [`density`]
  
  - providing frameshift potential (CRF score) for each given ORF [`CRF`]


```
./main_function.sh -h
```

A helper message will be shown: 

```
----------------------------------------------------------------------------------------------------
RiboWave : version 1.0 
This step is main function of RiboWave.                                               
Functions are provided including : predicting translated ORF, estimating reads density, predicting frameshift events.
----------------------------------------------------------------------------------------------------

Usage:
	 main_function.sh [-h] <job> -a P-site track -b ORF_list -o output_dir -s scripts_dir [-n output_name] [-p core]

Options:
	<job>	<string>  	(pvalue : provide the pvalue for each given ORF;		density : provide reads density for each given ORF;	CRF : provide frameshift potential for each given ORF)
	-a	<filename>	(psite track                                       )
	-b	<filename>	(ORF list                                          )
	-o	<directory>	(Output directory                                  )
	-s	<directory>	(Script directory                                  )
	-n	<string>  	(The name of the output file, default: test        )
	-p	<int>     	(The number of threads, default: 1                 )
	-h	          	(Help                                              )
----------------------------------------------------------------------------------------------------
```


It might take hours to perform the analysis if the input is large. It is **recommended** to specify the number of CPU cores through the `-p` option. 

Run `main_function.sh` on example:

```
scripts/main_function.sh     pvalue	 -a GSE52968/bedgraph/SRR1042853/final.psite -b annotation_yeast/final.ORFs -o GSE52968 -n SRR1042853 -s scripts -p 8;
scripts/main_function.sh     density	 -a GSE52968/bedgraph/SRR1042853/final.psite -b annotation_yeast/final.ORFs -o GSE52968 -n SRR1042853 -s scripts -p 8;
scripts/main_function.sh     CRF 	 -a GSE52968/bedgraph/SRR1042853/final.psite -b annotation_yeast/final.ORFs -o GSE52968 -n SRR1042853 -s scripts -p 8;
```

#### Input files:

- <P-site track\> : output from the previous step, containing the P-site track of transcripts of interest

- <ORF_list> : ORFs of interest ,eg : `final.ORFs`

#### Output files:

* _name_.feats1 	: the features of ORFs including chi-square P-value information. It may look like this :

Column | Explanation	
------------ | -------------
column1-column7 | basic information about the ORF
column8		| reads coverage within the ORF
column9		| P-value predicted by RiboWave
column10	| translational signal outside current ORF
column11	| reads intensity at start codon

**`result`** directory, including :

* _name_.95%.mx 	: the final translated product of RiboWave with translation initiation sites specified. It may look like this :

```
YBR073W_0_103_2874
YBR152W_0_406_873
YBR197C_0_292_651
```

* _name_.COV	: reads density (P-site/PF P-site) of given ORFs. It may look like this : 

Column | Explanation	
------------ | -------------
column1-column6 | basic information about the ORF
column7		| P-site density within the ORF
column8		| denoised PF P-site density within the ORF


* _name_.CRF  : ORFs that might experience reading frame translocation. It may look like this :

Column | Explanation	
------------ | -------------
column1 | ORF
column2		| the start of gap region 
column3		| the stop of gap region
column4		| reading frames after the shift ,eg: `0_7,2_2` where `0_7` refers to continuous seven PF P-sites of frame 0 followed by continuous two PF P-sites of frame 2.
column5		| the corresponding PF P-sites position after the shift ,eg : `1996,2152,2197,2437,2569,2596,2605;2712,2907` where `1996,2152,2197,2437,2569,2596,2605` corresponds to the exact position of `0_7` within the transcript. Discontinuity in the reading frame is separated by `;`.
column6		| F1 score describing the extent of change in the reading frame before and after the gap
column7 	| F2 score describing the extent of reading frame consistency before and after the gap
column8		| Fscore describing the potential of frameshift 



  

