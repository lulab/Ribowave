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
./create_annotation.sh -h
```

A helper message is shown:

```
----------------------------------------------------------------------------------------------------
RiboWave : version 1.0 
This step is set for the purpose of genome annotation.                                
Several of the output is necessary for the following steps.                         
----------------------------------------------------------------------------------------------------

Usage:
	 create_annotation.sh [-h] -G annotation.gtf -f fasta -o annotation_dir -s scripts_dir

Options:
	-G	<filename>	(annotation.gtf                                    )
	-f	<filename>	(genome fasta                                      )
	-o	<directory>	(Output annotation directory                       )
	-s	<directory>	(Script directory                                  )
	-h	          	(Help                                              )
----------------------------------------------------------------------------------------------------
```

Run `create_annotation.sh` on example:

```
scripts/create_annotation.sh -G annotation_yeast/Saccharomyces_cerevisiae.R64-1-1.90.gtf -f annotation_yeast/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa -o annotation_yeast -s scripts;
```

#### Input files:
- <annotation.gtf> : the annotation gtf should contain **start_codon** and **stop_codon** information,eg: `Saccharomyces_cerevisiae.R64-1-1.90.gtf` 
- <genome.fasta> : eg: `Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa` 

- <annotation_dir> :  the directory for all the annotation output

- <scripts_dir> 	: the directory of all the scripts in the package

#### Output files:

**`annotation`** directory, including :

* start_codon.bed 	: annotated start codon 

* final.ORFs 	: all ORFs, eg: `YAL001C_0_1_3480` where YAL001C refers to the transcript, 0 refers to the reading frame, 1 refers to the start site, 3480 refers to the stop codon.  

### 1. P-site determination

This step determine the P-site position for each read length by overlapping with the annotated start codon 

```
./P-site_determination.sh -h
```

A helper message is shown:

```
----------------------------------------------------------------------------------------------------
RiboWave : version 1.0 
This step is for determining the P-site position for each read length.                
                                                                                    
----------------------------------------------------------------------------------------------------

Usage:
	 P-site_determination.sh [-h] -i Ribo.bam -S start_codon.bed -o out_dir -s scripts_dir [-n study name]

Options:
	-i	<filename>	(Ribo-seq bam                                      )
	-S	<filename>	(start codon annotation file                       )
	-o	<directory>	(Output annotation directory                       )
	-s	<directory>	(Script directory                                  )
	-n	<string>  	(study name, defult:test                           )
	-h	          	(Help                                              )
----------------------------------------------------------------------------------------------------
```

Run `P-site_determination.sh` on example :

```
scripts/P-site_determination.sh  -i GSE52968/SRR1042853.sort.bam  -S annotation_yeast/start_codon.bed  -o GSE52968  -n SRR1042853  -s scripts;
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
./create_track_Ribo.sh -h
```

A helper message is shown:

```
----------------------------------------------------------------------------------------------------
RiboWave : version 1.0 
This step is the last step for data processing where P-site track is generated for each transcript.
                                                                                    
----------------------------------------------------------------------------------------------------

Usage:
	 create_track_Ribo.sh [-h] -i Ribo.bam -G exons.gtf -g genome_size -P psite.position -o out_dir -s scripts_dir [-n study name]

Options:
	-G	<filename>	(Ribo-seq bam                                      )
	-g	<filename>	(Genome size annotation file                       )
	-P	<filename>	(P-site position (offset+1)                        )
	-o	<directory>	(Output annotation directory                       )
	-s	<directory>	(Script directory                                  )
	-n	<string>  	(study name, defult:test                           )
	-h	          	(Help                                              )
----------------------------------------------------------------------------------------------------
```

Run `create_track_Ribo.sh` on example:

```
scripts/create_track_Ribo.sh  -i GSE52968/SRR1042853.sort.bam  -G annotation_yeast/exons.gtf  -g annotation_yeast/genome -P GSE52968/P-site/SRR1042853.psite1nt.txt  -o GSE52968 -n SRR1042853 -s scripts;
scripts/create_track_Ribo.sh  -i GSE52799/SRR1039770.sort.bam  -G annotation_fly/chrX.exons.gtf -g annotation_fly/genome -P GSE52799/P-site/SRR1039770.psite1nt.txt  -o GSE52799 -n SRR1039770 -s scripts;
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

- <P-site_position> :  the file listing the P-site position for each read length. This file can be found in the output of previous step,  _name_.psite1nt.txt

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

  - denoising [`denoise`]

  - providing predicted p value for each given ORF to identify its translation status [`pvalue`]
  
  - providing reads density (P-site/PF P-site) for each given ORF [`density`]
  
  - providing translation efficiency (TE) estimation for each given ORF [`TE`]
  
  - providing frameshift potential (CRF score) for each given ORF [`CRF`]
  

```
./Ribowave -h
```

A helper message will be shown: 

```
----------------------------------------------------------------------------------------------------
RiboWave : version 1.0 
This step is main function of RiboWave.                                               
Functions are provided including : predicting translated ORF, estimating reads density, estimating translation efficiency and predicting frameshift events.
----------------------------------------------------------------------------------------------------

Usage:
	 Ribowave [Options] -a P-site track -b ORF_list -o output_dir -s scripts_dir

Options:
	-P	                    	(providing P value for each ORF                    )
	-D	                    	(providing reads abundance for each ORF            )
	-F	                    	(predicting frameshift potential for each ORF      )
	-T	<int>  <RNA_FPKM>   	(estimating TE for each ORF                        )
	-a	<filename>          	(psite track                                       )
	-b	<filename>          	(ORF list                                          )
	-o	<directory>         	(Output directory                                  )
	-s	<directory>         	(Script directory                                  )
	-n	<string>            	(The name of the study, default: test              )
	-p	<int>               	(The number of threads, default: 1                 )
	-h	                    	(Help                                              )
----------------------------------------------------------------------------------------------------
```


It might take hours to perform the analysis if the input is large. It is **recommended** to specify the number of CPU cores through the `-p` option. 

Run `Ribowave` on example:

##### Denoise the P-site track

```
scripts/main_function.sh -a GSE52968/bedgraph/SRR1042853/final.psite -b annotation_yeast/final.ORFs -o GSE52968 -n SRR1042853 -s scripts -p 8;
```

##### Identifying translated ORF

```
scripts/main_function.sh -P -a GSE52968/bedgraph/SRR1042853/final.psite -b annotation_yeast/final.ORFs -o GSE52968 -n SRR1042853 -s scripts -p 8;
```

##### Estimating abundance 

```
scripts/main_function.sh -D -a GSE52968/bedgraph/SRR1042853/final.psite -b annotation_yeast/final.ORFs -o GSE52968 -n SRR1042853 -s scripts -p 8;
```

##### Estimating TE

```
scripts/main_function.sh -T 12378563 GSE52968/mRNA/SRR1042851.FPKM -a GSE52968/bedgraph/SRR1042853/final.psite -b annotation_yeast/final.ORFs -o GSE52968 -n SRR1042853 -s scripts -p 8;
```

##### Calculating frameshift potential 

```
scripts/main_function.sh -F -a GSE52968/bedgraph/SRR1042853/final.psite -b annotation_yeast/final.ORFs -o GSE52968 -n SRR1042853 -s scripts -p 8;
```

##### Multiple functions 

```
scripts/Ribowave -PDF -T 12378563 GSE52968/mRNA/SRR1042851.FPKM -a GSE52968/bedgraph/SRR1042853/final.psite -b annotation_yeast/final.ORFs -o GSE52968/Ribowave   -n SRR1042853 -s scripts -p 8;
```

#### Input files:

- <P-site track\> : output from the previous step, containing the P-site track of transcripts of interest

- <ORF_list> : ORFs of interest ,eg : `final.ORFs`

- <Ribo-seq total reads\> : the total number of Ribo-seq reads to calculate FPKM , eg: `12378563`

- <RNA FPKM\> : FPKM table. It may look like this :

	```
	YAL067W-A       0
	YAL064C-A       0.834264
	YAL066W 	0.452034
	YAL067C 	3.0719
	YAL064W-B       5.00558
	```

#### Output files:

* _name_.feats1 	: the features of ORFs including chi-square P-value information. It may look like this :

Column | Explanation	
------------ | -------------
Column1-Column4 | Basic information about the ORF
Column5		| Reads coverage within the ORF
Column6		| P-value predicted by RiboWave
Column7		| Translational signal outside current ORF
Column8		| Reads intensity at start codon

**`result`** directory, including :

* _name_.95%.mx 	: the final translated product of RiboWave with translation initiation sites specified ( **p.value < 0.05** ) . It may look like this :

```
YBR073W_0_103_2874
YBR152W_0_406_873
YBR197C_0_292_651
```

* _name_.density	: reads density ( **PF P-site** ) of given ORFs. It may look like this : 

Column | Explanation	
------------ | -------------
Column1-Column4 | Basic information about the ORF
Column5		| number of PF P-sites in transcript
Column6		| number of PF P-sites in given ORF
Column7		| Density of PF P-sites in given ORF


* _name_.TE   		: TE of given ORFs. It may look like this :

Column | Explanation	
------------ | -------------
Column1		| transcript 
Column2		| ORF 
Column3		| TE


* _name_.CRF.final  	: ORFs that might experience reading frame translocation. It may look like this :

Column | Explanation	
------------ | -------------
Column1 | ORF
Column2		| Start of gap region 
Column3		| Stop of gap region
Column4		| Reading frames after the shift ,eg: `0_7,2_2` where `0_7` refers to continuous seven PF P-sites of frame 0 followed by continuous two PF P-sites of frame 2.
Column5		| Corresponding PF P-sites position after the shift ,eg : `1996,2152,2197,2437,2569,2596,2605;2712,2907` where `1996,2152,2197,2437,2569,2596,2605` corresponds to the exact position of `0_7` within the transcript. Discontinuity in the reading frame is separated by `;`
Column6		| F1 score describing the extent of change in the reading frame before and after the gap
Column7 	| F2 score describing the extent of reading frame consistency before and after the gap
Column8		| Fscore describing the potential of frameshift 


## Example file
An example file is packed and found in [RiboWave](http://lulab.life.tsinghua.edu.cn/RiboWave/RiboWave_v1.0.tar.gz).

Enclosed in the RiboWave_v1.0.tar.gz, `run_Ribowave_yeast.sh` combines all the steps together into a pipeline. 

