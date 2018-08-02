# RiboWave 

RiboWave is a funtional Ribo-seq analysis tool to identify translated ORF based on Ribo-seq data.


The RiboWave workflow consists of:


- Pre-processing :
  - Create the annotation file for the subsequent analysis. [`create_annotation.sh`]
  - Determine the P-site position of Ribo-seq. [`P-site_determination.sh`]
  - Generate P-site tracks from Ribo-seq.[`create_track_Ribo.sh`]

- Main function :
  - Denoise [`Ribowave`]
  - Predict translated ORFs [`Ribowave`]
  - Estimate reads density for each given ORF [`Ribowave`]
  - Estimate frameshift potential for each given ORF [`Ribowave`]


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

This step scans for and annotates all putative ORFs as well as define annotated ORF according to the annotation file

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
scripts/create_annotation.sh -G annotation_fly/dmel-all-r6.18.gtf -f annotation_fly/dmel-all-chromosome-r6.18.fasta  -o annotation_fly  -s scripts;
```

#### Input files:
- <annotation.gtf>      : the annotation gtf should contain **start_codon** and **stop_codon** information,eg: `dmel-all-r6.18.gtf` 
- <genome.fasta>        : genome fasta ,eg: `dmel-all-chromosome-r6.18.fasta` 

- <annotation_dir>      :  the directory for all the annotation output

- <scripts_dir> 	: the directory of all the scripts in the package

#### Output files:

**`annotation`** directory, including :

* start_codon.bed 	: the bed file annotating start codon 

* final.ORFs 	: all identified ORFs, eg: `FBtr0300105_0_31_546` where `FBtr0300105` refers to the transcript, `0` refers to the reading frame relative to the start of transcript, `31` refers to the start site, `546` refers to the stop codon.  

### 1. P-site determination

This step determines the P-site position for each Ribo-seq reads length by overlapping with the annotated start codons from previous step

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
scripts/P-site_determination.sh -i GSE52799/SRR1039770.sort.bam -S annotation_fly/start_codon.bed -o GSE52799 -n SRR1039770 -s scripts;
```

#### Input files:
- <Ribo_bam> :  **secondary alignment removed**  to ensure one genomic position per aligned read and sorted

- **`annotation`**  : 
  - <start_codon.bed> : annotated start site `start_codon.bed`. It is generated in the `create_annotation.sh` step.

- <out_dir> 	: the directory of the output result, eg: `GSE52799`

- <study_name> 	: the name of all the output file, default: test. eg: `SRR1039770` 

- <scripts_dir>	: the directory of all the scripts in the package


#### Output files:
**`P-site`** directory, including :

* _name_.psite1nt.txt 	: the Ribo-seq reads length and its corresponding P-sites position(= offset + 1). It may look this this : 
  
  ```
  30	13
  ```
  
* _name_.psite.pdf 	: the **PDF** displaying the histogram of aggregated reads


### 2. Generating P-site track 

This step creats the P-site track for transcripts of interests using determined P-sites position from previous step.

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
	-i	<filename>	(Ribo-seq bam                                      )
	-G	<filename>	(exons gtf file                                    )
	-g	<filename>	(Genome size annotation file                       )
	-P	<filename>	(P-site position (offset+1)                        )
	-o	<directory>	(Output annotation directory                       )
	-s	<directory>	(Script directory                                  )
	-n	<string>  	(study name, defult:test                           )
	-h	          	(Help                                              )
----------------------------------------------------------------------------------------------------
```

Run `create_track_Ribo.sh` on example:
###### look at transcripts from chromosome X :

```
scripts/create_track_Ribo.sh -i GSE52799/SRR1039770.sort.bam -G annotation_fly/X.exons.gtf -g annotation_fly/genome -P GSE52799/P-site/SRR1039770.psite1nt.txt -o GSE52799 -n SRR1039770 -s scripts;
```

#### Input files:

- <Ribo_bam> 

- <exons.gtf> :  a gtf file for only the exons from transcripts of intersect, eg: `X.exons.gtf`

- <genome\> :  the file including all the chromosomes and its genome size.  **Noted: `genome` can be obtained by using `samtools faidx` function with the input of fasta file.** `genome` may look like this:
    
    ```
    2L	23513712
    2R	25286936
    3L	28110227
    3R	32079331 
    ```

- **`P-site`**  : 
	- <P-site_position> :  the file listing the P-site position for each read length. This file can be found in the output of previous step, eg: _name_.psite1nt.txt

- <out_dir> 	: the directory of the output result, eg: `GSE52799`

- <study_name> 	: the name of all the output file, default: test. eg: `SRR1039770` 

- <scripts_dir> : the directory of all the scripts in the package

#### Output files:

**`bedgraph/name`** directory, including :

* final.psite 	: P-site track at transcriptome wide. It may look like this : 
  
  ```
  FBtr0070533	0,0,0,0,0,0,0,0,0,0,0,0,6,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,1,0,0,0,0,1,0,2,1,0,0,0,0,0,0,4,8,0,0,3,0,5,12,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
  FBtr0073886	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,25,0,0,0,0,0,0,0,0,0,0,0,0,0
  FBtr0070604	0,0,0,0,0,0,0,0,0,0,0,0,59,6,0,1,0,0,2,6,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
  FBtr0070603	0,0,0,0,0,0,0,0,0,0,0,0,75,2,7,10,7,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,3,3,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
  ```
  
## Main function

### 3. RiboWave

This step can achieve multiple functions : 

  - denoising [`denoise`]

  - providing predicted p.value for each given ORF to identify its translation status [`pvalue`,`-P`]
  
  - providing reads density (P-site/PF P-site) for each given ORF [`density`,`-D`]
  
  - providing translation efficiency (TE) estimation for each given ORF [`TE`,`-T`]
  
  - providing frameshift potential (CRF score) for each given ORF [`CRF`,`-F`]
  

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
	-T	<int>  <RNA_FPKM>   	(estimating TE for each ORF.<int> requires total Ribo-seq sequence depth;<RNA_FPKM> requires RNA-seq FPKM)
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
mkdir -p Ribowave;
scripts/Ribowave  -a GSE52799/bedgraph/SRR1039770/final.psite -b annotation_fly/final.ORFs -o GSE52799/Ribowave -n SRR1039770 -s scripts -p 8;
```

##### Identifying translated ORF

```
mkdir -p Ribowave;
scripts/Ribowave -P -a GSE52799/bedgraph/SRR1039770/final.psite -b annotation_fly/final.ORFs -o GSE52799/Ribowave -n SRR1039770 -s scripts -p 8;
```

##### Estimating abundance 

```
mkdir -p Ribowave;
scripts/Ribowave -D -a GSE52799/bedgraph/SRR1039770/final.psite -b annotation_fly/final.ORFs -o GSE52799/Ribowave -n SRR1039770 -s scripts -p 8;
```

##### Estimating TE
IMPORTANT :
when estimating TE, user should input the **sequenced depth of Ribo-seq** and the **FPKM value from paired RNA-seq**

```
mkdir -p Ribowave;
scripts/Ribowave -T 9012445  GSE52799/mRNA/SRR1039761.RPKM -a GSE52799/bedgraph/SRR1039770/final.psite -b annotation_fly/final.ORFs -o GSE52799/Ribowave -n SRR1039770 -s scripts -p 8;
```

##### Calculating frameshift potential 
###### on annotated ORFs

```
mkdir -p Ribowave;
awk -F '\t' '$3=="anno"'  annotation_fly/final.ORFs  >   annotation_fly/aORF.ORFs;
scripts/Ribowave -F -a GSE52799/bedgraph/SRR1039770/final.psite -b annotation_fly/aORF.ORFs -o GSE52799/Ribowave -n SRR1039770 -s scripts -p 8;
```

##### Multiple functions 

```
mkdir -p Ribowave;
scripts/Ribowave -PD -T 9012445  GSE52799/mRNA/SRR1039761.RPKM -a GSE52799/bedgraph/SRR1039770/final.psite -b annotation_fly/final.ORFs -o GSE52799/Ribowave -n SRR1039770 -s scripts -p 8;
```

#### Input files:

- **`bedgraph/name`**:
	- <P-site track\> : output from the previous step, containing the P-site track of transcripts of interest, eg: `final.psite`

- <ORF_list> : ORFs of interest ,eg : `final.ORFs`. It is generated in the step of `create_annotation.sh`

- <Ribo-seq sequenced depth\> : the sequenced depth of Ribo-seq to calculate FPKM , eg: `9012445`

- <RNA FPKM\> : FPKM table. It may look like this :

	```
	FBtr0100871	22262
	FBtr0070604	18682
	FBtr0100231	14746.5
	FBtr0100874	14024.5
	FBtr0100864	11475.6
	```

#### Output files:
* _name_.PF_psite	: the denoised signal track(PF P-sites signal track) at transcriptome wide. It looks similar as the input `final psite`.

* _name_.feats1 	: the features of ORFs including chi-square P-value information. It may look like this :

Column | Explanation	
------------ | -------------
Column1-Column4 | Basic information about the ORF
Column5		| Reads coverage within the ORF
Column6		| P-value predicted by RiboWave
Column7		| Values estimating the relative abundance of PF P-sites outside of the studied ORF
Column8		| Reads intensity at the current start codon

**`result`** directory, including :

* _name_.95%.mx 	: RiboWave makes the prediction on the translation initiation sites and gives the final translated product output (**p.value < 0.05**) . It may look like this :

```
FBtr0070007_2_93_1028
FBtr0070008_1_128_943
FBtr0070025_2_135_1094
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
Column2		| Start of frameshift  
Column3		| Stop of frameshift
Column4		| PF P-sites' reading frames **after** the change point ,eg: `2_2,0_1` where `2_2` refers to continuous two PF P-sites of frame 2 followed by continuous one PF P-sites of frame 0.
Column5		| Relative position of PF P-sites **after** the shift ,eg : `1413,1440;1789` where `1413,1440` corresponds to the exact position of `2_2` within the transcript. Discontinuity in the reading frame is separated by `;`
Column6		| CRF score describing the potential of frameshift 


## Example file
An example file is packed and found in [here](http://lulab.life.tsinghua.edu.cn/RiboWave/RiboWave_v1.0.tar.gz).

Enclosed in the RiboWave_v1.0.tar.gz, `run_Ribowave_dmel.sh` combines all the steps together into a pipeline. 

## For any questions, please contact:

Zhiyu Xu [xanthexu18@gmail.com]
