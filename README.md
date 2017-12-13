# RiboWave 

RiboWave is a funtional Ribo-seq analysis tool to identify translated ORF based on Ribo-seq data.


The RiboWave workflow consists of:

* Create the annotation file for the subsequent analysis. [`create_annotation.sh`]

* Determine the P-site position of Ribo-seq. [`P-site_determination.sh`]

* Generate P-site tracks from Ribo-seq, dependent on P-sites calculation cutoffs.[`create_track_Ribo.sh`]

* Predicting translating ORFs [`main_function.sh`]

* Predict translated ORF and annotating translatome [`translated_protein_annotation.sh`]

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

## Before running 
It is **recommanded** to make a new directory and move the Ribo-seq bam file into that directory;

#### Input files:
- annotation_dir  : 
  - 1.annotation.gtf : the annotation gtf should contain ***start_codon*** and ***stop_codon*** information  `Saccharomyces_cerevisiae.R64-1-1.90.gtf` 
  - 2.genome.fasta `Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa` 

- scripts_dir 	: the directory of all the scripts in the package

### 0. Create annotation

This step scans for and annotates all putative ORFs 

```
Usage: ./create_annotation.sh <annotation_dir> <genome.gtf> <fasta> <scripts_dir>

Example: scripts/create_annotation.sh   annotation_yeast     annotation_yeast/Saccharomyces_cerevisiae.R64-1-1.90.gtf    annotation_yeast/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa   scripts;
```

#### Output files:
**`annotation`** directory, including :

* start_codon.bed 	: annotated start codon 

* final.ORFs 	: all ORFs 

## Workflow

#### Input files:
- Ribo_bam 	: **secondary alignment removed** and sorted

- annotation_dir  : 
  - 1.annotation directory with all ORFs scanned and annotated `final.ORFs` 
  - 2.annotated start site `start_codon.bed` 
  - 3.genome size `genome` looks like this:
    
    ```
    I 230218
    II 813184
    III 316620
    IV 1531933
    ```
    
  - 4.exon annotation gtf `exons.gtf` 

- out_dir 	: the directory of the output result, eg: `GSE52968`

- output_header 	: the header of all the output file, eg: `SRR1042853` 

- scripts_dir 	: the directory of all the scripts in the package


### 1. Determine P-site 

This step determine the P-site position for each read length by overlapping with the annotated start codon 

```
Usage: ./P-site_determination.sh <Ribo_bam> <start_codon.bed> <out_dir> <output_identifier> <scripts_dir>

Example: scripts/P-site_determination.sh   GSE52968/SRR1042853.sort.bam       annotation_yeast/start_codon.bed      GSE52968        SRR1042853         scripts;
```

#### Output files:
**`P-site`** directory, including :

* _identifier_.psite1nt.txt 	: the P-sites position (= offset + 1) for each length. It may look this this : 
  
  ```
  25	10
  26	11
  27	12
  28	13
  29	11
  30	12
  ```
  
* _identifier_.psite.pdf 	: the pdf displaying the histogram of aggregated reads


### 2. Generating P-site track 

This step creats the P-site track for transcripts of interests

```
Usage: ./create_track_Ribo.sh <Ribo_bam> <transcripts.exon.gtf> <genome> <out_dir> <output_identifier> <scripts_dir>

Example1: scripts/create_track_Ribo.sh      GSE52968/SRR1042853.sort.bam              annotation_yeast/exons.gtf       annotation_yeast/genome      GSE52968     SRR1042853         scripts;

Example2: scripts/create_track_Ribo.sh      GSE52799/SRR1039770.sort.bam       annotation_fly/chrX.exons.gtf       annotation_fly/genome  GSE52799     SRR1039770         scripts;
```

#### Output files:

**`bedgraph`** directory, including :

* final.psite 	: P-site track for each interested transcript. It may look like this : 
  
  ```
  YAL044C 3,0,0,3,10,0,0,0,0,2,0,0,36,0,0,4,0,0,0,6,0,12,0,7,32,0,0,6,7,9,19,2,5,28,0,0,0,0,0,0,0,0,0,0,4,0,0,0,24,0,1,34,0,1,9,2,0,8,0,0,0,0,0,38,0,4,33,0,10,24,0,8,2,0,6,16,0,0,2,0,0,4,0,0,0,2,0,
  YAL045C 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
  YAL046C 3,0,0,2,0,0,4,0,0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,1,0,0,16,0,0,6,0,0,0,0,0,12,0,0,19,0,0,0,0,0,5,0,1,5,0,0,0,0,1,6,0,0,0,0,0,14,0,0,2,0,0,7,0
  ```

### 3. RiboWave main function

This step predicts the translated ORF

```
Usage: ./create_track_Ribo.sh <transcript exon gtf> <genome> <out_dir> <output_header> <scripts_dir> <cores>

Example: scripts/main_function.sh          annotation_yeast/final.ORFs     GSE52968         SRR1042853         scripts     8;
```

#### Output files:

* _identifier_.feats1 	: the features of ORFs including chi-square P-value information

* _identifier_.COV	: reads coverage of ORFs


### 4. Identify translated ORF

This step incorporates all the information from each ORF and find ORFs that are predicted to be translated ( P-value < 0.05) 

```
Usage: ./translated_protein_annotation.sh <annotation_dir> <out_dir> <output_header> <scripts_dir>

Example: scripts/translated_protein_annotation.sh  annotation_yeast       GSE52968       SRR1042853    scripts;
```

#### Output files:

* _identifier_.mx 			: the combined information of all ORFs including chi-square P-value information and coverage information. It may look this this : 

  ```
  ORFID	TRtype	label	codon	getorf	orflonger	TRID	Start	Stop	psites	Start	Covfrm	CovORF	pvalue	upstream_3nt	startPsite
  YAL001C_0_1_3480	coding	unanno	ATG	1	3480	YAL001C	1	3480	67	0	0.0482758620689655	0.0192528735632184	1	0	0
  YAL001C_0_382_3480	coding	unanno	ATG	382	3480	YAL001C	382	3480	67	0	0.0435624394966118	0.0167796063246208	1	0.0207452938916635	0
  YAL001C_0_943_3480	coding	unanno	ATG	943	3480	YAL001C	943	3480	67	0	0.0378250591016548	0.0149724192277384	1	0.119861698040722	0
  ```
  
  

**`result`** directory, including :

* _identifier_.95%.mx 	: protein products predicted to be translated in the sample within the cutoff of P-value < 0.05

* _identifier_.95%.ORF_category : annotate ORFs in $header.95%.mx by the relative position of the annotated ORF and customize the output. It may look like this:

  ```
  FBtr0070007_2_93_1028	NC	anno	ATG	FBtr0070007_2_93_1028	anno
  FBtr0070052_2_255_905	NC	sharestop	ATG	FBtr0070052_2_57_905	truncated
  FBtr0071362_2_330_1952	NC	sharestop	ATG	FBtr0071362_2_609_1952	extended
  FBtr0070151_0_247_294	NC	unanno	ATG	FBtr0070151_1_668_4567	uORF
  FBtr0070379_2_3693_3761	NC	unanno	ATG	FBtr0070379_0_859_3693	internal_overlapped
  FBtr0347011_0_667_2583	NC	unanno	ATG	NA	unanno
  ```
  

