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

#### Arguments:
- annotation_dir  : 
  - 1.annotation.gtf : the annotation gtf should contain ***start_codon*** and ***stop_codon*** information  `Saccharomyces_cerevisiae.R64-1-1.90.gtf` 
  - 2.genome.fasta `Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa` 

- scripts_dir 	: the directory of all the scripts in the package

### create annotation

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

#### Arguments:
- Ribo_bam 	: **secondary alignment removed** and sorted

- annotation_dir  : 
  - 1.annotation directory with all ORFs scanned and annotated `final.ORFs` 
  - 2.annotated start site `start_codon.bed` 
  - 3.genome size `genome` 
  - 4.exon annotation gtf `exons.gtf` 

- out_dir 	: the directory of the output result, eg: `GSE52968`

- output_header 	: the header of all the output file, eg: `SRR1042853` 

- scripts_dir 	: the directory of all the scripts in the package


### Determine P-site 

This step determine the P-site position for each read length by overlapping with the annotated start codon 

```
Usage: ./P-site_determination.sh <Ribo_bam> <start_codon.bed> <out_dir> <output_identifier> <scripts_dir>

Example: scripts/P-site_determination.sh   GSE52968/SRR1042853.sort.bam       annotation_yeast/start_codon.bed      GSE52968        SRR1042853         scripts;
```

#### Output files:
**`P-site`** directory, including :

* _identifier_.psite1nt.txt 	: the P-sites position for each length

* _identifier_.psite.pdf 	: the pdf displaying the histogram of aggregated reads


### generating P-site track 

This step creats the P-site track for transcripts of interests

```
Usage: ./create_track_Ribo.sh <Ribo_bam> <transcripts.exon.gtf> <genome> <out_dir> <output_identifier> <scripts_dir>

Example: scripts/create_track_Ribo.sh      GSE52968/SRR1042853.sort.bam              annotation_yeast/exons.gtf       annotation_yeast/genome      GSE52968     SRR1042853         scripts;
```

#### Output files:

**`bedgraph`** directory, including :

* final.psite 	: P-site track for each interested transcript


### Predict ORF translation

This step predicts the translated ORF

```
Usage: ./create_track_Ribo.sh <transcript exon gtf> <genome> <out_dir> <output_header> <scripts_dir> <cores>

Example: scripts/main_function.sh          annotation_yeast/final.ORFs     GSE52968         SRR1042853         scripts     8;
```

#### Output files:

* _identifier_.feats1 	: the features of ORFs including chi-square P-value information

* _identifier_.COV	: reads coverage of ORFs


### Identify translated ORF

This step incorporates all the information from each ORF and find ORFs that are predicted to be translated ( P-value < 0.05) 

```
Usage: ./translated_protein_annotation.sh <annotation_dir> <out_dir> <output_header> <scripts_dir>

Example: scripts/translated_protein_annotation.sh  annotation_yeast       GSE52968       SRR1042853    scripts;
```

#### Output files:

* _identifier_.mx 			: the combined information of all ORFs including chi-square P-value information and coverage information

**`result`** directory, including :

* _identifier_.95%.mx 	: protein products predicted to be translated in the sample within the cutoff of P-value < 0.05

* _identifier_.95%.ORF_category : annotate ORFs in $header.95%.mx by the relative position of the annotated ORF and customize the output

