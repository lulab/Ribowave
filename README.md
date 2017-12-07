# RiboWave 

RiboWave is a funtional Ribo-seq analysis tools to identify translated ORF based on Ribo-seq data.

The RiboWave workflow consists of:

* Create the annotation file for the subsequent analysis. [`create_annotation.sh`]

* Determine the P-site position of Ribo-seq. [`create_track_Ribo.sh`]

* Generate P-site tracks from Ribo-seq, dependent on P-sites calculation cutoffs.[`create_track_Ribo.sh`]

* Predicting translating ORFs [`main_function.sh`]

* Predict the final protein product based on results from step3 [`translated_protein_annotation.sh`]

* Customize the output by annotating the relative position of translating ORFs. [`translated_protein_annotation.sh`]


## Requirements
### software
R 

bedtools v2.25.0 
### R packages
reshape

ggplot2

rhdf5

methods

wmtsa

## Before running 

It is **recommanded** to make a new directory and move the Ribo-seq bam file into that directory;


## workflow

#### Arguments:

- Ribo_bam 	: alignments of Ribo-seq followed by elimination of alignments flagged as secondary alignments to ensure one genomic position for one single 

- annotation_dir  : 
  - 1.annotation directory enclosed in the (`.zip`) file with all ORFs scanned `final.ORFs.annot3` 
  - 2.annotated start site `start_codon.bed` 
  - 3.genome size `genome` 
  - 4.exon annotation gtf `exons.gtf`

- out_dir 	: the directory of the output result, *eg: HCT116_test*

- output_header 	: the header of all the output file, *eg: HCT116* 

- scripts_dir 	: the directory of all the scripts in the package( enclosed in the `.zip` file)


### Determine P-site and create P-site track on transcript level 

```
Usage: ./create_track_Ribo.sh <Ribo_bam> <annotation_dir> <out_dir> <output_header> <scripts_dir>
```

This step determine the P-site position for each read length by overlapping with the annotated start codon and generate P-site track for each transcript.


#### Output files:

*$out_dir/P-site/$header.psite1nt.txt* 	: the P-sites position for each length

*$out_dir/P-site/$header.psite.pdf* 	: the pdf displaying the histogram of aggregated reads

*$out_dir/bedgraph/$header/final.psite* 	: the created P-sites track of each transcripts 


### Predict ORF translation

```
Usage: ./main_function.bash <annotation_dir> <out_dir> <output_header> <scripts_dir>
```

This step takes the information from the P-site track for each transcript and predict the translation status for each ORF.

#### Output files:

*$out_dir/out/$header.TR.psites.000* 	: the features of ORFs including chi-square P-value information

*$out_dir/out/$header.TR.psites.000.cv* 	: reads coverage of ORFs


### Identify translated ORF

```
Usage: ./translated_protein_annotation.sh <annotation_dir> <out_dir> <output_header> <scripts_dir>
```

This step incorporates all the information from each ORF and find ORFs that are predicted to be translated ( P-value < 0.05) 

#### Output files:

*$out_dir/$header.mx* 			: the combined information of all ORFs including chi-square P-value information and coverage information

*$out_dir/result/ATG/$header.95%.mx* 	: protein products predicted to be translated in the sample within the cutoff of P-value < 0.05

*$out_dir/result/ATG/$header.99%.mx* 	: protein products predicted to be translated in the sample within the cutoff of P-value < 0.01

*$out_dir/result/ATG/$header.95.ORF_category* : annotate ORFs in $header.95%.mx by the relative position of the annotated ORF and customize the output

*$out_dir/result/ATG/$header.99.ORF_category* : annotate ORFs in $header.99%.mx by the relative position of the annotated ORF and customize the output

