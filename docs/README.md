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
