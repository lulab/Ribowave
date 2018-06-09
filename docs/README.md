# Ribo-seq data analysis tool 

RiboWave analyses Ribosome profiling data (Ribo-seq). It utilizes wavelet transform to denoise the original signal by extracting 3-nt periodicity of ribosomes (i.e. signal frequency) and precisely locate their footprint. 

>Translation is dynamically regulated during cell development and stress response. In order to detect actively translated open reading frame (ORF) and dynamic cellular translation events, we have developed a computational method, RiboWave, to process ribosome profiling data. RiboWave utilizes wavelet transform to denoise the original signal by extracting 3-nt periodicity of ribosomes and precisely locate their footprint denoted as Periodic Footprint P-site (PF P-site). Such high-resolution footprint is found to capture the full track of actively elongating ribosomes, from which translational landscape can be explicitly characterized. 

>We show that RiboWave outperforms other tools in terms of both accuracy and usage when defining actively translating ORFs. Moreover, we show that PF P-site derived by RiboWave shows superior performance in characterizing the dynamics and complexity of cellular translatome by accurately estimating the abundance of protein levels, assessing differential translation and identifying dynamic translation frameshift. 



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

#### More detail can be inferred from [RiboWave Github](https://github.com/lulab/Ribowave).
