Introduction
============

This repository was created to outline the analyses described in the paper Cavassim et al., 2020 of which we extracted all the RefSeq data available and searched for genes coevolving with PRDM9. Multiple steps until the final results were executed and will be described in their chronological order here. 
Main analyzes are: 
0. Building confident calls for PRDM9
1. [RNAseq_analyses.md](./RNAseq_analyses.md) Contains scripts used to verify PRDM9 calls using RNA-seq data.
2. [Orthologous_search.md](./Orthologous_search.md) Contains scripts used to identify orthologs of candidate co-evolving meiosis genes across vertebrates (339 species). 
3. [PIC.Rmd](./PIC.Rmd) Contains R script used for the phylogenetic independent contrats (PIC) between PRDM9 and candidate genes.

* The folder *data* includes some of the data used for the PIC analyses. Other files are found the following repository https://doi.org/10.6084/m9.figshare.11672685

* RNAseq data of 2 fish species and 2 reptiles (fastq files):
Is Publicly available at the NCBI BioProject (accession no: PRJNA605699).

For more information or data access please do not hesitate to contact me! izabelcavassim at g mail . com
