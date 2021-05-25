Introduction
============

This repository was created to outline the analyses described in the paper Cavassim et al., 2020 of which we extracted all the RefSeq data available for vertebrates (July 2020) and searched for genes coevolving with PRDM9. Multiple steps until the final results were executed and will be described in their chronological order here. 
Main analyzes are: 

1. Building confident calls for PRDM9. The details are described in [Baker et al., 2017](https://elifesciences.org/articles/24133)
2. [RNAseq_analyses.md](./RNAseq_analyses.md) Contains scripts used to verify PRDM9 calls using RNA-seq data.
3. [Orthologous_search.md](./Orthologous_search.md) Contains scripts used to identify orthologs of candidate co-evolving meiosis genes across vertebrates (339 species). 
4. [PIC.Rmd](./PIC.Rmd) Contains R script used for the Phylogenetic Independent Contrats (PIC) between PRDM9 and candidate genes.

* The folder *data* includes some of the data used or generated for the PIC analyses. Other files are found the following repository https://doi.org/10.6084/m9.figshare.11672685

* RNAseq data of 2 fish species and 2 reptiles (fastq files):
Is Publicly available at NCBI under the [NCBI BioProject: PRJNA605699](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA605699).

For more information or data access please do not hesitate to contact me! izabelcavassim at g mail . com
