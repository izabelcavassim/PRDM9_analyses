Introduction
============

This repository was created to outline the analyses described in the paper Cavassim et al., 2020 of which we extracted all the RefSeq data available for vertebrates (July 2020) and searched for genes coevolving with PRDM9. Multiple steps until the final results were executed and will be described in their chronological order here. 
Main analyzes are: 

1. Building confident calls for PRDM9. The details are described in [Baker et al., 2017](https://elifesciences.org/articles/24133)
2. [RNAseq_analyses](./RNAseq_analyses.md) contains the script used to verify PRDM9 and candidate gene calls using RNAseq data.
3. [Orthologous search](./Orthologous_search.md) contains scripts used to identify orthologs of candidate co-evolving meiosis genes across vertebrates (339 species). 
4. [Phylogenetic tests](./PIC.Rmd) contains the R script used for the phylogenetic tests between two binary traits: PRDM9 completeness and the presence/absence of candidate genes (139).
5. [Conservation analyses](./Conservation_analyses.md) contains the scripts and files for the ZCWPW2 conservation of residues analyses.
6. [Codeml analyses](./scripts/PAML_codeml_analyses.py) contains the script for the dN/dS analyses. 


* The folder *Figures* includes the data and scripts used for reproducing each main figure of the paper.

* The folder *data* includes some of the data used or generated for the PIC analyses. Other files are found the following repository https://doi.org/10.6084/m9.figshare.11672685

* RNAseq data of 2 fish species and 2 reptiles (fastq files):
Is Publicly available at NCBI under the [NCBI BioProject: PRJNA605699](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA605699).

For more information or data access please do not hesitate to contact me! izabelcavassim at g mail . com
