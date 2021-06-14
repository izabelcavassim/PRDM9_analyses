Amino acid conservation analyses
============

We carried out a residue conservation analysis using an approach proposed by (Capra and Singh 2007), using the code score_conservation.py available at https://compbio.cs.princeton.edu/conservation/. This approach quantifies the Jensen–Shannon divergence between the amino acid distribution of the focal residue and a “background amino acid distribution.” The alignment of ZCWPW2 was produced using Clustal Omega (using default parameters) within MEGA (version 7, (Kumar et al. 2017; Kumar, Stecher, and Tamura 2016)). As recommended, the overall background amino acid distribution was drawn based on the BLOSUM62 amino acid substitution matrix provided by the software (Capra and Singh 2007). Any column of the gene sequence alignment with more than 30% gaps was ignored. A window size of 3 was used to incorporate information from sequential amino acids, as recommended by the default settings.



```{bash}
python score_conservation.py -a Homo_sapiens_ZCWPW2 -o Conservation_results.txt complete.ZCWPW2_alignment_Zach.fas
```

The file complete.ZCWPW2_alignment_Zach.fas is found in the data folder. 
The python script score_conservation.py, is found in the scripts folder.

