Introduction
============

This repository was created to outline the analyses described in the paper Cavassim et al., 2020 of which we extracted all the RefSeq data available and searched for genes coevolving with PRDM9.

Retrieving 100-way alignment sequences from UCSC
-----------------------------
We used the 100 vertebrate species amino acid alignments (100way project) from the multiz alignment available at the UCSC genome browser (Blanchette et al., 2004; Harris, 2007). The description of how the orthology was conducted for the 100-way alignment project is found [here](http://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg19&g=cons100way).
We imposed that exons with with a minimum threshold of sequence completeness of 20%. From the UCSC database we have a total of **201,663** exons (of the longest isoform).

The files "knownCanonical.exonAA.fa.gz" and "knownCanonical.exonNuc" were downloaded from [UCSC repository](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/multiz100way/alignments/) and exons were split in singular fasta files using a bash script.

``` bash
#!/bin/bash
i=1
fileName="gene_$i.fasta"
while read line ; do 
if [ "$line"  == ""  ] ; then
 ((++i))
 fileName="gene_$i.fasta"
else
 echo $line >> "$fileName"
fi
done < knownCanonical.exonAA.fa
``` 
Exons were further parsed in order to have headers matching with species tree headers (for that we used a metadata file metadata2.csv). Exons belonging to the same gene were combined into a singular fasta file using the script parsing_UCSC_data.py. This gave us a total of **21,678** genes to work with.

Downloading the Refseq protein sequences from NCBI
-----------------------------
First downloaded the species present in each of the vertabrates directories:

``` bash
# Mammalian vertebrate species
cat Extra_genomes_02_06_20_mammalian_vertebrates.txt | while read p;
  do   
    lftp -e "cls -1 > /tmp/list; exit" ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/"${p}"
    read -r  firstline< /tmp/list
    echo $firstline
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/"${p}""${firstline::-1}"/"${firstline::-1}_protein.faa.gz" 
    mv ${firstline::-1}_protein.faa.gz ${p%%/*}_protein.faa.gz
done

# Other vertebrate species
cat Extra_genomes_02_06_20_other_vertebrates.txt | while read p;
  do   
  lftp -e "cls -1 > /tmp/list; exit" ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/"${p}"
  read -r  firstline< /tmp/list
  echo $firstline
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/"${p}""${firstline::-1}"/"${firstline::-1}_protein.faa.gz"
  mv ${firstline::-1}_protein.faa.gz ${p%%/*}_protein.faa.gz  
done
```
