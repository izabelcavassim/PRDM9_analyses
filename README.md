Introduction
============

This repository was created to outline the analyses described in the paper Cavassim et al., 2020 of which we extracted all the RefSeq data available and searched for genes coevolving with PRDM9.

Retrieving 100-way alignment sequences from UCSC
-----------------------------
We used the 100 vertebrate species amino acid alignments (100way project) from the multiz alignment available at the UCSC genome browser (Blanchette et al., 2004; Harris, 2007), gene alignments with a minimum number of species (completeness) will be defined for this study. 

From the database we have a total of **201,663** exons (201663) (of the longest isoform).

The file "knownCanonical.exonAA.fa.gz" was retrieved from [UCSC](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/multiz100way/alignments/) and exons were split in singular fasta files using a bash script.
[description](http://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg19&g=cons100way)
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
