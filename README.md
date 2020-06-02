Introduction
============

This repository was created to outline the analyses described in the paper Cavassim et al., 2020 of which we extracted all the RefSeq data available and searched for genes coevolving with PRDM9.


Downloading the Refseq protein sequences
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
cat Extra_genomes_02_06_20_mammalian_vertebrates.txt | while read p;
  do   
  lftp -e "cls -1 > /tmp/list; exit" ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/"${p}"
  read -r  firstline< /tmp/list
  echo $firstline
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/"${p}""${firstline::-1}"/"${firstline::-1}_protein.faa.gz"
  mv ${firstline::-1}_protein.faa.gz ${p%%/*}_protein.faa.gz  
done
```
