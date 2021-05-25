Introduction
============

This document was created as a repository for the scripts used in the
analysis of the PRDM9 paper (link to the preprint XXXXXXXXX).

Data availability
===================
RNAseq data of 2 fish species and 2 reptiles (fastq files) will be publicly available at the NCBI BioProject (accession no: PRJNA605699).


Demultiplexing sequence reads
=============================

``` bash
#!/bin/sh                                                                                                 
#SBATCH --account=palab # The account name for the job.                                                   
#SBATCH --job-name=Demultiplexing  # The job name.                                                        
#SBATCH -c 1     # The number of cpu cores to use.                                           
#SBATCH -t 2:00:00      # The time the job will take to run                                              
#SBATCH --mem-per-cpu 1gb                                                                                  
#SBATCH -o /moto/palab/users/mc4719/RNA_seq_pipeline/data/Demultiplexed_PRDM9_samples/bcl2fastq_Run1.out  
#SBATCH -e /moto/palab/users/mc4719/RNA_seq_pipeline/data/Demultiplexed_PRDM9_samples/bcl2fastq_Run1.err                                                                     
export PATH=$PATH:/mc4719/opt/gcc-7.3.0/lib64/
export LD_LIBRARY_PATH=/mc4719/opt/gcc-7.3.0/lib64/
export PATH=/moto/palab/users/mc4719/RNA_seq_pipeline/bcl2fastq2-v2.20.0/bin:$PATH

# First run                                                                                                                                                                  
bcl2fastq --runfolder-dir /moto/palab/users/crh2152/projects/PRDM9_RNAseq/190517_NB551405_0099_AHKCG5AFXY --output-dir /moto/palab/users/mc4719/RNA_seq_pipeline/data/Demultiplexed_PRDM9_samples --sample-sheet /moto/palab/users/mc4719/RNA_seq_pipeline/data/indices_RNAseq_libraries_demultiplex.csv                                                
# Second run                                                                                                                                                                 
bcl2fastq --runfolder-dir /moto/palab/users/crh2152/projects/PRDM9_RNAseq/05_28_2017_CRH/190528_NB551405_0102_AHKKKFAFXY --output-dir /moto/palab/users/mc4719/RNA_seq_pipeline/data/Demultiplexed_PRDM9_samples_2 --sample-sheet /moto/palab/users/mc4719/RNA_seq_pipeline/data/indices_RNAseq_libraries_demultiplex.csv
```

FastQC of the reads
===================

``` bash
#!/bin/sh                                                                                                 
#SBATCH --account=palab # The account name for the job.                                                   
#SBATCH --job-name=Fastqc  # The job name.                                                        
#SBATCH -c 1                  # The number of cpu cores to use.                                           
#SBATCH -t 2:00:00      # The time the job will take to run                                              
#SBATCH --mem-per-cpu 1gb

module load anaconda/2-5.3.1
source activate prdm9

# First round of sequencing
for i in /moto/palab/users/mc4719/RNA_seq_pipeline/data/Demultiplexed_PRDM9_samples/Combined_lanes/*.gz; do fastqc $i --outdir=/moto/palab/users/mc4719/RNA_seq_pipeline/data/Demultiplexed_PRDM9_samples/Combined_lanes/fastQC_stats/; done       

# Second round of sequencing
for i in /moto/palab/users/mc4719/RNA_seq_pipeline/data/Demultiplexed_PRDM9_samples_2/Combined_lanes/*.gz; do fastqc $i --outdir=/moto/palab/users/mc4719/RNA_seq_pipeline/data/Demultiplexed_PRDM9_samples_2/Combined_lanes/fastQC_stats/; done
```

Combining multiple lanes of two different runs
==============================================

``` bash
# Forward reads
for i in `cat ./ID`; do cat 'Sample_'$i'_L001_R1_001.fastq.gz' 'Sample_'$i'_L002_R1_001.fastq.gz' 'Sample_'$i'_L003_R1_001.fastq.gz' 'Sample_'$i'_L004_R1_001.fastq.gz' > Combined_lanes/'Sample_'$i'_R1_001.fastq.gz'; done

# Reverse reads
for i in `cat ./ID`; do cat 'Sample_'$i'_L001_R2_001.fastq.gz' 'Sample_'$i'_L002_R2_001.fastq.gz' 'Sample_'$i'_L003_R2_001.fastq.gz' 'Sample_'$i'_L004_R2_001.fastq.gz' > Combined_lanes/'Sample_'$i'_R2_001.fastq.gz'; done
```

Changing fastq files to Trinity format
======================================

``` bash
if files are .gz:

Delete everything after the space in the header:
for i in R1_001.fastq.gz*; do zcat $i | sed '/^@/ s/ .*//' > trinity_1_$i; done

Include the "/1" at the end of the headers:
for i in trinity_1_*; do awk '{ if (NR%4==1) { print $1""$2"/1" } else { print } }' $i > new_$i; done

in the reverse reads you change the second step ({ print $1""$2"/2" }, instead of { print $1""$2"/1" }:
for i in trinity_2_*; do awk '{ if (NR%4==1) { print $1""$2"/2" } else { print } }' $i > new_$i; done
```

Denovo assembly using Trinity (04\_trinity\_assembly.sh)
========================================================

``` bash
#!/bin/bash 
#                                                                                                                                                                                          
#SBATCH --account=palab              # The account name for the job.                                                                                                                       
#SBATCH --job-name=trimming_reads    # The job name. 
#SBATCH -N 1
#SBATCH -n 10
#SBATCH -p general                   # may consider running on a bigmem node for large dataset
#SBATCH -e trinity_%A.err            # File to which STDERR will be written
#SBATCH -o trinity_%A.out           # File to which STDOUT will be written
#SBATCH -J trinity_%A               # Job name
#SBATCH --mem=192000                 # Memory requested
#SBATCH --time=3-23:00:00              # Runtime in D-HH:MM:SS
#SBATCH --mail-type=ALL              # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=izabelcavassim@gmail.com # Email to send notifications to

# $1 = comma-separated list of R1 file first batch
# $2 = comma-separated list of R1 file second batch
# $3 = comma-separated list of R2 file first batch
# $4 = comma-separated list of R2 file second batch
# $5 = name of output directory Trinity will create to store results. This must include Trinity in the name, otherwise the job will terminate

module load anaconda/2-5.3.1
source activate prdm9

Trinity --seqType fq --SS_lib_type FR --max_memory 100G --min_kmer_cov 1 --trimmomatic --CPU 32 --left $1,$2 --right $3,$4 --output $5
```

Parallel submission
===================

``` bash
in_directory_reads="/moto/palab/users/mc4719/RNA_seq_pipeline/data/Demultiplexed_PRDM9_samples/Combined_lanes/"
in_directory_reads_2="/moto/palab/users/mc4719/RNA_seq_pipeline/data/Demultiplexed_PRDM9_samples_2/Combined_lanes/"

module load anaconda/2-5.3.1
source activate prdm9
for i in `cat /moto/palab/users/mc4719/RNA_seq_pipeline/data/Demultiplexed_PRDM9_samples/ID`;
do
    out_directory_reads="/moto/palab/users/mc4719/RNA_seq_pipeline/data/Assembled_reads_Trinity_Sample_"$i"/"
   
    # Forward reads first sequencing
    R1_1="$in_directory_reads""trinity_1_Sample_"$i"_R1_001.fastq"
    
    # Reverse reads first sequencing
    R2_1="$in_directory_reads""trinity_2_Sample_"$i"_R2_001.fastq"
    
    # Forward reads second sequencing
    R1_2="$in_directory_reads_2""trinity_1_Sample_"$i"_R1_001.fastq"

    # Reverse reads second sequencing
    R2_2="$in_directory_reads_2""trinity_2_Sample_"$i"_R2_001.fastq"
    dir="$out_directory_reads"
    
    # Printing the sample file
    sbatch 04_trinity_assembly.sh $R1_1 $R1_2 $R2_1 $R2_2 $out_directory_reads
done
```
Blasting PRDM9 and other meiotic genes against each assembly
===================
``` python
from Bio import SeqIO
import pandas as pd
from Bio.Blast.Applications import NcbitblastnCommandline
import re
import subprocess
import numpy as np
from sys import argv
import glob
import os

def parse_fasta(filename):
	file = open(filename, 'r').read() 
	file_separe = file.split('>') 
	file_separe.remove('')
	parse_dict = {}
	header = []
	for entry in file_separe:
		seq = entry.splitlines()
		header = seq[0].split(" ")[0]
		seq = ''.join(seq[1:])
		parse_dict[header] = seq
	return parse_dict


def blasting(db_name = None, evalue = 0.000001, query_name = None, gene_name = None):
	
	#print 'Blasting the sequences against reference query' 

	blastdb_cmd = 'makeblastdb -in {} -dbtype nucl'.format(db_name)

	# Creating a home database
	DB_process = subprocess.Popen(blastdb_cmd,
				      shell=True,
				      stdin=subprocess.PIPE,
				      stdout=subprocess.PIPE,
				      stderr=subprocess.PIPE)
	DB_process.wait()

	blastx_cline = NcbitblastnCommandline(cmd="/moto/palab/users/mc4719/RNA_seq_pipeline/softwares/ncbi-blast-2.9.0+/bin/tblastn", query=str(query_name), db=str(db_name), evalue=evalue, outfmt = 6)
	out, err = blastx_cline()
	list_out = re.split('\n|\t', out)
	del list_out[-1]
	#print(list_out)
	blast_df = pd.DataFrame(np.array(list_out).reshape(len(list_out) // 12,12), columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qends','start','send', 'evalue', 'bitscore'])

	# Changing the type of the data frame:
	blast_df[['pident','length','mismatch','gapopen','qstart','qends','start','send', 'evalue', 'bitscore']] = blast_df[['pident','length','mismatch','gapopen','qstart','qends','start','send', 'evalue', 'bitscore']].apply(pd.to_numeric)

	parsed_assembly = parse_fasta(db_name)
		
	sequences = list()

	# Extracting the sequences that match with PRDM9
	for i in blast_df['sseqid'].tolist(): 
		sequences.append(parsed_assembly[i])
		
	blast_df['sequences'] = sequences
	
	db_name_basename = os.path.basename(db_name)

	# Saving the blast results in csv file:
	name = 'blast_results_{}_{}.csv'.format(gene_name, db_name_basename)
	blast_df.to_csv(name, sep = '\t')

	return(blast_df)

assemblies = [(f) for f in glob.glob("/moto/palab/users/mc4719/RNA_seq_pipeline/data/Combined_assemblies_TRINITY/*/*.fasta")] 
print(assemblies)


genes = ["MEI4", "MRE11A", "RAD50", "HORMAD1", "WDR61", "MEI1", "HORMAD1"]

for j in genes:
	for i in assemblies:
		print(i)
		query = "/moto/palab/users/mc4719/RNA_seq_pipeline/data/Query_REC_genes/{}_main_isoform.fasta".format(j)
		t = blasting(db_name = i, query_name=query, gene_name = j)
		print(t)
``` 
Estimating transcript abundance (RSEM)
===================

``` bash
#!/bin/sh                                                                                                                     #SBATCH --account=palab # The account name for the job.                                                                      #SBATCH --job-name= transcripts_quantification  # The job name.                                                               #SBATCH -N 1                  # The number of cpu cores to use.                                                             #SBATCH --time=2:00:00       # The time the job will take to run                                                             #SBATCH --mem=120G                                                                                                          #SBATCH -e RSEM_%A.err            # File to which STDERR will be written                                                    #SBATCH -o RSEM_%A.out           # File to which STDOUT will be written    
#SBATCH --mail-type=ALL              # Type of email notification- BEGIN,END,FAIL,ALL                                      #SBATCH --mail-user=izabelcavassim@gmail.com # Email to send notifications to                                                                                              
                                               
assembly_prefix=$(basename $5 |sed 's/.fasta//g')

module load anaconda/2-5.3.1
source activate prdm9

# Adding RSEM to the path
PATH=$PATH:/moto/palab/users/mc4719/RNA_seq_pipeline/softwares/RSEM
             
fasta_directory="/moto/palab/users/mc4719/RNA_seq_pipeline/data/Combined_assemblies_TRINITY/"
output_directory="/moto/palab/users/mc4719/RNA_seq_pipeline/data/Transcripts_abundance/"

#NOTE: if your assembled your reads with Trinity by providing lists of left and right reads                                                                                                                               
/moto/palab/users/mc4719/RNA_seq_pipeline/softwares/trinityrnaseq-Trinity-v2.8.5/util/align_and_estimate_abundance.pl --seqType fq \
	--left $1  \
	--right $2  \
	--transcripts $3  \
	--output_dir $4 \
        --thread_count 1 \
	--est_method RSEM \
        --aln_method bowtie \
        --prep_reference \
	--trinity_mode \
```
Parallel submission of transcript abundance (RSEM)
===================

```bash
in_directory_reads="/moto/palab/users/mc4719/RNA_seq_pipeline/data/Demultiplexed_PRDM9_samples/Combined_lanes/"
in_directory_reads_2="/moto/palab/users/mc4719/RNA_seq_pipeline/data/Demultiplexed_PRDM9_samples_2/Combined_lanes/"
in_directory_reads_3="/moto/palab/users/mc4719/RNA_seq_pipeline/data/Demultiplexed_PRDM9_samples/Combined_lanes_runs/"

module load anaconda/2-5.3.1
source activate prdm9
for i in `cat /moto/palab/users/mc4719/RNA_seq_pipeline/data/Demultiplexed_PRDM9_samples/ID`;
do
    out_directory_reads="/moto/palab/users/mc4719/RNA_seq_pipeline/data/Assembled_reads_Trinity_Sample_"$i"/"
    assembly_output="/moto/palab/users/mc4719/RNA_seq_pipeline/data/Combined_assemblies_TRINITY/Assembled_reads_Trinity_Sample_"$i"/"
    assembly_trinity="/moto/palab/users/mc4719/RNA_seq_pipeline/data/Combined_assemblies_TRINITY/Assembled_reads_Trinity_Sample_"$i"/Assembled_reads_Trinity_Sample_"$i".fasta"

    #/moto/palab/users/mc4719/RNA_seq_pipeline/data/Assembled_reads_Trinity_Sample_A01_S1
    # Forward reads first sequencing
    R1_1="$in_directory_reads""trinity_1_Sample_"$i"_R1_001.fastq"
    
    # Reverse reads first sequencing
    R2_1="$in_directory_reads""trinity_2_Sample_"$i"_R2_001.fastq"
    
    # Forward reads second sequencing
    R1_2="$in_directory_reads_2""trinity_1_Sample_"$i"_R1_001.fastq"

    # Reverse reads second sequencing
    R2_2="$in_directory_reads_2""trinity_2_Sample_"$i"_R2_001.fastq"
    #dir="$out_directory_reads"

    # Concatenating files 
    #cat $R1_1 $R1_2 > $in_directory_reads_3"trinity_1_Sample_"$i"_R1_001_combined.fastq"
    #cat $R2_1 $R2_2  > $in_directory_reads_3"trinity_1_Sample_"$i"_R2_001_combined.fastq"
    
    R1=$in_directory_reads_3"trinity_1_Sample_"$i"_R1_001_combined.fastq"
    R2=$in_directory_reads_3"trinity_1_Sample_"$i"_R2_001_combined.fastq"

    # Printing the sample file
    sbatch 09_transcript_quantification.sh $R1 $R2 $assembly_trinity $assembly_output
done
