from Bio import SeqIO
import pandas as pd
from Bio.Blast.Applications import NcbitblastnCommandline
import re
import subprocess
import numpy as np
from sys import argv
import glob
import os

canonical_exon_dir = argv[1] # /Users/PM/Desktop/knownCanonical.exonNuc_genes/
species_metadata = argv[2] # /Users/PM/Dropbox/PHD/Recombination_columbia/Recombination_analysis/Paper_folder/Species_transitions.csv
genes_rec_metadata = argv[3] # '/Users/PM/Dropbox/PHD/Recombination_columbia/Recombination_analysis/data/Recombination_decode_results.txt'
gene_names_conversion = argv[4] # "/Users/PM/Dropbox/PHD/Recombination_columbia/Recombination_analysis/data/gene_names_updated.txt"
query_directory = argv[5] # "/Users/PM/Dropbox/PHD/Recombination_columbia/Recombination_analysis/data/Query_REC_genes/
genomes_file = argv[6] # "/Users/PM/Dropbox/PHD/Recombination_columbia/Recombination_analysis/Paper_folder/"
query_directory_prot = argv[7] # Query_REC_genes_prot_all_seqs

# This script blast recombination gene candidates againts a WGS database. IF the genome of interest is .fa, 
# then we do a blastx of the genome againts the query, 
# If the genome is .fna then we do a blastn. The query is always the human sequence.
# The inability to detect a sequence with a e-value of 10e-10, implies that the gene is not present in the genome
# Later we concatennate the sequences into fasta file and we make matrices of present and absence

def parse_fasta(filename):
	file = open(filename, 'r').read() #opening and reading the fasta file, putting it in a object called file
	file_separe = file.split('>') #spliting each entry by the > 
	file_separe.remove('')
	parse_dict = {}
	header = []
	for entry in file_separe:
		seq = entry.splitlines()
		#print(seq)
		header = seq[0].split(" ")[0] #these are the first elements of the list
		seq = ''.join(seq[1:]) #joining the sequences 
		parse_dict[header] = seq
	return parse_dict

def writing_fasta_dict(parse_dict,filename, direc):
	# Saving as a fasta file:
	with open(direc+filename+".fasta", 'w') as f:
		for names, seq in parse_dict.items():
			if names == "Homo_sapiens":
				f.write('>{}\n'.format(names))
				f.write('{}\n\n'.format(seq))

def writing_fasta_dict_prot(parse_dict,filename, direc_prot):
	# Saving as a fasta file:
	with open(direc_prot+filename+".fasta", 'w') as f:
		for names, seq in parse_dict.items():
			f.write('>{}_{}\n'.format(filename, names))
			f.write('{}\n\n'.format(seq))


def blasting(db_name = None, evalue = 1e-5, query_name = None, gene_name = None, make_database=True):
	
	#print 'Blasting the sequences against reference query' 

	#print 'Making a blast data set'
	#if make_database ==True:
	# 	blastdb_cmd = 'makeblastdb -in {} -dbtype prot'.format(db_name)

	# 	# Creating a home database
	# 	DB_process = subprocess.Popen(blastdb_cmd,
	# 							  shell=True,
	# 							  stdin=subprocess.PIPE,
	# 							  stdout=subprocescs.PIPE,
	# 							  stderr=subprocess.PIPE)
	# 	DB_process()

	blastx_cline = NcbitblastnCommandline(cmd="blastp", query=str(query_name), db=str(db_name), evalue=evalue, outfmt = 6)
	print(blastx_cline)
	out, err = blastx_cline()
	list_out = re.split('\n|\t', out)
	del list_out[-1]
	#print(list_out)
	blast_df = pd.DataFrame(np.array(list_out).reshape(len(list_out) // 12,12), columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qends','start','send', 'evalue', 'bitscore'])

	print("here are the blast results")
	print(blast_df)
	# Changing the type of the data frame:
	blast_df[['pident','length','mismatch','gapopen','qstart','qends','start','send', 'evalue', 'bitscore']] = blast_df[['pident','length','mismatch','gapopen','qstart','qends','start','send', 'evalue', 'bitscore']].apply(pd.to_numeric)
	parsed_assembly = parse_fasta(db_name)
		
	sequences = list()

	# Extracting the sequences that match with PRDM9
	for i in blast_df['sseqid'].tolist(): 
		sequences.append(parsed_assembly[i])
		
	blast_df['sequences'] = sequences
	
	blast_df['Gene'] = os.path.basename(gene_name)
	return(blast_df[:5])


def searching_for_genes_transcriptomic_data(genes, species_dir):

	## Looking at the extra genomes:
	#assemblies = [(f) for f in glob.glob("{species_dir}*.fa".format(species_dir=species_dir))]
	#print(assemblies)

	# Creating the index for the genomes
	#!/usr/bin/env python
	#for file in assemblies:
	#subprocess.call("makeblastdb -in " + species_dir +" -dbtype prot", shell=True)

	t = list()
	for j in genes:
		print(j)
		query = "{j}".format(j=j)
		try:
			t.append(blasting(db_name = species_dir, query_name=query, gene_name = j))
			print(t)
		except:
			pass
						#print(j)
	# Writing the dataframe
	blast_df = pd.concat(t)
	assembly = os.path.basename(species_dir)
	name = "/home/mica16/MouseLemur/scripts/Blast_results/"+'blast_results_recombination_genes_test_{}.csv'.format(assembly)
	blast_df.to_csv(name, sep = '\t')

def reading_candidates(genes_rec_metadata):
	candidates = pd.read_csv(genes_rec_metadata, sep ="\t")
	return(candidates)

# Genes conversions from transcript name to gene name:
def genes_converter(gene_names_conversion):
	gene_transcript = pd.read_csv(gene_names_conversion, sep="\t", header = 0, names= ["transcript_id", "gene_name"])
	
	dict_gene_transcript = dict()

	for index, row in gene_transcript.iterrows():
		name1 = row['gene_name']
		name2 = row['transcript_id']
		if name1 not in dict_gene_transcript:
			dict_gene_transcript[name1] = [name2]
		else:
			dict_gene_transcript[name1].append(name2)
	return(dict_gene_transcript)
	

def extracting_rec_query_genes(genes, candidates, converter, dir_exon_dir, dir_query_rec_genes, dir_query_rec_genes_prot):
	#print(genes)
	for i in range(len(genes)):
		if genes[i] not in ["intergenic", "inversion", "MEI4", "NBS1", "BHMG1"]:
			try:
				#print([genes[i], converter[genes[i]]])
				transcripts_IDs = converter[genes[i]]
			except:
				alias = candidates.iloc[i,8]
				transcripts_IDs = converter[alias]
			for g in transcripts_IDs:
				g = dir_exon_dir+g+".fasta" #dir_exon_dir
				try:
					gene = parse_fasta(g)
					#writing_fasta_dict(gene,genes[i], dir_query_rec_genes)
					#writing_fasta_dict_prot(gene, genes[i], dir_query_rec_genes_prot)
				except:
					pass

candidates = reading_candidates(genes_rec_metadata)
converter = genes_converter(gene_names_conversion)
genes = candidates["Gene"].tolist()

extracting_rec_query_genes(genes=genes, candidates= candidates, converter=converter, dir_exon_dir=canonical_exon_dir, dir_query_rec_genes=query_directory, dir_query_rec_genes_prot=query_directory_prot)

genes_directory = [(f) for f in glob.glob(canonical_exon_dir+"*.fasta")]
species_df = pd.read_csv(species_metadata, sep=";")
species = species_df["Species"].tolist()
genes_found = [(f) for f in glob.glob(query_directory + "*.fasta")]

searching_for_genes_transcriptomic_data(genes=genes_found, species_dir=genomes_file)

