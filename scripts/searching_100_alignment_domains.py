# First parse the hitdata file 
# Extract the coordinates of the superfamily for each gene
# Compare it againts the 

# Looking at the domains (zinc fingers) of the zcwpw1
import glob
import pandas as pd
from sys import argv
import numpy as np
from collections import OrderedDict 
import os
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
from ete3 import Tree
from io import StringIO
from Bio import Phylo
# Constructing a tree
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import DistanceCalculator
import glob

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

def writing_fasta_dict(parse_dict,filename):
	gene = os.path.basename(filename)
	print(gene)
	# Saving as a fasta file:
	with open(filename+".fasta", 'w') as f:
		for names, seq in parse_dict.items():
				if len(seq) >=1:
					#print(seq)
					#f.write('>{}_{}\n'.format(names, gene))
					f.write('>{}\n'.format(names))
					f.write('{}\n\n'.format(seq[0]))

# Reading domains:
domains_place = pd.read_csv("/Users/PM/Dropbox/PHD/Recombination_columbia/Recombination_analysis/data/Query_REC_genes_prot_all_seqs/hitdata_PRDM9_homo_sapiens_domains.txt", sep="\t", comment="#")
print(domains_place)

results = []

for index, row in domains_place.iterrows():

	gene_ = row['Query']
	gene_ = gene_.split(">")[1]
	gene_ = gene_.split("_")[0]
	superfamily = row['Short name']
	accession = row['Accession']
	superfamily = superfamily + "_" + accession
	species = row['Query']
	from_ = row['From']
	to_ = row['To']
	#print([from_, to_])

	# parsing fasta file
	average_completness = list()
	gene_parsed = parse_fasta("/Users/PM/Dropbox/PHD/Recombination_columbia/Recombination_analysis/data/Query_REC_genes_prot_all_seqs/{}.fasta".format(gene_))
	#print(gene_)
	for species in gene_parsed:
		#print(species)
		sequence = gene_parsed[species][from_:to_]

		len_seq = len(sequence)
		len_ns = sequence.count('-')
		#print(sequence)
		average_completness.append(len_ns/len_seq)
	average_completness = np.max(average_completness)
	results.append([gene_ , superfamily, average_completness, len(gene_parsed.keys())])

dataframe_results = pd.DataFrame(results) 
dataframe_results.columns = ['Gene', 'Domain', 'Average_completeness', 'Num_seqs']
dataframe_results.to_csv("Candidates_incompleteness_domain_test.csv", index=False)
print(dataframe_results)
test = dataframe_results[['Gene','Num_seqs', 'Average_completeness']].groupby('Gene').mean()
test = test.sort_values('Num_seqs')
test.to_csv("Candidates_max_incompleteness_domain_test.csv", index=True)
print(test[(test['Num_seqs']>=100) & (test['Average_completeness']<=0.05)].to_string())
print(test.to_string())
	#if ("PWWP" in domain) or ("MSH6_like" in domain):
	#if ("zf-CW" in domain):
		#print(">{species}".format(species=species))
		#print(zcwpw2[species][0][from_:to_])
