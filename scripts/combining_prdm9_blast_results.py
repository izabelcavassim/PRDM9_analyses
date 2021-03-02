# Combining the presence absence dataframes
# Making a 0 and 1 matrix
# Making a fasta file of each sequence
# python combining_prdm9_candidates.py /Users/PM/Desktop/Blast_results/ 
import glob
import pandas as pd
from sys import argv
import numpy as np

def pres_abs(lists, default):
	results = []
	for i in default:
		if i in lists:
			results.append(1)
		else:
			results.append(0)
	return(results)


dir_dfs = argv[1]

dfs = [(f) for f in glob.glob(dir_dfs+"blast_results_recombination_genes*.csv")]
print(dfs)

# Homo sapiens 
all_genes = pd.read_csv(dir_dfs+"blast_results_recombination_genes"+"_Homo_sapiens.fa.csv", sep="\t")
candidates = sorted(all_genes["Gene"].tolist())
candidates_names = [s.replace('.fasta', '') for s in candidates]


dictionary_presence_absence = pd.DataFrame(columns=candidates_names)
print(dictionary_presence_absence)

dir_genes = {}
df_genes = []
species_names = []
for f in dfs:
	species = f.split("blast_results_recombination_genes_")[1]
	species = species.split(".fa.csv")[0]
	species_names.append(species)
	t = pd.read_csv(f, sep="\t")
	gene_list = t["Gene"].tolist()
	dir_genes[species] = sorted(gene_list)

	# creating dataframe:
	df = pd.DataFrame(columns=candidates_names)
	print(len(gene_list))
	print(len(candidates))
	t = pres_abs(gene_list, candidates)
	print(t)

	df.loc[0] = t
	df_genes.append(df)

all_com = pd.concat(df_genes)
all_com['Species'] = species_names
all_com = all_com.set_index('Species')