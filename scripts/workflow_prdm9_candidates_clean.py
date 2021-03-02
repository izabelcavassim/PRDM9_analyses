######################################################################################################
# This workflow was created in order to run the analysis of branch lentghs across the 100-way alignments
# alignments were assesed at UCSC and further details on the analysis are found in the README file
#
#######################################################################################################

from gwf import Workflow
import glob
import os
import pandas as pd
gwf = Workflow()

def creating_fasta_files(ucsc_file_dir,transcript_out, types):
	inputs=[f'{ucsc_file_dir}']
	outputs=[f'/home/mica16/MouseLemur/parsed_UCSC_100_way_canonical_genes/knownCanonical.exon'+ f'{types}_transcripts/status.txt']
	options = {
	'memory':'8g',
	'cores':'8',
	'walltime':'12:00:00',
	'account': 'NChain'
	}
	spec = '''
	python split.py {ucsc_file_dir} {transcript_out}
	'''.format(ucsc_file_dir=ucsc_file_dir, transcript_out=transcript_out)
	print(spec)
	print(outputs)
	return inputs, outputs, options, spec
	
def parsing_UCSC_fasta_files(in_directory, out_directory_1, out_directory_2):
	inputs = []
	outputs = [f'{out_directory_2}'+'uc031tks.1.fasta']
	options = {
	'memory':'8g',
	'cores':'8',
	'walltime':'04:00:00',
	'account': 'NChain'
	}
	spec = f'python parsing_UCSC_data.py {in_directory} {out_directory_1} {out_directory_2}'
	print(spec)
	print(outputs)
	return inputs, outputs, options, spec


def subsetting_fasta_files(canonical_genes_directory, subset_canonical_genes, subset_tree_dir):
	inputs = []
	outputs = [f'{subset_tree_dir}']
	options = {"memory": "2g", "walltime":"05:00:00"}
	spec='''
	python /home/mica16/MouseLemur/scripts/codeml/subsetting_fasta_rec_genes.py {canonical_genes_directory} {subset_canonical_genes} {subset_tree_dir}
	'''.format(canonical_genes_directory=canonical_genes_directory, subset_canonical_genes=subset_canonical_genes, subset_tree_dir=subset_tree_dir)
	print(spec)
	print(outputs)
	return inputs, outputs, options, spec

def presence_absence_genes(canonical_dir, species_data, genes_rec_metadata, genes_names_conversion, query_directory, genomes_directory, query_directory_prot):
	inputs = [genomes_directory]
	outputs = ["/home/mica16/MouseLemur/scripts/Blast_results/"+"blast_results_recombination_genes_test_{}.csv".format(os.path.basename(genomes_directory))] #blast_results_recombination_genes_Salmo_salar.fa.csv 
	options = {"memory": "6g","walltime":"24:00:00"}
	spec='''                                                                                                                                
	source activate prdm9
	python blasting_prdm9_candidates.py {canonical_dir} {species_data} {genes_rec_metadata} {genes_names_conversion} {query_directory} {genomes_directory} {query_directory_prot}                                                                                                           
	'''.format(canonical_dir=canonical_dir, species_data=species_data, genes_rec_metadata=genes_rec_metadata, genes_names_conversion=genes_names_conversion, query_directory=query_directory, genomes_directory=genomes_directory, query_directory_prot=query_directory_prot)
	print(spec)
	return inputs, outputs, options, spec

def combinig_top_hits(dir_dfs, dir_out, output_blast_list, genes_hits, gene):
	inputs = output_blast_list
	#outputs = [(f'{dir_out}{f}') for f in genes_hits]
	outputs = [dir_out+"cdd_results_{}.csv".format(gene)]
	options = {"memory": "6g","walltime":"4:00:00"}
	spec='''                                                                                                                                
	source activate prdm9
	python Combining_prdm9_5_top_hits_blast_cdd.py {dir_dfs} {dir_out} {gene}                                                                                                      
	'''.format(dir_dfs=dir_dfs, dir_out=dir_out, gene=gene)
	print(spec)
	return inputs, outputs, options, spec

def preparing_Bayesfiles(blast_results, dir_out, genes_hits):
	inputs = [blast_results]
	outputs = [f'{dir_out}ZNF84.txt']
	print(outputs)
	options = {"memory": "6g","walltime":"4:00:00"}
	spec='''
	source activate prdm9                                                                                                                               
	python Separating_presence_absence_files_Bayestraits.py {blast_results} {dir_out}                                                                                                   
	'''.format(blast_results=blast_results, dir_out=dir_out)
	print(spec)
	return inputs, outputs, options, spec  

def preparing_Bayesfiles_2(pairs_presence_absence_file, species_tree, gene, output_dir):
	inputs = [f'{dir_out}{gene}.txt']
	outputs = [f'{output_dir}{gene}/pruned_tree.newick', f'{output_dir}{gene}/reference_tree.newick', f'{output_dir}{gene}/trait_table.tab']
	print(outputs)
	options = {"memory": "6g","walltime":"4:00:00"}
	# format_tree_and_trait_table.py -i zcwpw2_prdm9.txt -t Species_used_13_01_2020_replaced_names.nwk -n -o zcwpw2  
	spec='''
	source activate picrust1                                                                                                                               
	format_tree_and_trait_table.py -i {pairs_presence_absence_file} -t {species_tree} -n -o {output_dir}{gene}                                                                                         
	'''.format(pairs_presence_absence_file=pairs_presence_absence_file, species_tree=species_tree, output_dir=output_dir, gene=gene)
	print(spec)
	return inputs, outputs, options, spec  

def preparing_Bayesfiles_pruned(pairs_presence_absence_file, gene, output_dir):
	inputs = [f'{dir_out}{gene}.txt']
	outputs = [f'{output_dir}{gene}/pruned_tree.newick', f'{output_dir}{gene}/reference_tree.newick', f'{output_dir}{gene}/trait_table.tab']
	print(outputs)
	tree = pairs_presence_absence_file[:-4]
	tree = tree+".nw"
	options = {"memory": "6g","walltime":"4:00:00"}
	# format_tree_and_trait_table.py -i zcwpw2_prdm9.txt -t Species_used_13_01_2020_replaced_names.nwk -n -o zcwpw2  
	spec='''
	source activate picrust1                                                                                                                               
	format_tree_and_trait_table.py -i {pairs_presence_absence_file} -t {tree} -n -o {output_dir}{gene}                                                                                         
	'''.format(pairs_presence_absence_file=pairs_presence_absence_file, tree=tree, output_dir=output_dir, gene=gene)
	print(spec)
	return inputs, outputs, options, spec  

def Bayestrait_run_models(gene, output_dir,all_runs_dir, name):
	gene_dir=f'{output_dir}{gene}'
	inputs = [f'{gene_dir}/pruned_tree.newick', f'{gene_dir}/reference_tree.newick', f'{gene_dir}/trait_table.tab']
	print(inputs)
	outputs = [f'{all_runs_dir}{gene}_dependent_model_with_restrictions.txt_results_{name}', f'{all_runs_dir}{gene}_independent_model_with_restrictions.txt_results_{name}']
	options = {"memory": "6g","walltime":"2:00:00"}
	spec='''
	for FILE in /home/mica16/MouseLemur/scripts/Bayestraits_candidates/Models_Bayes_traits/*.txt; 
	do basename "$FILE"; 
	f="$(basename -- $FILE)"; 
	/home/mica16/MouseLemur/software/BayesTraitsV3.0.2-Linux/BayesTraitsV3 {gene_dir}/pruned_tree.newick  {gene_dir}/trait_table.tab < $FILE > {all_runs_dir}{gene}_$f"_results_{name}"; 
	done
	'''.format(gene_dir=gene_dir, gene=gene, all_runs_dir=all_runs_dir,name=name)
	print(spec)
	return inputs, outputs, options, spec 

def Bayestrait_parse_old(gene, bayes_results_dir, name):
	inputs = [f'{bayes_results_dir}{gene}_dependent_model_with_restrictions.txt_results', f'{bayes_results_dir}{gene}_independent_model_with_restrictions.txt_results']
	print(inputs)
	outputs = [f'results_Bayestraits_{name}.txt']
	options = {"memory": "6g","walltime":"2:00:00"}
	spec='''
	for i in {bayes_results_dir}*results; do IFS='_' read -r -a array <<< "$i"; python parsing_calc_lik.py "$array[0]"; done > results_Bayestraits_{name}.txt
	'''.format(bayes_results_dir=bayes_results_dir, name=name)
	print(spec)
	return inputs, outputs, options, spec              

def presence_absence_genes_wgs(canonical_dir, species_data, genes_rec_metadata, genes_names_conversion, query_directory, genomes_directory):
	inputs = [genomes_directory]
	outputs = ["/home/mica16/MouseLemur/scripts/Blast_results_wgs/"+"blast_results_recombination_genes_{}_wgs.csv".format(os.path.basename(genomes_directory))] #blast_results_recombination_genes_Salmo_salar.fa.csv
	print(outputs)
	options = {"memory": "100g","walltime":"02:00:00"}
	spec='''                                                                                                                                
	source activate prdm9
	python blasting_candidates_wgs.py {canonical_dir} {species_data} {genes_rec_metadata} {genes_names_conversion} {query_directory} {genomes_directory}                                                                                                           
	'''.format(canonical_dir=canonical_dir, species_data=species_data, genes_rec_metadata=genes_rec_metadata, genes_names_conversion=genes_names_conversion, query_directory=query_directory, genomes_directory=genomes_directory)
	print(spec)
	return inputs, outputs, options, spec

def Bayestrait_parse_tries(gene, bayes_results_dir):
	inputs = [f'{bayes_results_dir}{gene}_dependent_model_with_restrictions.txt_results_49']
	outputs = [f'/home/mica16/MouseLemur/scripts/Bayestraits_tries_results/results_{gene}_H0_Ha3.csv', f'/home/mica16/MouseLemur/scripts/Bayestraits_tries_results/results_{gene}_H0_Ha2.csv', f'/home/mica16/MouseLemur/scripts/Bayestraits_tries_results/results_{gene}_H0_Ha1.csv',  f'/home/mica16/MouseLemur/scripts/Bayestraits_tries_results/{gene}_hist_plot.png']
	options = {"memory": "1g","walltime":"00:20:00"}
	spec='''
	source activate py36
	python parsing_Bayestraits_50_tries.py {gene} {bayes_results_dir}
	'''.format(bayes_results_dir=bayes_results_dir, gene=gene)
	print(spec)
	print(inputs)
	print(outputs)
	return inputs, outputs, options, spec     

######### Parsing UCSC data ##########
# Aminoacid files
ucsc_file_dir="/home/mica16/MouseLemur/parsed_UCSC_100_way_canonical_genes/alignments/alignments/knownCanonical.exonAA/knownCanonical.exonAA.fa"
transcript_out="/home/mica16/MouseLemur/parsed_UCSC_100_way_canonical_genes/alignments/alignments/knownCanonical.exonAA/" 
gwf.target_from_template("splitting_AA_files", creating_fasta_files(ucsc_file_dir, transcript_out, types = 'AA'))

# Nucleotide files
ucsc_file_dir="/home/mica16/MouseLemur/parsed_UCSC_100_way_canonical_genes/alignments/alignments/knownCanonical.exonNuc/knownCanonical.exonNuc.fa"
transcript_out="/home/mica16/MouseLemur/parsed_UCSC_100_way_canonical_genes/alignments/alignments/knownCanonical.exonNuc/"
gwf.target_from_template("splitting_Nuc_files", creating_fasta_files(ucsc_file_dir, transcript_out, types = 'Nuc'))

# First the AA files: putting every transcript into a fasta file and combining it into genes
indict1 = "/home/mica16/MouseLemur/parsed_UCSC_100_way_canonical_genes/alignments/alignments/knownCanonical.exonAA/"
out1 = "/home/mica16/MouseLemur/parsed_UCSC_100_way_canonical_genes/knownCanonical.exonAA_transcripts/"
out2 = "/home/mica16/MouseLemur/parsed_UCSC_100_way_canonical_genes/knownCanonical.exonAA_genes/"

gwf.target_from_template("parsing_canonical_AA", parsing_UCSC_fasta_files(in_directory=indict1, out_directory_1=out1, out_directory_2=out2))

# Second the nucleotide files: putting every transcript into a fasta file and combining it into genes
indict1 = "/home/mica16/MouseLemur/parsed_UCSC_100_way_canonical_genes/alignments/alignments/knownCanonical.exonNuc/"
out1 = "/home/mica16/MouseLemur/parsed_UCSC_100_way_canonical_genes/knownCanonical.exonNuc_transcripts/"
out2 = "/home/mica16/MouseLemur/parsed_UCSC_100_way_canonical_genes/knownCanonical.exonNuc_genes/"
gwf.target_from_template("parsing_canonical_Nuc", parsing_UCSC_fasta_files(in_directory=indict1, out_directory_1=out1, out_directory_2=out2))


# Looking at all the Refseq sequences available
## Blasting genes againts the extra genomes
canonical_dir="/home/mica16/MouseLemur/parsed_UCSC_100_way_canonical_genes/knownCanonical.exonAA_genes/" #### HERE THE QUERY IS PROTEIN, WE NEED TO ADAPT THAT TO THE DATATYPE
species_metadata="/home/mica16/MouseLemur/parsed_UCSC_100_way_canonical_genes/Metadata/Species_transitions_340.csv"
genes_rec_metadata = "/home/mica16/MouseLemur/parsed_UCSC_100_way_canonical_genes/Metadata/Recombination_decode_results_update_04_06_20.txt"
genes_names_conversion = "/home/mica16/MouseLemur/parsed_UCSC_100_way_canonical_genes/Metadata/gene_names_updated.txt"
query_directory = "/home/mica16/MouseLemur/Query_sequences_Homo_sapiens/"
query_directory_prot = "/home/mica16/MouseLemur/Query_REC_genes_prot_all_seqs/"
genomes_directory = "/home/mica16/MouseLemur/Refseq_genomes/"

# Looking at a subset of species
assemblies = [(f) for f in glob.glob("{genomes_directory}*.faa".format(genomes_directory=genomes_directory))]
print(assemblies)
for i in assemblies:
	assembly_basename = os.path.basename(i)
	gwf.target_from_template("Presence_absence_rec_genes"+"_"+assembly_basename, presence_absence_genes(canonical_dir, species_metadata, genes_rec_metadata, genes_names_conversion, query_directory, i, query_directory_prot))


# Extracting all the tophits to individual files
dir_dfs = "/home/mica16/MouseLemur/scripts/Blast_results/"
dir_out = "/home/mica16/MouseLemur/scripts/Fasta_parsed_candidates/"
blast_results = [(f) for f in glob.glob("/home/mica16/MouseLemur/scripts/Blast_results/*._protein.faa.csv")]
homo_sapiens_genes = pd.read_csv("{dir_dfs}blast_results_recombination_genes_test_Homo_sapiens_protein.faa.csv".format(dir_dfs=dir_dfs), sep="\t")
genes = homo_sapiens_genes["Gene"].tolist()
genes = set(genes)
#for gene in genes:
#        gwf.target_from_template("parsing_fasta_top5_hits_{gene}".format(gene=gene), combinig_top_hits(dir_dfs=dir_dfs, dir_out=dir_out, output_blast_list=blast_results, genes_hits=genes, gene=gene))


# Preparing files to Bayestraits
blast_file = "/home/mica16/MouseLemur/scripts/blast_results_presence_absence_genes_all_340.csv"
blast_file = "blast_results_presence_absence_genes_all_340_superfamily_domains_identified_test_specific_0.0001_30_06_20.csv"

dir_out = "/home/mica16/MouseLemur/scripts/Bayestraits_candidates/"
gwf.target_from_template("preparing_Bayestrait", preparing_Bayesfiles(blast_results=blast_file, dir_out=dir_out, genes_hits=genes))

Bayestraits_results = "/home/mica16/MouseLemur/scripts/Bayestraits_results/"
species_tree = "/home/mica16/MouseLemur/scripts/species_tree_339_changed_replacements_to_old_names.nwk"
genes_files = [(f) for f in glob.glob("/home/mica16/MouseLemur/scripts/Bayestraits_candidates/*.txt")]
for g in genes_files:
	gene_name = os.path.basename(g)
	gene_name = gene_name[:-4]
	if gene_name == "MSH5-SAPCD1":
		gene_job = "MSH5_SAPCD1"
		gwf.target_from_template(f'preparing_Bayestrait_{gene_job}', preparing_Bayesfiles_2(pairs_presence_absence_file=g, species_tree=species_tree, gene=gene_name, output_dir=dir_out))
		gwf.target_from_template(f'Running_Bayestrait_{gene_job}', Bayestrait_run_models(gene=gene_name, output_dir=dir_out, all_runs_dir=Bayestraits_results))
	else:
		gwf.target_from_template(f'preparing_Bayestrait_{gene_name}', preparing_Bayesfiles_2(pairs_presence_absence_file=g, species_tree=species_tree, gene=gene_name, output_dir=dir_out))
		gwf.target_from_template(f'Running_Bayestrait_{gene_name}', Bayestrait_run_models(gene=gene_name, output_dir=dir_out, all_runs_dir=Bayestraits_results))

### Running the prunned version
Bayestraits_results_pruned = "/home/mica16/MouseLemur/scripts/Bayestraits_results_pruned/"
dir_out = "/home/mica16/MouseLemur/scripts/Bayestraits_candidates_pruned/"
genes_files = [(f) for f in glob.glob("/home/mica16/MouseLemur/scripts/Bayestraits_candidates_pruned/*.txt")]


# # commenting for now
# for i in range(1,50):
# 	for g in genes_files:
# 		gene_name = os.path.basename(g)
# 		gene_name = gene_name[:-4]
# 		if gene_name == "MSH5-SAPCD1":
# 			gene_job = "MSH5_SAPCD1"
# 			if i == 1:
# 				gwf.target_from_template(f'preparing_Bayestrait_{gene_job}_pruned', preparing_Bayesfiles_pruned(pairs_presence_absence_file=g, gene=gene_name, output_dir=dir_out))
# 			gwf.target_from_template(f'Running_Bayestrait_{gene_job}_pruned_{i}', Bayestrait_run_models(gene=gene_name, output_dir=dir_out, all_runs_dir=Bayestraits_results_pruned, name=i))
# 		else:
# 			if i == 1:
# 				#gwf.target_from_template(f'preparing_Bayestrait_{gene_name}_pruned', preparing_Bayesfiles_pruned(pairs_presence_absence_file=g, gene=gene_name, output_dir=dir_out))
# 			#gwf.target_from_template(f'Running_Bayestrait_{gene_name}_pruned_{i}', Bayestrait_run_models(gene=gene_name, output_dir=dir_out, all_runs_dir=Bayestraits_results_pruned, name=i))

# #gwf.target_from_template(f'Parsing_Bayestrait_all_pruned', Bayestrait_parse(gene='ZNF84', bayes_results_dir=Bayestraits_results, name="50_tries_epsilon_zero_pruned"))

# Parsing all genes all tries (100)
for g in genes_files:                                                                                                       
	gene_name = os.path.basename(g)                                                                                     
	gene_name = gene_name[:-4]
	if gene_name != "MSH5_SAPCD1":
		gwf.target_from_template(f'Parsing_Bayestrait_{gene_name}', Bayestrait_parse_tries(gene=gene_name, bayes_results_dir=Bayestraits_results_pruned))


#genomes_directory = "/home/mica16/MouseLemur/Extra_genomes/fna_files/"
## Blasting genes againts the whole genonme sequences
#assemblies = [(f) for f in glob.glob("{genomes_directory}*.fna".format(genomes_directory=genomes_directory))]
#for i in assemblies:
#        assembly_basename = os.path.basename(i)
#
#         gwf.target_from_template("Presence_absence_wgs"+"_"+assembly_basename, presence_absence_genes_wgs(canonical_dir, species_metadata, genes_rec_metadata, genes_names_conversion, query_directory, i))

# Running proteinortho
#gwf.target_from_template('orthology_detection', proteinortho_run(directory = '/home/mica16/MouseLemur/Extra_genomes', project_name = 'orthologous_37_species'))
#print(inputs)

