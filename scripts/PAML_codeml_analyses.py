# Python pipeline to run codeml using the Biopython library ()
import sys
from Bio.Phylo.PAML import codeml

options = sys.argv[1]
print(options)
# Instead of using a control file we can change the arguments within python 
# Files used in these analyzes are found in: the data folder

working_directory = "/Users/PM/Desktop/COLUMBIA_files/codeml_analysis/zcwpw2_canids"

# Analyses with only platypus (files are found in: )
if options == "null_platypus":

	gene_file = f'{working_directory}/Alignments/ZCWPW2_alignment_platypus.fas'
	tree_file = f'{working_directory}/Trees/ZCWPW2_platypus.nwk'
	out_file = f'{working_directory}/Results/results_null_zcwpw2_platypus.txt'
	options = f'{working_directory}/Models/ctl_branchmodel0_null_zcwpw2.txt'

	# Add the null tree in here
	cml = codeml.Codeml(alignment = gene_file, tree = tree_file,
                    out_file = out_file, working_dir=working_dir)

	cml.read_ctl_file(options)
	cml.print_options()
	cml.run(verbose = True)

if options == "alternative_platypus":
	
	gene_file = f'{working_directory}/Alignments/ZCWPW2_alignment_platypus.fas'
	tree_file = f'{working_directory}/Trees/alternative_ZCWPW2_platypus.nwk'
	out_file = f'{working_directory}/Results/results_alternative_zcwpw2_platypus.txt'
	#out_file = f'{working_directory}/Results/results_alternative_zcwpw2_platypus_fixed_omega_1.txt'
	
	options = f'{working_directory}/Models/ctl_branchmodel2_alternative_zcwpw2.txt'
	#options = f'{working_directory}/Models/ctl_branchmodel2_alternative_zcwpw2_fixed_omega_branch.txt'

	# Add the null tree in here
	cml = codeml.Codeml(alignment = gene_file, tree = tree_file,
                    out_file = out_file, working_dir=working_dir)

	cml.read_ctl_file(options)
	cml.print_options()

	cml.run(verbose = True)