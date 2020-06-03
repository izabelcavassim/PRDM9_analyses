# 1. Create a parser for the aminoacid sequences from UCSC (using maybe gwf?)
# 2. Run aaml in the protein sequences (should I select genes that are found in a certain percentage of the species? (80% of the species))
# 3. Run RERconverge
# 4. Interpret analysis and do enrichment analysis
# 5. Run codeml for the top candidates to see if the gene is under constraint relaxation of if the gene is under positive selection.

import glob
from collections import OrderedDict
import pandas as pd
from collections import OrderedDict
from sys import argv

alignment_directory = argv[1]
parsed_transcripts_directory = argv[2]
parsed_genes_directory = argv[3]


# 2. Parsing fasta files:
# Parsing the fasta files:
present_in_baker = ['Balaenoptera acutorostrata scammoni', 'Ovis aries', 'Python bivittatus', 'Pogona vitticeps', 'Echinops telfairi', 'Marmota marmota marmota', 'Ursus maritimus', 'Equus przewalskii', 'Myotis lucifugus', 'Camelus bactrianus', 'Dasypus novemcinctus', 'Poecilia formosa', 'Heterocephalus glaber', 'Loxodonta africana', 'Macaca mulatta', 'Nannospalax galili', 'Rousettus aegyptiacus', 'Poecilia reticulata', 'Macaca fascicularis', 'Camelus dromedarius', 'Tursiops truncatus', 'Latimeria chalumnae', 'Mandrillus leucophaeus', 'Bison bison bison', 'Xiphophorus maculatus', 'Gorilla gorilla gorilla', 'Chrysemys picta bellii', 'Fundulus heteroclitus', 'Rattus norvegicus', 'Gekko japonicus', 'Fukomys damarensis', 'Macaca nemestrina', 'Microtus ochrogaster', 'Microcebus murinus', 'Ochotona princeps', 'Colobus angolensis palliatus', 'Larimichthys crocea', 'Leptonychotes weddellii', 'Mus musculus', 'Pan troglodytes', 'Myotis davidii', 'Bubalus bubalis', 'Chinchilla lanigera', 'Felis catus', 'Protobothrops mucrosquamatus', 'Mustela putorius furo', 'Camelus ferus',
                    'Sorex araneus', 'Propithecus coquereli', 'Dipodomys ordii', 'Cricetulus griseus', 'Monodelphis domestica', 'Sus scrofa', 'Galeopterus variegatus', 'Octodon degus', 'Panthera tigris altaica', 'Ictidomys tridecemlineatus', 'Pan paniscus', 'Haplochromis burtoni', 'Pongo abelii', 'Equus caballus', 'Bos taurus', 'Lipotes vexillifer', 'Eptesicus fuscus', 'Oryzias latipes', 'Pteropus vampyrus', 'Sarcophilus harrisii', 'Tupaia chinensis', 'Ailuropoda melanoleuca', 'Cercocebus atys', 'Danio rerio', 'Erinaceus europaeus', 'Poecilia latipinna', 'Myotis brandtii', 'Vicugna pacos', 'Acinonyx jubatus', 'Orcinus orca', 'Otolemur garnettii', 'Poecilia mexicana', 'Miniopterus natalensis', 'Nomascus leucogenys', 'Pundamilia nyererei', 'Pteropus alecto', 'Orycteropus afer afer', 'Homo sapiens', 'Carlito syrichta', 'Bos mutus', 'Chlorocebus sabaeus', 'Rhinopithecus roxellana', 'Peromyscus maniculatus bairdii', 'Papio anubis', 'Ovis aries musimon', 'Mesocricetus auratus', 'Chrysochloris asiatica', 'Jaculus jaculus', 'Oryctolagus cuniculus', 'Ornithorhynchus anatinus']

def fasta_parser(filename, metadata_dict, tuple_in=False):
    file = open(filename, 'r').read()
    file_separe = file.split('>')
    try:
        file_separe.remove('')
    except:
        file_separe.remove('\n')
        pass

    missing_rate = 0
    dict_fasta = OrderedDict()
    for entry in file_separe:
        seq = entry.splitlines()  # list

        header = seq[0]

        # Stripping header, and leaving only the fasta name.
        header_list = header.split('_')
        #print(header_list)

        gene_ucsc_id = header_list[0]
        fasta_name = header_list[1]
        exon_position = header_list[2]

        # Also replacing the fasta name by the Scientific_tree_name based on the metadata file
        header = metadata_dict['Scientific_tree_name'][fasta_name]

        sequences = ''.join(seq[1::])

        # count the number of missing data per gene
        seq_len = len(sequences)
        miss_len = sequences.count('-')
        if seq_len == miss_len:
            missing_rate += 1
            dict_fasta[header] = sequences
        else:
            dict_fasta[header] = sequences
    print(dict_fasta)
    completeness = 1 - float(missing_rate)/100
    completeness = str(round(completeness, 3))
    gene_name = gene_ucsc_id + '_'  + exon_position + '_' + str(seq_len) + '_' + completeness
    if tuple_in == True:
        return((gene_name, dict_fasta))
    else:
        return(dict_fasta)

def parse_fasta_2(filename):
    file = open(filename, 'r').read() #opening and reading the fasta file, putting it in a object called file
    gene = filename.split("/")
    gene = gene[len(gene)-1]
    file_separe = file.split('>') #spliting each entry by the > 
    print(file_separe)
    file_separe.remove('')
    parse_dict = {}
    header = []
    for entry in file_separe:
        seq = entry.splitlines()
        header = seq[0] #these are the first elements of the list 
        seq = ''.join(seq[1:]) #joining the sequences 
        parse_dict[header] = seq

    print('This genes %s contain %d members' % (gene, len(parse_dict.keys())))
    return parse_dict

metadata = pd.read_csv('/home/mica16/MouseLemur/scripts/codeml/metadata2.csv', sep=";", header=0)
metadata = metadata.set_index('Fasta_name')
metadata_dict = metadata.to_dict()

def write_fasta(tuple_fasta, directory, counts):
    (fasta_name, dict_fasta) = tuple_fasta
    fasta_name = fasta_name + '_' + str(counts)+".fasta"

    with open(directory+"{}".format(fasta_name), 'w') as f:

        for header, sequences in dict_fasta.items():
            f.write('>{}\n'.format(header))
            chunks = [sequences[i:i+60] for i in range(0, len(sequences), 60)]
            for i in chunks:
                f.write('{}\n'.format(i))
        f.close()


# Finding all alignments:
def parsing_all_alignments(metadata_dict=metadata_dict):
    assemblies = [(f) for f in glob.glob(alignment_directory+"*.fasta")]
    #assemblies = [alignment_directory+"gene_1.fasta", alignment_directory+"gene_2.fasta", alignment_directory+"gene_3.fasta", alignment_directory+"gene_4.fasta", alignment_directory+"gene_5.fasta"]
    counts = 0
    for i in assemblies:
        print(i)
        test = fasta_parser(i, metadata_dict, tuple_in=True)
        counts += 1
        print(counts)
        write_fasta(test, directory=parsed_transcripts_directory, counts = counts)

parsing_all_alignments()

def getKey(item):
    return item[0]

def combining_exons():
     dict_parsed_files = [(f) for f in glob.glob(parsed_transcripts_directory+"*")]
     # creating a dictionary where keys are gene names and values are a list of tuples of exon position followed by the protein sequence
     combining_exons = OrderedDict()
     for i in dict_parsed_files:
         ucsc_id = i.split('/')
         ucsc_id_gene = ucsc_id[6]
         gene_id = ucsc_id_gene.split('_')[0]
         exon_number = ucsc_id_gene.split('_')[1]
         print(ucsc_id_gene)
         #print(exon_number)
         
         if gene_id not in combining_exons:
             combining_exons[gene_id] = {}
             combining_exons[gene_id][int(exon_number)] = parse_fasta_2(i)
         else:
             combining_exons[gene_id][int(exon_number)] = parse_fasta_2(i)

     print(combining_exons[gene_id].keys())
     #exons = list()
     for gene in combining_exons:
         d = {}
         for key, value in combining_exons[gene].items():
             for key2, value2 in value.items():
                 if key2 not in d:
                     d[key2] = [value2]
                 else:
                     d[key2].append(value2)
         with open(parsed_genes_directory+"{}.fasta".format(gene), 'w') as f:
             for header, sequences in d.items():
                 sequences = ''.join(sequences)
                 len_seq = len(sequences)
                 len_ns = sequences.count('-')
                 if len_ns <= (0.80)*len_seq:
                     f.write('>{}\n'.format(header))
                     #chunks = [sequences[i:i+60] for i in range(0, len(sequences), 60)]
                     #for i in chunks:
                     f.write('{}\n'.format(sequences))
                 else:
                     pass
     
combining_exons()
