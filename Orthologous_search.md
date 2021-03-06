Blast candidates against each species protein sequence
-----------------------------
The blast exercise here is similar to that for detecting PRDM9 in the de novo transcriptomes. The script is very similar and is found here:
* [blasting_prdm9_candidates.py](https://github.com/izabelcavassim/PRDM9_analyses/blob/master/scripts/blasting_prdm9_candidates.py)

This script would give dataframes of the blast result for each species. The following script would build one presence and absence matrix for all species and genes.
* [combining_prdm9_blast_results.py](https://github.com/izabelcavassim/PRDM9_analyses/blob/master/scripts/combining_prdm9_blast_results.py)


Identification of the domain structure for each protein sequence
-----------------------------
To characterize the domain architecture for each sequence in each species we made use of the Conserved Domain Database(Marchler-Bauer et al., 2005). There are two ways one can submit their fasta sequences to the database. 
* Through python request library. The following code exemplifies how one could do it: 

[CDD_submission.py](https://github.com/izabelcavassim/PRDM9_analyses/blob/master/scripts/CDD_submission.py)

* Or through the CDD [website](https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi) using their Batch CD-Search.  
The Batch CD-Search accepts only protein sequences. The maximal number of queries per request is 4000.
The following code would split your concatenated fasta file (with all the sequences you want to describe domains: all_genes_combined.fasta) in subset fasta files of 3900 size:
``` bash
awk -v size=3900 -v pre=Batch_candidates -v pad=5 '
   /^>/ { n++; if (n % size == 1) { close(fname); fname = sprintf("%s.%0" pad "d", pre, n) } }
   { print >> fname }
' all_genes_combined.fasta
```
You would need to submit each of your subsets manually in the website. 

Building the phylogenetic tree
-----------------------------
To evaluate evidence for co-evolution between PRDM9 and the meiosis related genes among vertebrates, we need a high quality phylogenetic tree. Evolutionary information on species divergence times and relationships were obtained using the TimeTree resource [http://timetree.org/; Kumar et al., 2017, accessed date]. 32 out of the 379 species we considered were not present in this database. For these species, we used information from a close relative present in the database to determine their appropriate placement in the phylogenetic tree as follow:
```
Cebus imitator ----------> Cebus capucinus
Rhincodon typus (replaced with Ginglymostoma cirratum)
Paramormyrops kingsleyae (replaced with Paramormyrops gabonensis)
Oreochromis aureus (replaced with Oreochromis tanganicae)
Cynoglossus semilaevis (replaced with Cynoglossus lingua)
Paralichthys olivaceus (replaced with Paralichthys dentatus)
Etheostoma spectabile (replaced with Etheostoma caeruleum)
Etheostoma cragini (replaced with Etheostoma pallididorsum)
Notolabrus celidotus (replaced with Notolabrus gymnogenis)
Periophthalmus magnuspinnatus (replaced with Periophthalmus argentilineatus)
Thalassophryne amazonica (replaced with Opsanus tau)
Tachysurus fulvidraco (replaced with Leiocassis longirostris)
Astyanax mexicanus (replaced with Hollandichthys multifasciatus)
Microcaecilia unicolor (replaced with Microcaecilia sp. PZ-2009)
Python bivittatus (replaced with Python molurus)
Apteryx mantelli (replaced with Apteryx australis)
Chiroxiphia lanceolata (replaced with Chiroxiphia caudata)
Gopherus evgoodei (replaced with Gopherus agassizii)
Tupaia chinensis (replaced with Tupaia glis)
Corvus cornix (no substitute found) ------------> Corvus_woodfordi
Astatotilapia calliptera (no substitute found) ------------> Copadichromis_virginalis
Oreochromis niloticus (no substitute found) -------------> Pungu_maclareni
Strigops habroptila (no substitute found) ---------------> Probosciger_aterrimus
Tinamus guttatus (no substitute found) ----------------> Tinamus_major
Poecilia formosa (no substitute found) ----------------> Poecilia_butleri
Kryptolebias marmoratus (no substitute found) ----------------> Nothobranchius_virgatus
Oryzias melastigma (no substitute found) -------------------> Oryzias_luzonensis
Neophocaena asiaeorientalis (no substitute found) --------------> Phocoena_phocoena
Austrofundulus limnaeus (no substitute found) ----------------> Nothobranchius_fuscotaeniatus
Nematolebias whitei (no substitute found) -----------> Nothobranchius thierryi
Hippoglossus stenolepis (no substitute found) -----------------> Hippoglossoides_elassodon
Poecilia latipinna (no substitute found) --------------------> Poecilia_vivipara
```

To plot the phylogenetic tree we used the software itool https://itol.embl.de/. The final species tree with original names is found at: species_tree_379_classification_09_03_21.nwk

Prunning the zero branches
-----------------------------
Because polytomies (zero branches along the tree) can affect the phylogenetic independent constrast analyses, we pruned the tree such that we mantained the same number of PRDM9 transitions. The pruning script is found at: remove.zero.branches.R 

