# *Bradyrhizobium ottawaense*: How pan-genomic research has suddenly reveal strain with freakish host specific compatibility
<br>`выходные данные статьи`
<br><ссылка>
<p align="center">
  <img width="1200" height="900" src="https://github.com/blood097/Scientific_materials/blob/8cfd4c6360cda5c359535f348f1414eacabf3675/Pan-genome_Bradyrhizobium/abstract.png">
</p>

**Data**
<br>`Gene Presence Absence.txt` - gene presence and absence .Rtab file lists each gene and which samples it is present in
<br>`Number of Conserved Genes.txt` - .Rtab file on how the core genome varies as genomes are added (in random orders)
<br>`Number of Genes in Pan Geneome.txt` - .Rtab file on how the pan-genome varies as genomes are added (in random orders)
<br>`FASTA-to-Tabular pan-genome reference.tabular` - .tabular file contains a single representative nucleotide sequence from each of the clusters in the pan-genome
<br>`FastANI.tabular` - .tabular file contains query genome, reference genome, ANI value, count of bidirectional fragment mappings, and total query fragments
<br>`parsnp.nwk` - Newick format tree created using the binary presence and absence of accessory genes
<br>`IQ-TREE_Report.iqtree` -  IQ-TREE report file with full result of the run
<br>`emapper.annotations.xlsx` -  results from eggNOG mapper annotation phase in .xlsx format

**Analysis**
<br>`Pan-genome_post_roary.ipynb` - Jupyter Notebook for ploting presence/absence matrix against the tree, pangenome pie chart and mathematically evaluation of core genome size and pan-genome structure
<br>`FastANI_plotter.ipynb` - Jupyter Notebook for FastANI data analysis and visualisation
<br>`Roary_unique_gene_COG.ipynb` - Jupyter Notebook for COG category summary visualisation
<br>`rhierbaps.R` - R script for implementation of hierBAPS algorithm according to Tonkin-Hill *et al*. 2018

**Accessory**
<br>`COG one letter code descriptions.xlsx` - guide to COG one-letter codes categorize protein functions into main groups
___
Created by Kirichek Evgenii - e.kirichek@arriam.ru
