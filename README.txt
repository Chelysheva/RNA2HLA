####################################################################
#
#Title:
#RNA2HLA: HLA-based quality control of RNA-seq datasets
#
#Synopsis:
#Tool extracts the HLA types of I and II classes from all the files in the folder containing raw RNA-seq data (paired on single-end).
#The alleles are then cross-compared between the RNA-seq samples to identify the common source of the samples based on HLA types (4 digital resolution).
#
#Release: 1.0
#
#Author: Irina Chelysheva, 2020 (c)
#Oxford Vaccine Group, Department of Paediatrics, University of Oxford
#Contact: irina.chelysheva@paediatrics.ox.ac.uk
#
#Usage: 
#python RNA2HLA.py -f /folder_with_raw_RNAseq_data [-r /global_name_of_run]* [-p <int>]** [-3 <int>]*** [-c <float>]**** [-g <int>]*****
#*optional: to be used as a prefix for all output files
#**optional: number of parallel search threads for bowtie optional (default: 6)
#***optional: trim int bases from the low-quality end of each read
#****optional: confidence level for HLA-typing (default: 0.05)
#*****optional: number of HLA genes to be included for typing (default: 5, maybe increased to 6 - adding DQB1)
#folder should contain raw RNA-seq samples, single- or paired-end or both types, in a compressed or not compressed formats
#
#Dependencies:
#1) RNA2HLA is a python script (available in two versions: Python 2 and Python 3).
#2) Dependent python scripts single_end.py and paired_end.py must be located in the same folder.
#3) bowtie must be reachable by the command "bowtie".
#4) R must be installed.
#5) Index files must be located in the folder "references".
#6) Packages: biopython (developed with V1.58), numpy (1.3.0).
#
#Output:
#The final output - overall comparison matrix in csv format, which crosscompares all RNA-seq samples in the folder.
#Individual outputs in txt format produced for each RNA-seq sample in the folder (classes I and II are written in one file): 
#1)<sampleID>.bowtielog.txt - file with statistics of HLA mapping; 
#2) <sampleID>.ambiguity.txt - reports typing ambuigities (if more than one solution for an allele possible);
#3) <sampleID>.expression.txt - RPKM expression of HLA;
#4) <sampleID>.HLAgenotype4digits.txt - 4 digital HLA type.
#
#Version history:
#1.0: initial tool
#
#References:
#Boegel, Sebastian; Loewer, Martin; Schaefer, Michael; Bukur, Thomas; Graaf, Jos de; Boisguerin, Valesca et al. (2013): HLA typing from RNA-Seq sequence reads. In: Genome Med 4 (12), S. 102. DOI: 10.1186/gm403.
#
#License:
#GNU General Public License v3.0
#Copyright (c) 2020 Irina Chelysheva

####################################################################

