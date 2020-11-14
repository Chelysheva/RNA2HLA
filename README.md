# RNA2HLA

HLA-based quality control of RNA-seq datasets

## Synopsis

Tool extracts the HLA types of I and II classes from all the files in the folder containing raw RNA-seq data (paired- on single-end).
The alleles are then cross-compared between the RNA-seq samples to identify the common source of the samples based on HLA types (4 digital resolution).


#### Release: 1.0

#### Author 
Irina Chelysheva, 2019-2020 (c)\
Oxford Vaccine Group, Department of Paediatrics, University of Oxford\
[Contact](irina.chelysheva@paediatrics.ox.ac.uk)

## Usage
```python RNA2HLA.py -f /raw_RNAseq_data_folder [-r /global_name_of_run] [-p <int>] [-3 <int>] [-c <float>] [-g <int>]```\
\
```-f``` is required for running RNA2HLA. Folder should contain raw RNA-seq samples, single- or paired-end or both types, in a compressed or not compressed formats.\
\
Optional parameters:
- ```-r``` to be used as a prefix for all output files
- ```-p``` number of parallel search threads for bowtie optional (default: 6)
- ```-3``` trim <int> bases from the low-quality end of each read
- ```-c``` confidence level for HLA-typing (default: 0.05)
- ```-g``` number of HLA genes to be included for typing (default: 5, may be increased to 6 - adding DQB1)


## Dependencies

1) RNA2HLA is a Python script (available in two versions: for Python 2 and Python 3 (coming soon)).
2) Dependent Python scripts single_end.py and paired_end.py and function scripts in R and Python must be located in the same folder.
3) [bowtie](http://bowtie-bio.sourceforge.net/index.shtml) must be reachable by the command ```bowtie```.
4) R must be installed.
5) Index files must be located in subfolder /references.
6) Packages: [biopython](https://github.com/biopython/biopython) (developed with V1.58), [numpy](https://github.com/numpy/numpy) (≥1.3.0).

## Output
The final output - overall comparison matrix in csv format, which crosscompares all RNA-seq samples in the folder.

Individual outputs in txt format produced for each RNA-seq sample in the folder (classes I and II are written in one file):
1) <sampleID>.bowtielog.txt - file with statistics of HLA mapping; 
2) <sampleID>.ambiguity.txt - reports typing ambuigities (if more than one solution for an allele possible);
3) <sampleID>.expression.txt - RPKM expression of HLA;
4) <sampleID>.HLAgenotype4digits.txt - 4 digital HLA type.
  
## Version history
1.0: initial tool

## References
Boegel S, Löwer M, Schäfer M, et al. HLA typing from RNA-Seq sequence reads. *Genome Med.* 2012;4(12):102. Published 2012 Dec 22. doi:10.1186/gm403

## License
[MIT](https://choosealicense.com/licenses/mit/)
