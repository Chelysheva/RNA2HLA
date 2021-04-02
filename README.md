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
```$ python RNA2HLA.py -f /raw_RNAseq_data_folder [-r /global_name_of_run] [-p <int>] [-3 <int>] [-c <float>] [-g <int>]```\
\
```-f``` is required for running RNA2HLA. Folder should contain raw RNA-seq samples, single- or paired-end or both types, in a compressed or not compressed formats.\
\
Optional parameters:
- ```-r``` to be used as a prefix for all output files
- ```-p``` number of parallel search threads for bowtie (default: 6)
- ```-3``` trim <int> bases from the low-quality end of each read
- ```-c``` confidence level for HLA-typing (default: 0.05)
- ```-g``` number of HLA genes to be included for typing (default: 5, may be increased to 6 - adding DQB1)


## Dependencies

1) RNA2HLA is a Python script (available in two versions: for Python 2 and Python 3 (coming soon)).
2) All the dependencies provided within RNA2HLA depository (Python scripts single_end.py and paired_end.py, function scripts in R and Python, HLA class I and II databases) must be downloaded and located in the same folder.
3) Index files must be downloaded and located in subfolder /references.
4) Ther easiest way to run RNA2HLA is to create a [conda](https://github.com/conda/conda) environment using RNA2HLA_env.yml file provided:\
```$ conda env create -f RNA2HLA_env.yml``` \
And activate it:\
```$ source activate RNA2HLA_env``` or ```$ conda activate RNA2HLA_env``` (depends on the conda version)

**Update from 2.04.2021:**
One user reported an error while trying to create an environment from the original yml file (this error does not appear in most cases). If you experience an error, please, use an alternative environment file RNA2HLA_env_alt.yml instead.

Otherwise:\
4a) [bowtie](http://bowtie-bio.sourceforge.net/index.shtml) must be reachable by the command ```bowtie``` (developed with version 1.1.2)\
4b) R must be installed.\
4c) Packages: [biopython](https://github.com/biopython/biopython) (developed with 1.76), [numpy](https://github.com/numpy/numpy) (developed with 1.16.6, !caused an error for some users, in those cases - 1.15 is preferable), [pandas](https://github.com/pandas-dev/pandas) (developed with 0.24.2)

## Output
The final output - overall comparison matrix in csv format, which cross-compares all RNA-seq samples in the given folder.

Individual outputs in txt format produced for each RNA-seq sample in the folder (classes I and II are written in one file):
1) <sampleID>.bowtielog.txt - file with statistics of HLA mapping; 
2) <sampleID>.ambiguity.txt - reports typing ambuigities (if more than one solution for an allele possible based on the expression and HLA databases);
3) <sampleID>.expression.txt - RPKM expression of HLA;
4) <sampleID>.HLAgenotype4digits.txt - 4 digital HLA type.

## Limitation
In the case of studying a particular population with prior knowledge of the low HLA allele diversity, RNA2HLA should not be used as a QC, but only as a convenient study-wide HLA-typing method. One can refer to the [Allele Frequency Net Database](http://allelefrequencies.net/pop6001a_gsb.asp) and discover HLA diversity of particular population through the interactive map. The populations with less than 50 of total known alleles should be considered as of low diversity.

## Version history
1.0: initial tool

## Citations - RNA2HLA
Please, cite the following publication, if you are using RNA2HLA in your research: 
[Irina Chelysheva, Andrew J Pollard, Daniel O’Connor, RNA2HLA: HLA-based quality control of RNA-seq datasets, Briefings in Bioinformatics, 2021](https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bbab055/6184409)

## License
[MIT](https://choosealicense.com/licenses/mit/)
