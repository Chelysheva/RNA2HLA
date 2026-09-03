# RNA2HLA

HLA-based quality control of RNA-seq datasets

## Synopsis

Tool extracts the HLA types of I and II classes from all the files in the folder containing raw RNA-seq data (paired- or single-end).
The alleles are then cross-compared between the RNA-seq samples to identify the common source of the samples based on HLA types (4 digital resolution).

#### Releases: v1.0 - original tool; v1.1 - script creating heatmap output is added; v1.2 - Python 3.7 port with compatibility and reproducibility fixes

#### Author

Dr. Irina Chelysheva, 2019-2026 (c)\
Senior Lecturer in Human Genomics, Oxford Brookes University\
Oxford Vaccine Group, Department of Paediatrics, University of Oxford\
[Contact](i.chelysheva@brookes.ac.uk)

## Usage

```$ python RNA2HLA.py -f /raw_RNAseq_data_folder [-r /global_name_of_run] [-p <int>] [-3 <int>] [-c <float>] [-g <int>]```

Optional parameters:

- ```-r``` to be used as a prefix for all output files
- ```-p``` number of parallel search threads for bowtie (default: 6)
- ```-3``` trim bases from the low-quality end of each read
- ```-c``` confidence level for HLA-typing (default: 0.05)
- ```-g``` number of HLA genes to be included for typing (default: 5, may be increased to 6 - adding DQB1)

## Dependencies

1) RNA2HLA is a Python 3 script (v1.2, developed and validated with Python 3.7). The Python 2 version is archived in release v1.1.
2) All the dependencies provided within RNA2HLA repository (Python scripts single_end.py and paired_end.py, function scripts in R and Python, HLA class I and II databases) must be downloaded and located in the same folder.
3) Index files must be downloaded and located in subfolder /references.
4) The easiest way to run RNA2HLA is to create a [conda](https://github.com/conda/conda) environment using the RNA2HLA_env.yml file provided:

```$ conda env create -f RNA2HLA_env.yml```

And activate it:

```$ conda activate RNA2HLA_py37_env```

Otherwise:\
4a) [bowtie](http://bowtie-bio.sourceforge.net/index.shtml) (1.x) must be reachable by the command ```bowtie``` (validated with versions 1.1.2 and 1.3.1)\
4b) R must be installed (validated with R 4.4.3; only base stats functions are used).\
4c) Packages: [biopython](https://github.com/biopython/biopython) (validated with 1.79), [numpy](https://github.com/numpy/numpy) (validated with 1.21.6), [pandas](https://github.com/pandas-dev/pandas) (validated with 1.3.5)

## Output

The final output - overall comparison matrix in csv format, which cross-compares all RNA-seq samples in the given folder.

Individual outputs in txt format produced for each RNA-seq sample in the folder (classes I and II are written in one file):
1) .bowtielog.txt - file with statistics of HLA mapping;
2) .ambiguity.txt - reports typing ambiguities (if more than one solution for an allele possible based on the expression and HLA databases);
3) .expression.txt - RPKM expression of HLA;
4) .HLAgenotype4digits.txt - 4 digital HLA type.

A heatmap can be created from the overall comparison matrix csv file using the R script heatmap_HLA_identity_comparison.R

## What is new in v1.2

- **Python 3.7 port.** The tool was ported from Python 2.7 and validated side-by-side against the original on real RNA-seq data (6 whole-blood samples from 3 participants x 2 time points, SRP276081, plus a full-depth check): HLA calls, expression values, ambiguity reports and the sample-comparison matrix are reproduced; bowtie mapping logs are byte-identical between versions.
- **bowtie >= 1.2 compatibility fix.** Previous versions crash with bowtie >= 1.2 because the parser assumed a fixed first line of the bowtie log (bowtie 1.2+ prepends a deprecation warning). The log is now parsed by content, and the tool works with both bowtie 1.1.2 and 1.3.1 with identical results.
- **Deterministic tie-breaking.** When two allele groups receive exactly the same number of reads, the choice between them was arbitrary in previous versions (dependent on Python 2 dict hash order). Ties are now resolved deterministically (lexicographic order). Note: in the rare event of an exact read-count tie, v1.2 can therefore report a different (but reproducible) allele than v1.0/v1.1.
- **Confidence formatting** is numerically identical to the Python 2 output (12 significant digits).

## Limitation

In the case of studying a particular population with prior knowledge of the low HLA allele diversity, RNA2HLA should not be used as a QC, but only as a convenient study-wide HLA-typing method. One can refer to the [Allele Frequency Net Database](http://allelefrequencies.net/pop6001a_gsb.asp) and discover HLA diversity of particular population through the interactive map. The populations with less than 50 of total known alleles should be considered as of low diversity.

Note: HLA class II second-allele calls are sensitive to sequencing depth; low-confidence calls (p > 0.05) are excluded from the comparison matrix.

## Version history

1.0: initial tool\
1.1: script creating heatmap output is added\
1.2: Python 3.7 port; bowtie >= 1.2 compatibility fix; deterministic tie-breaking

## Citations - RNA2HLA

Please, cite the following publication, if you are using RNA2HLA in your research:
[Irina Chelysheva, Andrew J Pollard, Daniel O'Connor, RNA2HLA: HLA-based quality control of RNA-seq datasets, Briefings in Bioinformatics, 2021](https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bbab055/6184409)

## License

[MIT](https://choosealicense.com/licenses/mit/)
