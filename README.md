# BAD: Bold Assertor of Difficulty
[![image](https://img.shields.io/pypi/v/bad-phylo.svg)](https://pypi.python.org/pypi/bad-phylo)
[![image](https://img.shields.io/conda/vn/conda-forge/bad-phylo.svg)](https://anaconda.org/conda-forge/bad-phylo)
[![image](https://img.shields.io/badge/License-GPL3-yellow.svg)](https://opensource.org/licenses/GPL-3-0)

## Description

BAD is a python tool for predicting the difficulty of phylogenetic placements. BAD uses RAxML-NG outputs as input. It requires a RAxML-NG installation.
It was trained on empirical datasets from TreeBASE and can use both AA and DNA data. The output of BAD is a score between 0 (easy) and 1 (difficult). 
BAD provides an explanation of its prediction using the Shapley values implementation SHAP ([Github](https://github.com/shap/shap), [Paper](https://proceedings.neurips.cc/paper_files/paper/2017/file/8a20a8621978632d76c43dfd28b67767-Paper.pdf)).


## Installation
### Using pip
```
pip install bad-phylo
```
## Usage Example
A simple command line call of BAD looks like this:
```
bad -msa /test/example.fasta -tree /test/example.bestTree -model /test/example.bestModel -query /test/query.fasta -o test_bad 
```
This command will use the MSA and query file in fasta format, and the best tree inferred with RAxML-NG as well as the model.
It will compute features from all four data sources and predict the placement difficulties for each taxon in the query file.
All output files will be stored in an output folder called *test_bad* in the current directory. 
BAD will summarize the explanations for the prediction in the command line. For further details, please look at the SHAP summary plots or the *bad.log* file in the output folder.

Before interpreting the explanations provided by BAD, please make sure you know how to properly interpret Shapley values.
Easy to understand introduction to Shapley values: [https://christophm.github.io/interpretable-ml-book/shapley.html](https://christophm.github.io/interpretable-ml-book/shapley.html)

Please keep in mind that BAD requires an installation of RAxML-NG. By default, it uses the command ```raxml-ng```. 
If your RAxML-NG installation is not part of the PATH variable, you can specify the path to the RAxML-NG binary file with the parameter ```-raxmlng PATH_TO_RAXMLNG```.
### References
* A. M. Kozlov, D. Darriba, T. Flouri, B. Morel, and A. Stamatakis (2019) 
**RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference** 
*Bioinformatics*, 35(21): 4453–4455. 
[https://doi.org/10.1093/bioinformatics/btz305](https://doi.org/10.1093/bioinformatics/btz305)

* S. M. Lundberg and S.-I. Lee. A unified approach to interpreting model predictions. In
Proceedings of the 31st International Conference on Neural Information Processing Systems,
NIPS’17, page 4768–4777. Curran Associates Inc., 2017. ISBN 9781510860964. [https://proceedings.neurips.cc/paper_files/paper/2017/file/8a20a8621978632d76c43dfd28b67767-Paper.pdf](https://proceedings.neurips.cc/paper_files/paper/2017/file/8a20a8621978632d76c43dfd28b67767-Paper.pdf)
