# bmDA

The implementation of benchmarking the performance of differential abundance (DA) testing methods including 5 clustering-free methods:

1. [Testing for differential abundance in mass cytometry data (Cydar)](https://www.nature.com/articles/nmeth.4295)
2. [Detection of differentially abundant cell subpopulations in scRNA-seq data (DAseq)](https://www.pnas.org/doi/abs/10.1073/pnas.2100293118)
3. [Quantifying the effect of experimental perturbations at single-cell resolution (MELD)](https://www.nature.com/articles/s41587-020-00803-5)
4. [Differential abundance testing on single-cell data using k-nearest neighbor graphs (Milo)](https://www.nature.com/articles/s41587-021-01033-z)
5. [Co-varying neighborhood analysis identifies cell populations associated with phenotypes of interest from single-cell transcriptomics (CNA)](https://www.nature.com/articles/s41587-021-01066-4),

and a clustering-based method, *Louvain*.

## Dependencies

To run this benchmarking codes, it needs to install a list of R and Python packages. The R packages needed are:

- argparse
- SingleCellExperiment
- scran
- [DAseq](https://github.com/KlugerLab/DAseq)
- [miloR](https://github.com/MarioniLab/miloR)
- tibble
- dplyr
- tidyverse
- igraph
- [cydar](http://bioconductor.org/packages/cydar)
- pdist
- reshape2

The Python packages needed are

- [MELD](https://github.com/KrishnaswamyLab/MELD)
- [cna](https://github.com/immunogenomics/cna)
- scanpy
- [graphtools](https://github.com/KrishnaswamyLab/graphtools)
- scikit-learn
- multianndata

## Data

- Synthetic datasets are available under the `data` directory.
- Real dataset are available at this [link](https://drive.google.com/drive/folders/15wWFD5FMe0VdzN1pUnaUUpQ17OXkeebH?usp=sharing).

## Usage

*Note:* Our implementation can only be used on a cluster with Slurm job scheduler since we need to run thousands of jobs.

The benchmarking scripts are all located in the `bin` drectory.

```text
bin
├── bm_parameter.sh
├── bm_runtime.sh
├── bm_syn_real.sh
└── make_bm_data.sh
```

To run a benchmarking job, use the following command:

```sh
bash bm_{the script}.sh
```


## Acknowledgement

Our implementation is inspired by the repo https://github.com/MarioniLab/milo_analysis_2020.
