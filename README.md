### myoRDNA
`myoRDNA` is an R package developed by MuscleOmicsLab
to facilitate the analysis of rDNA copy number and DNA
methylation using high-throughput sequencing data.

The package provides a collection of functions for extracting
rDNA-specific data, estimating copy number, calculating
methylation levels and identifying differential methylation
between experimental conditions.

## Installation
Given that this package is currently under development installation must 
performed using devtools.
Please implement the following commands in the R terminal:
```
install.packages("devtools")
```
Once installation has completed
```
# Load devtools
library(devtools)

# Install the package from GitHub
install_github("MuscleOmicsLab/myoRDNA", branch = "main")

# Load the package
library(myoRDNA)
```

## Link to Useful Resources

# Manuscript highlighting examples of myoRDNA
```
Vaughan, D., Wood, N. and Seaborne, R.A., 2026. The ribosomal DNA landscape of 
mammalian muscle during acute and chronic physiological stress. bioRxiv, 
pp.2026-08.

Bibtex for latex users:

@article{vaughan2026ribosomal,
  title={The ribosomal DNA landscape of mammalian muscle during acute and 
  chronic physiological stress},
  author={Vaughan, Daniel and Wood, Nathanael and Seaborne, Robert AE},
  journal={bioRxiv},
  pages={2026--08},
  year={2026},
  publisher={Cold Spring Harbor Laboratory}
}
```

# Data used in the user guide and vignette:
```
Oyabu M, Ohira Y, Fujita M, Yoshioka K et al. Dnmt3a overexpression disrupts 
skeletal muscle homeostasis, promotes an aging-like phenotype, and reduces 
metabolic elasticity. iScience 2025 Apr 18;28(4):112144. 
```
[GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE262342)

## Contributing
We welcome contributions in any form: suggestions, issues, bugfixes. 
Pull requests should be made to the development branch.
