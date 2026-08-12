### myoRDNA
A toolkit for analysing rDNA methylation data from muscle and used in unison with methylKit.

## Installation
Given that this package is currently under development installation must performed using devtools.
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

## Functionality
```
absolute_rDNA_CN <- function(
    Meth_object,
    num_samples,
    autosomes = TRUE,
    rDNA_chr,
    diploid = TRUE)
```
 This function calculates the **absolute rDNA copy number** for each sample in a `Meth_object` (e.g., a list of methylation data or coverage data).

 ```
 relative_rDNA_CN <- function(
    Meth_object,
    num_samples,
    autosomes = TRUE,
    rDNA_chr,
    diploid = TRUE)
 ```
 This function calculates the **relative rDNA copy number** for each sample in a `Meth_object` (e.g., a list of methylation data or coverage data).

 ```
 get_unique_chr_names <- function(
    Meth_object,
    num_samples)
 ```
 This function extracts the unique chromosome names from the first sample in a `Meth_object` (e.g., a list of methylation data) and checks if all samples have the same chromosomes.
