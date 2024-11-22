# AllocateMate #

AllocateMate (Version 2) is a R package (R Core Team 2016) comprised of two primary functions: allocate.mate.ped and allocate.mate.H

AllocateMate (Version 1) is available at https://github.com/mghamilton/AllocateMateV1Archive

### What is AllocateMate ###

* AllocateMate generates a mating list for a set of parents.
* The mating list can be generated i) to minimise the average inbreeding coefficient (F) of families generated or ii) according to assortative mating principles.

* AllocateMate requires the following CRAN packages: lpSolveAPI, dplyr, AGHmatrix

### To install AllocateMate in R ###

*   install.packages(c("lpSolveAPI", "dplyr", "AGHmatrix"))

*   install.packages("devtools")
*   library(devtools)
*   install_github("mghamilton/AllocateMate")
*   library(AllocateMate)
*   help ("allocate.mate.ped")
*   help ("allocate.mate.H")

### Contact details ###

* <matthewhamilton75@gmail.com>

### References ###

* R Core Team (2016) R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria
