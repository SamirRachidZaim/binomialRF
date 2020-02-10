# binomialRF
binomialRF R package

The binomialRF package is an R package that provides a feature selection algorithm to be used in randomForest classifiers. Treating each tree as a binomial stochastic process in a random forest, binomialRF determines a feature’s importance by how often they are selected in practice vs. as expected by random chance.

# Installing from CRAN

The binomialRF R package is on CRAN, and you can install as follows: 

```
install.packages('binomialRF')
```
The CRAN version will always be the most stable release. 

# Installing from GitHub

To install experimental updates from the binomialRF , install it from GitHub directly, follow the code instructions below!

```
install.packages("devtools")


# The following dependencies might need to be installed
# manually if they're not installed by devtools. 

install.packages(c("ggplot2", "randomForest", "data.table","rlist", "correlbinom"))
devtools::install_github("SamirRachidZaim/binomialRF")
library(binomialRF)
```

These GitHub updates and features are experimental and will not be available in the CRAN version until the next, stable release is pushed. 


# References: 

The main manuscript is included as a preprint in bioRxiv: https://doi.org/10.1101/681973, and has also been submitted for consideration at Frontiers in Genetics. 

