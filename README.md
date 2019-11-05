# binomialRF
binomialRF R package

The binomialRF package is an R package that provides a feature selection algorithm to be used in randomForest classifiers. Treating each tree as a binomial stochastic process in a random forest, binomialRF determines a feature’s importance by how often they are selected in practice vs. as expected by random chance.


# Installing from GitHub

To install from GitHub directly, follow the code instructions below!

```
install.packages("devtools")


# The following dependencies might need to be installed
# manually if they're not installed by devtools. 

install.packages(c("ggplot2", "randomForest", "data.table","rlist", "correlbinom"))
devtools::install_github("SamirRachidZaim/binomialRF")
library(binomialRF)
```

# Installing from CRAN

Alternatively, the binomialRF R package has been submitted to CRAN, and, upon availability, you will be able to install as follows: 

```
install.packages('binomialRF')
```

# References: 

The main manuscript is included as a preprint in bioRxiv: https://doi.org/10.1101/681973, and has also been submitted for consideration at Frontiers in Genetics. 

