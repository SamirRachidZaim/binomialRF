# binomialRF R Package

The binomialRF package is an R package that provides a feature selection algorithm to be used in randomForest classifiers. Treating each tree as a quasi binomial stochastic process in a random forest, binomialRF determines a feature’s importance by how often they are selected in practice vs. as expected by random chance. Given that trees are co-dependent as they subsample the same data, a theoretical adjustment is made using a generalization of the binomial distribution that adds a parameter to model correlation/association between trials. 

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

The main manuscript can be found on BMC Bioinformatics: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03718-9.
