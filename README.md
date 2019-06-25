# binomialRF
binomialRF R package

The binomialRF package is an R package that provides a feature selection algorithm to be used in randomForest classifiers. Treating each tree as a binomial stochastic process in a random forest, binomialRF determines a feature’s importance by how often they are selected in practice vs. as expected by random chance.


# Installing from GitHub

To install from GitHub directly, follow the code instructions below!

```
install.packages("devtools")


# The following dependencies might need to be installed
# manually if they're not installed by devtools. 

install.packages(c("ggplot2", "randomForest", "data.table","rlist"))
devtools::install_github("SamirRachidZaim/binomialRF")
library(binomialRF)
```

# Installing from BioConductor

Alternatively, the binomialRF R package has been submitted to BioConductor, and, upon availability, you will be able to install as follows: 

```
BiocManager::install("binomialRF")

```
