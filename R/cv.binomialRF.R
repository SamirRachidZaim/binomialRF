#' random forest feature selection based on binomial exact test
#'
#' \code{cv.binomialRF} is the cross-validated form of the \code{binomialRF}, where K-fold crossvalidation is conducted to assess the feature's significance. Using the \code{cvFolds}=K parameter, will result in a K-fold cross-validation where the data is 'chunked' into K-equally sized groups and then the averaged result is returned.
#' @param X design matrix
#' @param y class label
#' @param cvFolds how many times should we perform cross-validation
#' @param fdr.threshold fdr.threshold for determining which set of features are significant
#' @param fdr.method how should we adjust for multiple comparisons (i.e., \code{p.adjust.methods} =c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"))
#' @param ntrees how many trees should be used to grow the \code{randomForest}? (Defaults to 5000)
#' @param percent_features what percentage of L do we subsample at each tree? Should be a proportion between (0,1)
#' @param keep.rf should we keep the randomForest object?
#'
#' @references Zaim, SZ. binomialRF: A novel randomForest feature selection algorithm. ArXiv, 2019.
#'
#' @return a data.frame with 4 columns: Feature Name, cross-validated average for Frequency Selected, CV Median (Probability of Selecting it randomly), CV Median(Adjusted P-value based on \code{fdr.method}), and averaged number of times selected as signficant.

cv.binomialRF <- function(X,y, cvFolds=5, fdr.threshold=.05, fdr.method='BH', ntrees=1000, percent_features=.5, keep.rf =F){
  require(randomForest)
  require(data.table)

  if(class(ntrees) != 'numeric' | class(percent_features)!= "numeric" | class(fdr.threshold)!= 'numeric'){
    stop("Error: threshold, ntrees, and percent_features should be numeric inputs")
  } else if( percent_features >1 | percent_features <0){
    stop("percent_features is outside the acceptable (0-1) range")
  } else if(ntrees <2){
    stop('L must be a positive integer >1')
  } else if(!fdr.method %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")){
    stop('Please select acceptable fdr method from ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")')
  } else if(class(keep.rf) != 'logical'){
    stop('keep.rf must be a boolean value. Set to T or F')
  } else if(fdr.threshold >1 | fdr.threshold <0){
    stop("fdr.threshold is outside the acceptable (0-1) range")
  }

  chunks = nrow(X)/cvFolds
  percent_features= seq(0.1, 1, length.out    = cvFolds)

  cv.bigMat = sapply(1:cvFolds, function(i) max(SMARTVIS::binomialRF(X[(((i-1)*chunks)+1): ((i)*chunks),],factor(y[(((i-1)*chunks)+1): ((i)*chunks)]),
                                           fdr.threshold, fdr.method, ntrees, percent_features[i] , keep.rf = F)$weight))

  cv.best.index = which.max(cv.bigMat)
  cv.bigMat = data.frame(OOB.Error = 1/cv.bigMat)
  cv.bigMat$PercentFeatures = percent_features

  # cv.bigMat$cv.Pvalue[cv.bigMat$cv.Pvalue==0] <- '<.0001'
  # cv.bigMat$cv.AdjPvalue[cv.bigMat$cv.AdjPvalue==0] <- '<.0001'

  return(cv.bigMat)

}

