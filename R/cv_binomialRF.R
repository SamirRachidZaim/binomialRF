#' random forest feature selection based on binomial exact test
#'
#' \code{cv.binomialRF} is the cross-validated form of the \code{binomialRF}, where K-fold crossvalidation is conducted to assess the feature's significance. Using the \code{cvFolds}=K parameter, will result in a K-fold cross-validation where the data is 'chunked' into K-equally sized groups and then the averaged result is returned.
#' 
#' @param X design matrix
#' @param y class label
#' @param cvFolds how many times should we perform cross-validation
#' @param fdr.threshold fdr.threshold for determining which set of features are significant
#' @param fdr.method how should we adjust for multiple comparisons (i.e., \code{p.adjust.methods} =c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"))
#' @param ntrees how many trees should be used to grow the \code{randomForest}? (Defaults to 5000)
#' @param keep.both should we keep the naive binomialRF as well as the correlated adjustment
#'
#' @references Zaim, SZ; Kenost, C.; Lussier, YA; Zhang, HH. binomialRF: Scalable Feature Selection and Screening for Random Forests to Identify Biomarkers and Their Interactions, bioRxiv, 2019.
#'
#' @return a data.frame with 4 columns: Feature Name, cross-validated average for Frequency Selected, CV Median (Probability of Selecting it randomly), CV Median(Adjusted P-value based on \code{fdr.method}), and averaged number of times selected as signficant.
#'
#' @examples
#' set.seed(324)
#'
#' ###############################
#' ### Generate simulation data
#' ###############################
#'
#' X = matrix(rnorm(1000), ncol=10)
#' trueBeta= c(rep(10,5), rep(0,5))
#' z = 1 + X %*% trueBeta
#' pr = 1/(1+exp(-z))
#' y = as.factor(rbinom(100,1,pr))
#'
#' ###############################
#' ### Run cross-validation
#' ###############################
#'

.cv_binomialRF <- function(X,y, cvFolds=5, fdr.threshold=.05,  fdr.method='BY', ntrees=2000, keep.both =FALSE){
  requireNamespace('randomForest')
  requireNamespace('data.table')
  requireNamespace('stats')

  if(!is.numeric(ntrees)  | !is.numeric(fdr.threshold)){
    stop("Error: threshold, ntrees, and percent_features should be numeric inputs")
  } else if(ntrees <2){
    stop('ntrees must be a positive integer >1')
  } else if(!fdr.method %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")){
    stop('Please select acceptable fdr method from ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")')
  } else if(!is.logical(keep.both)){
    stop('keep.both must be a boolean value. Set to T or F')
  } else if(fdr.threshold >1 | fdr.threshold <0){
    stop("fdr.threshold is outside the acceptable (0-1) range")
  }

  if(!is.data.frame(X)){
    X = data.frame(X)
  }

  chunks = nrow(X)/cvFolds
  percent_features= seq(0.1, 1, length.out    = cvFolds)

  cv.bigMat = sapply(1:cvFolds, function(i) max(binomialRF(X[(((i-1)*chunks)+1): ((i)*chunks),],factor(y[(((i-1)*chunks)+1): ((i)*chunks)]),
                                           fdr.threshold, fdr.method, ntrees, percent_features[i] , keep.both = FALSE)))

  # cv.best.index = which.max(cv.bigMat)
  # cv.bigMat = data.frame(OOB.Error = 1/cv.bigMat)
  # cv.bigMat$PercentFeatures = percent_features

  # cv.bigMat$cv.Pvalue[cv.bigMat$cv.Pvalue==0] <- '<.0001'
  # cv.bigMat$cv.AdjPvalue[cv.bigMat$cv.AdjPvalue==0] <- '<.0001'

  return(cv.bigMat)

}

