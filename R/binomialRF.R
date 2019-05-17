#' random forest feature selection based on binomial exact test
#'
#' \code{binomialRF} is the R implementation of the feature selection algorithm by (Zaim 2019)
#'
#' @param X design matrix
#' @param y class label
#' @param fdr.threshold fdr.threshold for determining which set of features are significant
#' @param fdr.method how should we adjust for multiple comparisons (i.e., \code{p.adjust.methods} =c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"))
#' @param ntrees how many trees should be used to grow the \code{randomForest}? (Defaults to 5000)
#' @param percent_features what percentage of L do we subsample at each tree? Should be a proportion between (0,1)
#' @param keep.rf should we keep the randomForest object?
#'
#' @references Zaim, SZ. binomialRF: A novel randomForest feature selection algorithm. ArXiv, 2019.
#'
#' @return a data.frame with 4 columns: Feature Name, Frequency Selected, Probability of Selecting it randomly, Adjusted P-value based on \code{fdr.method}

binomialRF <- function(X,y , fdr.threshold=.05, fdr.method='BH', ntrees=5000, percent_features=.2, keep.rf =F){
  require(randomForest)

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

  if(class(X)=='matrix'){
    X = data.frame(X)
  }

  nms <- names(X)
  L = ncol(X)

  if(percent_features < 1/L){
    percent_features = 2/L
  }

  m = ceiling(percent_features*ncol(X)[1])
  rf.object <- randomForest(X,y, ntree = ntrees, mtry=m, keep.forest = T, keep.inbag = T)

  p = calculateBinomialP(L,percent_features )

  ## obtain root-node splitting variables
  main.effects <- rf.object$forest$bestvar[1,]

  ## tabulate and calculate p-vals
  tabular_main.effects = as.data.frame(table(main.effects), stringsAsFactors = F)
  colnames(tabular_main.effects) <- c('Variable','Frequency')
  tabular_main.effects$Variable <- nms[as.numeric(tabular_main.effects$Variable)]

  tabular_main.effects$Pvalue <- as.numeric(sapply(tabular_main.effects$Freq, function(x) binom.test(x, n= ntrees, p, alternative='greater')$p.value))
  tabular_main.effects$AdjPvalue <- p.adjust(tabular_main.effects$Pvalue, method = fdr.method)
  tabular_main.effects$Significant <- tabular_main.effects$AdjPvalue < fdr.threshold
  tabular_main.effects = tabular_main.effects[order(tabular_main.effects$Freq, decreasing = T),]
  tabular_main.effects$weight = 1/colMeans(rf.object$err.rate)[1]

  if(keep.rf){
    return(list(binomRF = tabular_main.effects,
                RF.object = rf.object,
                OutOfBagError = colMeans(rf.object$err.rate)[1]))
  } else{
   return(FeatureSelection=tabular_main.effects)
  }
}
