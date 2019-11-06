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
#' @param user_cbinom_dist insert either a pre-specified correlated binomial distribution or calculate one via the R package \code{correlbinom}.
#'
#' @references Zaim, SZ; Kenost, C.; Lussier, YA; Zhang, HH. binomialRF: Scalable Feature Selection and Screening for Random Forests to Identify Biomarkers and Their Interactions, bioRxiv, 2019.
#'
#' @return a data.frame with 4 columns: Feature Name, Frequency Selected, Probability of Selecting it randomly, Adjusted P-value based on \code{fdr.method}
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
#' ### Run binomialRF
#' ###############################
#'
#' binom.rf <-binomialRF(X,factor(y), fdr.threshold = .05,
#'                       ntrees = 2000,percent_features = .5,
#'                       fdr.method = 'BY',user_cbinom_dist=NULL)
#'
#' print(binom.rf)
#' @export


binomialRF <- function(X,y , fdr.threshold=.05, fdr.method='BY', ntrees=2000, percent_features=.5, keep.rf =FALSE, user_cbinom_dist=NULL){

  if(!is.numeric(ntrees)  | !is.numeric(percent_features)| !is.numeric(fdr.threshold)){
    stop("Error: threshold, ntrees, and percent_features should be numeric inputs")
  } else if( percent_features >1 | percent_features <0){
    stop("percent_features is outside the acceptable (0-1) range")
  } else if(ntrees <2){
    stop('L must be a positive integer >1')
  } else if(!fdr.method %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")){
    stop('Please select acceptable fdr method from ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")')
  } else if(!is.logical(keep.rf)){
    stop('keep.rf must be a boolean value. Set to T or F')
  } else if(fdr.threshold >1 | fdr.threshold <0){
    stop("fdr.threshold is outside the acceptable (0-1) range")
  }
  
  if(is.null(user_cbinom_dist)){
    if(ntrees==500 & ncol(X)==10){
      cbinom_dist = pmf_list$prob0.1$pmf_N500_Rho63
    } else if(ntrees==1000 & ncol(X)==10){
      cbinom_dist = pmf_list$prob0.1$pmf_N1000_Rho63
    } else if(ntrees==2000 & ncol(X)==10){
      cbinom_dist = pmf_list$prob0.1$pmf_N2000_Rho63
    } else if(ntrees==500 & ncol(X)==100){
      cbinom_dist = pmf_list$prob0.1$pmf_N500_Rho63
    } else if(ntrees==1000 & ncol(X)==100){
      cbinom_dist = pmf_list$prob0.1$pmf_N1000_Rho63
    } else if(ntrees==2000 & ncol(X)==100){
      cbinom_dist = pmf_list$prob0.1$pmf_N2000_Rho63
    } else if(ntrees==500 & ncol(X)==1000){
      cbinom_dist = pmf_list$prob0.1$pmf_N500_Rho63
    } else if(ntrees==1000 & ncol(X)==1000){
      cbinom_dist = pmf_list$prob0.1$pmf_N1000_Rho63
    } else if(ntrees==2000 & ncol(X)==1000){
      cbinom_dist = pmf_list$prob0.1$pmf_N2000_Rho63
    } else {
      stop('The correlated binomial distribution is outside of the pre-specified parameters. Please re-calculate it using the correlbinom R package.')
    }
    
  } else {
    cbinom_dist = user_cbinom_dist
  }

  if(!is.data.frame(X)){
    X = data.frame(X)
  }

  nms <- names(X)
  L = ncol(X)

  if(percent_features < 1/L){
    percent_features = 2/L
  }

  m = ceiling(percent_features*ncol(X)[1])
  
  ### need only grow 2 terminal nodes 
  ## for main effects to speed up computation
  
  rf.object <- randomForest(X,y, ntree = ntrees, mtry=m, keep.forest = TRUE, keep.inbag = TRUE,  replace=F , maxnodes = 2 )

  p = calculateBinomialP(L,percent_features )

  ## obtain root-node splitting variables
  main.effects <- data.table(rf.object$forest$bestvar[1,])
  main.effects[ , value:=T]
  main.effects[ , row.num:= 1:ntrees]
  
  rc.main.effects <- t(dcast(main.effects, idvar = 'V1', V1 ~ row.num))
  rc.main.effects <- data.frame(rc.main.effects[-1, ])
  rc.main.effects[is.na(rc.main.effects)] <- 0
  
  #### GIVES YOU Negative Log-Likelihood 
  #### for probaiblity ditrsitrubion
  
  cor.binomRF = data.table(melt(rc.main.effects))
  cor.binomRF = cor.binomRF[, list(freq=sum(value)), by='variable']
  
  pmf <- cbinom_dist/ sum(cbinom_dist)
  round(pmf,6)
  cmf <- cumsum(pmf)
  cor.binomRF$significance <- 1-cmf[cor.binomRF$freq]
  cor.binomRF$adjSignificance <- p.adjust(cor.binomRF$significance, method = fdr.method)
  
  cor.binomRF$Variable <- as.character(cor.binomRF$variable)
  
  cor.binomRF <- cor.binomRF[order(cor.binomRF$freq, decreasing=T)]
  
    if(keep.rf){
      return(list(binomRF = cor.binomRF,
                  RF.object = rf.object,
                  OutOfBagError = colMeans(rf.object$err.rate)[1]))
    } else{
     return(FeatureSelection=cor.binomRF)
    }
}
