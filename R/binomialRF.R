#' random forest feature selection based on binomial exact test
#'
#' \code{binomialRF} is the R implementation of the feature selection algorithm by (Zaim 2019)
#' 
#' @usage binomialRF(X,y, fdr.threshold = .05,fdr.method = 'BY',
#'                       ntrees = 2000, percent_features = .5,
#'                       keep.both=FALSE, user_cbinom_dist=NULL,
#'                       sampsize=round(nrow(X)*.63))
#'                       
#' @param X design matrix
#' @param y class label
#' @param fdr.threshold fdr.threshold for determining which set of features are significant
#' @param fdr.method how should we adjust for multiple comparisons (i.e., \code{p.adjust.methods} =c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"))
#' @param ntrees how many trees should be used to grow the \code{randomForest}? 
#' @param percent_features what percentage of L do we subsample at each tree? Should be a proportion between (0,1)
#' @param keep.both should we keep the naive binomialRF as well as the correlated adjustment
#' @param user_cbinom_dist insert either a pre-specified correlated binomial distribution or calculate one via the R package \code{correlbinom}.
#' @param sampsize how many samples should be included in each tree in the randomForest
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
#' rho = 0.33
#' ntrees = 500
#' cbinom = correlbinom::correlbinom(rho, successprob =  calculateBinomialP(10, .5), trials = ntrees, 
#'                                precision = 1024, model = 'kuk')
#'
#' binom.rf <-binomialRF(X,y, fdr.threshold = .05,fdr.method = 'BY',
#'                       ntrees = 500,percent_features = .5,
#'                       keep.both=FALSE, user_cbinom_dist=cbinom,
#'                       sampsize=round(nrow(X)*rho))
#'
#' print(binom.rf)
#' @export


binomialRF <- function(X,y , fdr.threshold=.05, fdr.method='BY', ntrees=2000, percent_features=.5, keep.both =FALSE,  user_cbinom_dist=NULL, sampsize=round(nrow(X)*.63)){

  pmf_list = binomialRF::pmf_list
  if(!is.numeric(ntrees)  | !is.numeric(percent_features)| !is.numeric(fdr.threshold)){
    stop("Error: threshold, ntrees, and percent_features should be numeric inputs")
  } else if( percent_features >1 | percent_features <0){
    stop("percent_features is outside the acceptable (0-1) range")
  } else if(ntrees <2){
    stop('L must be a positive integer >1')
  } else if(!fdr.method %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")){
    stop('Please select acceptable fdr method from ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")')
  } else if(!is.logical(keep.both)){
    stop('keep.both must be a boolean value. Set to T or F')
  } else if(fdr.threshold >1 | fdr.threshold <0){
    stop("fdr.threshold is outside the acceptable (0-1) range")
  } 
  
  if(is.null(user_cbinom_dist)){
    if(ntrees==2000 & ncol(X)==10){
      cbinom_dist = pmf_list$prob0.1$pmf_N2000_Rho63
    } else if(ntrees==2000 & ncol(X)==100){
      cbinom_dist = pmf_list$prob0.01$pmf_N2000_Rho63
    } else if(ntrees==2000 & ncol(X)==1000){
      cbinom_dist = pmf_list$prob0.001$pmf_N2000_Rho63
    } else {
      stop("Please specify a correlated binomial distribution. See the pmf_list object for some default options, otherwise refer to the correlbinom R package and recalculate a valid distribution.")
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
    percent_features = 0.3
  }

  m = ceiling(percent_features*ncol(X)[1])
  
  ### need only grow 2 terminal nodes 
  ## for main effects to speed up computation
  
  rf.object <- randomForest::randomForest(X,y, ntree = ntrees, mtry=m, keep.forest = TRUE, keep.inbag = TRUE,  replace=F , maxnodes = 2 , sampsize=sampsize)

  p = calculateBinomialP(L,percent_features )

  ## obtain root-node splitting variables
  main.effects <- data.table::data.table(rf.object$forest$bestvar[1,])
  main.effects$value <- T
  main.effects$row.num <- 1:ntrees

  ## obtain the matrix of trees / per-variable
  rc.main.effects <- t(data.table::dcast(main.effects, idvar = 'V1', V1 ~ row.num))
  col.names = nms[rc.main.effects[1,]]
  
  rc.main.effects <- data.frame(rc.main.effects[-1, ])
  colnames(rc.main.effects) <- col.names
  rc.main.effects[is.na(rc.main.effects)] <- 0
  
  #### GIVES YOU Negative Log-Likelihood 
  #### for probaiblity ditrsitrubion
  
  cor.binomRF = data.table::data.table(data.table::melt(rc.main.effects))
  cor.binomRF = data.frame(variable = dimnames(table(cor.binomRF))$variable,
                           freq = table(cor.binomRF)[,2])

  binomRF <- cor.binomRF
  
  ### calculate correlation adjusted binomial
  pmf <- cbinom_dist/ sum(cbinom_dist)
  round(pmf,6)
  cmf <- cumsum(pmf)
  cor.binomRF$significance <- 1-cmf[cor.binomRF$freq]
  cor.binomRF$adjSignificance <- stats::p.adjust(cor.binomRF$significance, method = fdr.method)
  cor.binomRF$variable <- as.character(cor.binomRF$variable)
    
  cor.binomRF <- cor.binomRF[order(cor.binomRF$freq, decreasing=T),]
  
  ### calculate naive binomial
  binomRF$significance <- sapply(binomRF$freq, function(zzz) stats::binom.test(zzz, n=ntrees, p=p, alternative = 'greater')$p.value)
  binomRF$adjSignificance <- stats::p.adjust(binomRF$significance, method = fdr.method)
  binomRF <- binomRF[order(binomRF$freq, decreasing=T),]
    

  
    if(keep.both){
      return(list(cor.binomRF=cor.binomRF,
                  binomRF =binomRF))
       } else{
         return(cor.binomRF=cor.binomRF)
       }
}
