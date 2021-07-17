#' random forest feature selection based on binomial exact test
#'
#' \code{k_binomialRF} is the R implementation of the interaction feature selection algorithm by (Zaim 2019). \code{k_binomialRF} extends the \code{binomialRF} algorithm by searching for k-way interactions.
#'
#' @param X design matrix
#' @param y class label
#' @param fdr.threshold fdr.threshold for determining which set of features are significant
#' @param fdr.method how should we adjust for multiple comparisons (i.e., \code{p.adjust.methods} =c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"))
#' @param ntrees how many trees should be used to grow the \code{randomForest}? (Defaults to 5000)
#' @param percent_features what percentage of L do we subsample at each tree? Should be a proportion between (0,1)
#' @param K for multi-way interactions, how deep should the interactions be?
#' @param cbinom_dist user-supplied correlated binomial distribution
#' @param sampsize user-supplied sample size for random forest
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
#' y = rbinom(100,1,pr)
#'
#' ###############################
#' ### Run interaction model
#' ###############################
#' 
#' require(correlbinom)
#' 
#' rho = 0.33
#' ntrees = 500
#' cbinom = correlbinom(rho, successprob =  calculateBinomialP_Interaction(10, .5,2), 
#'                                trials = ntrees, precision = 1024, model = 'kuk')
#'
#' k.binom.rf <-k_binomialRF(X,y, fdr.threshold = .05,fdr.method = 'BY',
#'                       ntrees = ntrees,percent_features = .5,
#'                       cbinom_dist=cbinom,
#'                       sampsize=round(nrow(X)*rho))
#'
#'
#' 
#'
#' @export

k_binomialRF <- function(X,y , fdr.threshold=0.05, fdr.method='BY', ntrees=2000, percent_features=0.3, K=2, cbinom_dist=NULL,
                         sampsize = nrow(X) *.4){

  if(!is.numeric(ntrees)  | !is.numeric(percent_features)| !is.numeric(fdr.threshold)){
    stop("Error: threshold, ntrees, and percent_features should be numeric inputs")
  } else if( percent_features >1 | percent_features <0){
    stop("percent_features is outside the acceptable (0-1) range")
  } else if(ntrees <2){
    stop('L must be a positive integer >1')
  } else if(!fdr.method %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")){
    stop('Please select acceptable fdr method from ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")')
  } else if(fdr.threshold >1 | fdr.threshold <0){
    stop("fdr.threshold is outside the acceptable (0-1) range")
  } else if( !is.numeric(K) | K< 2){
    stop('K should be an integeer denoting the level of interaction. Make K > 2')
  } else if(is.null(cbinom_dist)){
    stop('user_cbinom_dist is NULL, please calculate the correlated distribution using the correlbinom package')
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
  rf.object <- randomForest::randomForest(X,y, ntree = ntrees, mtry=m, keep.forest = TRUE, keep.inbag = TRUE)

  p = calculateBinomialP_Interaction(L,percent_features, K )

  climbToRoot <- function(final.node, Tree, K) {

    Tree = as.data.frame(Tree)
    path <- c(final.node)

    while(final.node != 1){
      if(final.node %% 2 ==0){
        new.parent <- which(Tree[,'left daughter'] == final.node)
        path <- c(path, new.parent)
        final.node <- new.parent

      } else if(final.node %% 2 != 0){
        new.parent <- which(Tree[,'right daughter'] == final.node)
        path <- c(path, new.parent)
        final.node <- new.parent
      }

    }

    if(length(path)>K){
      new.path <- rev(path)[-1]

      return(new.path)

    } else {
      return(path)
    }
  }

  findKInteractions <- function(Tree, K){
    Tree= as.data.frame(Tree)

    final.nodes <- (2^(K-1)):(2^K-1)

    interactions <- lapply(final.nodes, function(x) nms[sort(Tree[climbToRoot(final.node=x, Tree,K),'split var'])])
    pruned.interactions <- lapply(interactions, function(x) if(length(x)==K){return(x)})


    return(pruned.interactions)
  }

  ## obtain root-node splitting variables
  interaction.list <- lapply(1:ntrees, function(zzz) findKInteractions(randomForest::getTree(rf.object, k = zzz), K))

  interaction.list <- rlist::list.flatten(interaction.list)
  interaction.list <- do.call(rbind, interaction.list)

  interaction.list <- data.table::data.table(interaction.list)

  if(nrow(interaction.list) >0 ){
    interaction.list$Interaction <- do.call(paste, c(interaction.list, sep=" | "))
    
    sum.list<- data.frame(table(interaction.list$Interaction))
    colnames(sum.list) <- c('Interaction','Frequency')
    
    ### calculate correlation adjusted binomial
    pmf <- cbinom_dist/ sum(cbinom_dist)
    cmf <- cumsum(pmf)
    sum.list$significance <- 1-cmf[sum.list$Frequency]
    
    sum.list$adjSignificance <- stats::p.adjust(sum.list$significance, method = fdr.method)
    sum.list$Interaction <- as.character(sum.list$Interaction)
    
    sum.list <- sum.list[order(sum.list$Frequency, decreasing=T),]
    
    return(FeatureSelection=sum.list)
  } else {
    print('found no interactions')
  }


}

