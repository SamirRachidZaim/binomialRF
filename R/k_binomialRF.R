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
#' @param keep.rf should we keep the randomForest object?
#' @param K for multi-way interactions, how deep should the interactions be?
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


.k_binomialRF <- function(X,y , fdr.threshold=0.05, fdr.method='BY', ntrees=2000, percent_features=0.3, keep.rf =FALSE, K=2){

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
  } else if( !is.numeric(K) | K< 2){
    stop('K should be an integeer denoting the level of interaction. Make K > 2')
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
  rf.object <- randomForest(X,y, ntree = ntrees, mtry=m, keep.forest = TRUE, keep.inbag = TRUE)

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


  #interaction.list <- interaction.list[interaction.list$V1!='X0',]

  interaction.list$Interaction <- do.call(paste, c(interaction.list, sep=" | "))

  sum.list <- interaction.list[, .N , by =interaction.list$Interaction]
  sum.list$Pvalue <- as.numeric(sapply(sum.list$N, function(x) binom.test(x, n= ntrees, p, alternative='greater')$p.value))
  sum.list$Adj.Pval <- stats::p.adjust(sum.list$Pvalue, method = fdr.method)
  sum.list$Significant <- sum.list$Adj.Pval < fdr.threshold
  sum.list <- sum.list[order(sum.list$Adj.Pval, decreasing = FALSE),]

  if(keep.rf){
    return(list(binomRF = sum.list,
                RF.object = rf.object,
                OutOfBagError = colMeans(rf.object$err.rate)[1]))
  } else{
   return(FeatureSelection=sum.list)
  }
}

