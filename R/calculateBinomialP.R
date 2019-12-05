#' calculate the probability, p, to conduct a binomial exact test
#'
#' \code{calculateBinomialP} returns a probability of randomly selecting a feature as the root node in a decision tree. This is a generic function that is called internally in \code{binomialRF} but that may also be called directly if needed. The arguments \code{...}
#' should be, L= Total number of features in X, and percent_features= what percent of L is subsampled in the \code{randomForest} call.
#'
#' @param L the total number of features in X. Should be a positive integer >1
#' @param percent_features what percentage of L do we subsample at each tree? Should be a proportion between (0,1)
#'
#' @return If L is an integeter returns a probability value for selecting predictor Xj randomly
#'
#' @examples
#' calculateBinomialP(110, .4)
#' calculateBinomialP(13200, .5)
#' @export


calculateBinomialP <- function(L, percent_features){

  if(!is.numeric(L) | !is.numeric(percent_features)){
    stop("Error: L or percent_features not numeric inputs")
  } else if( percent_features >1 | percent_features <0){
    stop("percent_features is outside the acceptable (0-1) range")
  } else if(L <2){
      stop('L must be a positive integer >1')
  }

  m = floor(L * percent_features)
  prod.vector = sapply(1:m, function(x) (L-x)/(L-(x-1)))
  (1-prod(prod.vector))*(1/m)
}

#' calculate the probability, p, to conduct a binomial exact test
#'
#' \code{calculateBinomialP} returns a probability of randomly selecting a feature as the root node in a decision tree. This is a generic function that is called internally in \code{binomialRF} but that may also be called directly if needed. The arguments \code{...}
#' should be, L= Total number of features in X, and percent_features= what percent of L is subsampled in the \code{randomForest} call.
#'
#' @param L the total number of features in X. Should be a positive integer >1
#' @param percent_features what percentage of L do we subsample at each tree? Should be a proportion between (0,1)
#' @param K interaction level
#'
#' @return If L is an integeter returns a probability value for selecting predictor Xj randomly
#'
#'
calculateBinomialP_Interaction <- function(L, percent_features, K=2){

  if(!is.numeric(L) | !is.numeric(percent_features)){
    stop("Error: L or percent_features not numeric inputs")
  } else if( percent_features >1 | percent_features <0){
    stop("percent_features is outside the acceptable (0-1) range")
  } else if(L <2){
    stop('L must be a positive integer >1')
  }

  m = floor(L * percent_features)

  prod.vector = lapply(1:K, function(K) sapply(1:m, function(x) ((L)-x)/((L)-(x-1))))

  p = prod(sapply(1:K, function(x)  (1-prod(prod.vector[[x]]))*1/m ))

  return(p)

}

