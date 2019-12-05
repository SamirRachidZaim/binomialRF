#' a binomialRF_ensemble visualization and analysis for feature selection
#'
#' \code{binomialRF_ensemble} experimental section on plotting feature selection
#'
#' @param candidateModels takes as input a list of candidate models by specifying which predictors you want to consider for each model as a string of vector names
#' @param X the design matrix X of predictors.
#' @param y the outcome y as a factor variable.
#' @param ntrees how many trees should be used to grow the \code{randomForest}? (Defaults to 5000)
#' @param PLOT should a graphical summary be released (boolean)
#' @return a ggplot2 object comparison of different models

#' @examples
#' #' set.seed(324)
#' require(randomForest)
#' require(ggplot2)
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
#' ### Run ensemble
#' ###############################
#'
#' candidateModels <- list(
#'       m1=c('X1','X4','X6','X9'),
#'       m2=c(paste('X',c(1,2,3,4,6,9),sep='')),
#'       m3=c(paste("X", 3:10,sep='')),
#'       m4= c(paste("X", 5:9,sep='')),
#'       m5= c(paste("X", c(4,5,8,10),sep='')),
#'       m6= c(paste("X",1:10,sep=''))
#'       )
#'

.binomialRF_ensemble <- function(candidateModels, X, y, ntrees=2000, PLOT=F){

  if( !is.data.frame(X) ){
    X = data.frame(X)
  }
  if(!is.factor(y)){
    y = factor(y)
  }
  dim.X = ncol(X)

  f <- function(i, candidateModels, X,y){
    
    tmp_cbinom_dist <- .fast_correlbinom(rho = .63,trials = ntrees, successprob = 1/length(candidateModels[[i]]), precision = 1024)
    
    d= binomialRF(X[, candidateModels[[i]]], factor(y), ntrees = ntrees, user_cbinom_dist =tmp_cbinom_dist)
    d$model = names(candidateModels)[i]
    return(d)
  }

  candidateList = lapply(1:length(candidateModels), function(i) f(i, candidateModels, X,y) )
  numModels <- length(candidateList)

  new.err.mat <- do.call(rbind, candidateList)

  models <- unlist(sapply(1:length(candidateList), function(X) rep(paste('m',X,sep=''), nrow(candidateList[[X]])))  )

  new.err.mat$model <- models
  new.err.mat$Norm.Weight = new.err.mat$weight / sum(unique(new.err.mat$weight))
  new.err.mat <- new.err.mat[order(new.err.mat$Norm.Weight,decreasing=TRUE), ]



  # err.mat$Significant = as.numeric(err.mat$Significant)
  # err.mat$Norm.Weight = err.mat$weight / sum(unique(err.mat$weight))
  # err.mat$Significant = as.logical(err.mat$Significant)
  # err.mat= err.mat[err.mat$Variable!='',]
  # w= round(unique(err.mat$Norm.Weight),3)
  # err.mat$Variable <- factor(err.mat$Variable, levels = names(data.frame(X))[c(1:10)])


  # err.mat_weighted <- err.mat[order(err.mat$Norm.Weight,decreasing=TRUE), ]
  # err.mat_weighted$model <- factor(err.mat_weighted$model, levels = unique(err.mat_weighted$model))
  # err.mat_weighted$Variable <- factor(err.mat_weighted$Variable, levels = names(data.frame(X)))

  # plt = ggplot(data=err.mat_weighted, aes(x=model, y=Variable,  fill=Significant, width=Norm.Weight)) +
  #   geom_tile(inherit.aes = F , aes(x=model,y=Variable,fill=Significant, width=Norm.Weight), position = position_identity())+
  #   theme_minimal()+ labs(title='Feature Selection by Model',
  #                         x='Likeliest Candidate Model (OOB Error)',
  #                         y='Feature',
  #                         fill="Significant")+
  #   theme(
  #     plot.title = element_text(color="red", size=20, face="bold.italic"),
  #     axis.title.x = element_text(color="blue", size=14, face="bold"),
  #     axis.title.y = element_text(color="#993333", size=14, face="bold")
  #   )

  if(PLOT){
    
    plt = ggplot(data=new.err.mat, aes(x=new.err.mat$model, y=new.err.mat$Variable,  fill=new.err.mat$Significant)) +
      geom_tile(inherit.aes = FALSE , aes(x=new.err.mat$model,y=new.err.mat$Variable,fill=new.err.mat$Significant), position = position_identity())+
      theme_minimal()+ labs(title='Feature Selection by Model',
                                              x='Likeliest Candidate Model (OOB Error)',
                                              y='Feature',
                                              fill="Significant")+
      theme(
        plot.title = element_text(color="red", size=20, face="bold.italic"),
        axis.title.x = element_text(color="blue", size=14, face="bold"),
        axis.title.y = element_text(color="#993333", size=14, face="bold")
      )
    
    print(plt)
    
  }

  err.mat <- new.err.mat

  err.mat=data.table::data.table(err.mat)
  # data.table::setkey(err.mat, model)


  err.mat$Significant=as.numeric(err.mat$Significant)
  err.mat$Significant.Weight = err.mat$Significant
  w= round(unique(err.mat$Norm.Weight),3)


  err.mat = data.table::dcast(err.mat, formula = Variable ~ model)

  prop.selected = function(row.vals){
    n1 = sum(row.vals>0,na.rm = TRUE)
    n2 = sum(!is.na(row.vals))
    n1/n2
  }

  err.mat$Prop.Selected = sapply(1:nrow(err.mat), function(x) prop.selected(err.mat[x, 2:ncol(err.mat)]))
  err.mat <- err.mat[, c('Variable','Prop.Selected', paste('m',1:numModels, sep='')), with=FALSE]
  err.mat[, 2:length(err.mat)] <- round(err.mat[, 2:ncol(err.mat)],3)
  err.mat = err.mat[order(err.mat$Variable),]

  err.mat <- err.mat[order(err.mat$Prop.Selected, decreasing = TRUE),]
  err.mat[is.na(err.mat)] <- ''


  # ### organize table
  # a = data.frame(t(e(c('OOB_Weight','',w)))); colnames(a) <- c('Variable','Prop.Selected', paste('m','',1:numModels, sep=''))
  # b = data.frame(t(e(c(rep('',(numModels+2)))))); colnames(b) <- c('Variable', 'Prop.Selected', paste('m',1:numModels, sep=''))
  # err.mat =rbind(err.mat,b,a )
  # err.mat = err.mat[, c(1:2, order(a,decreasing = T)), with=F]
  # err.mat = err.mat[, unique(names(err.mat)), with=F]


  return(err.mat)
}
