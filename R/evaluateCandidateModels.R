#' a binomialRf visualization for feature selection
#'
#' \code{evaluateCandidateModels} experimental section on plotting feature selection
#'
#' @param candidateModels takes as input a list of candidate models by specifying which predictors you want to consider for each model as a string of vector names
#' @param X the design matrix X of predictors.
#' @param y the outcome y as a factor variable.
#'
#' @return a ggplot2 object comparison of different models

#' @examples
#' candidateModels <- list(
#'       m1=c('X1','X4','X6','X9'),
#'       m2=c(paste('X',c(1,2,3,4,6,9),sep='')),
#'       m3=c(paste("X", 3:10,sep='')),
#'       m4= c(paste("X", 5:9,sep='')),
#'       m5= c(paste("X", c(4,5,8,10),sep='')),
#'       m6= c(paste("X",1:10,sep=''))
#'       )
#'
#' evaluateCandidateModels(candidateModels, X,y)

evaluateCandidateModels <- function(candidateModels, X, y, ntrees=20000, percent_features=.2){
  library(dplyr)
  library(ggplot2)
  library(data.table)

  X = data.frame(X)
  y = factor(y)
  dim.X = ncol(X)

  f <- function(i, candidateModels, X,y){
    d= binomialRF(X[, candidateModels[[i]]], factor(y), percent_features = percent_features, ntrees = ntrees)
    d$model = names(candidateModels)[i]
    return(d)
  }

  candidateList = lapply(1:length(candidateModels), function(i) f(i, candidateModels, X,y) )
  numModels <- length(candidateList)

  err.mat= data.frame(Variable= character(dim.X*numModels),
                      Significant=logical(dim.X*numModels),
                      weight=numeric(dim.X*numModels),
                      model=t(do.call(cbind, lapply(1:numModels, function(x) t(rep(paste('m',x,sep=''),dim.X))))),
                      stringsAsFactors = F)


  i=0
  for(model in candidateList){
    a = ((i*10)+1)
    err.mat[a:(a+nrow(model)-1),1:3] <- model[, c('Variable','Significant','weight')]
    i=i+1
  }

  err.mat$Significant = as.numeric(err.mat$Significant)
  err.mat$Norm.Weight = err.mat$weight / sum(unique(err.mat$weight))
  err.mat$Significant = as.logical(err.mat$Significant)
  err.mat= err.mat[err.mat$Variable!='',]
  w= round(unique(err.mat$Norm.Weight),3)
  err.mat$Variable <- factor(err.mat$Variable, levels = names(data.frame(X))[c(1:10)])


  err.mat_weighted <- err.mat[order(err.mat$Norm.Weight,decreasing=TRUE), ]
  err.mat_weighted$model <- factor(err.mat_weighted$model, levels = unique(err.mat_weighted$model))
  err.mat_weighted$Variable <- factor(err.mat_weighted$Variable, levels = names(data.frame(X)))

  plt = ggplot(data=err.mat_weighted, aes(x=model, y=Variable,  fill=Significant, width=Norm.Weight)) +
    geom_tile(inherit.aes = F , aes(x=model,y=Variable,fill=Significant, width=Norm.Weight), position = position_identity())+
    theme_minimal()+ labs(title='Feature Selection by Model',
                          x='Likeliest Candidate Model (OOB Error)',
                          y='Feature',
                          fill="Significant")+
    theme(
      plot.title = element_text(color="red", size=20, face="bold.italic"),
      axis.title.x = element_text(color="blue", size=14, face="bold"),
      axis.title.y = element_text(color="#993333", size=14, face="bold")
    )


  err.mat=data.table::data.table(err.mat)
  setkey(err.mat, model)

  print(plt)

  err.mat$Significant=as.numeric(err.mat$Significant)
  err.mat$Significant.Weight = err.mat$Significant* err.mat$Norm.Weight


  err.mat = dcast(err.mat, formula = Variable ~ model)


  prop.selected = function(row.vals){
    n1 = sum(row.vals>0,na.rm = T)
    n2 = sum(!is.na(row.vals))
    n1/n2
  }

  err.mat$Prop.Selected = sapply(1:dim.X, function(x) prop.selected(err.mat[x, 2:(numModels-1)]))
  err.mat <- err.mat[, c('Variable','Prop.Selected', paste('m',1:numModels, sep='')), with=F]
  err.mat[, 2:length(err.mat)] <- round(err.mat[, 2:ncol(err.mat)],3)
  err.mat = err.mat[order(err.mat$Variable),]


  ### organize table
  a = data.frame(t(data.table(c('OOB_Weight','',w)))); colnames(a) <- c('Variable','Prop.Selected', paste('m','',1:numModels, sep=''))
  b = data.frame(t(data.table(c(rep('',(numModels+2)))))); colnames(b) <- c('Variable', 'Prop.Selected', paste('m',1:numModels, sep=''))
  err.mat =rbind(err.mat,b,a )
  err.mat = err.mat[, c(1:2, order(a,decreasing = T)), with=F]
  err.mat = err.mat[, unique(names(err.mat)), with=F]

  return(err.mat)
}
