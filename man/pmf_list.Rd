\name{pmf_list}
\docType{data}
\alias{pmf_list}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A prebuilt distribution for correlated binary data
%%  ~~function to do ... ~~
}
\description{
This data contains probability mass functions (pmf's) for correlated binary data for various parameters. The sum of correlated exchangeable binary data is a generalization of the binomial distribution that deals with correlated trials. The correlation in decision trees occurs as the subsampling and bootstrapping step in random forests touch the same data, creating a co-dependency. This data contains some pre-calculated distributions for random forests with 500, 1000, and 2000 trees with 10, 100, and 1000 features. For more distributions, they can be calculated via the correlbinom R package.
}
\usage{
pmf_list
}
\format{A list of lists}

\references{
Witt, Gary. "A Simple Distribution for the Sum of Correlated, Exchangeable Binary Data." Communications in Statistics-Theory and Methods 43, no. 20 (2014): 4265-4280.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datasets}% use one of  RShowDoc("KEYWORDS")
