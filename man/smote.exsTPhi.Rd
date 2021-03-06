% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/smoteTPhi.R
\name{smote.exsTPhi}
\alias{smote.exsTPhi}
\title{Generate synthetic cases for smoteTPhi}
\usage{
smote.exsTPhi(data, tgt, N, k, dist, p, pc)
}
\arguments{
\item{data}{the rare cases (the minority "class" cases)}

\item{tgt}{the column nr of the target variable}

\item{N}{the percentage of over-sampling to carry out;}

\item{k}{the number of nearest neighours}

\item{dist}{the distance function used for the neighours computation}

\item{p}{an integer used when a "p-norm" distance is selected}

\item{pc}{relevance function}
}
\value{
SMOTEd Examples
}
\description{
The result of the function is a (N-1)*nrow(data) set of synthetically generated examples with rare values on the target
}
