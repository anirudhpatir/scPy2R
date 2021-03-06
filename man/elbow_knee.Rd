% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{elbow_knee}
\alias{elbow_knee}
\title{elbow_knee}
\usage{
elbow_knee(x, addpca = 5)
}
\arguments{
\item{x}{float: the standard deviation of each principle component}

\item{addpca}{int: the PCs to add to those considered}
}
\value{
int: Number of PCs to consider
}
\description{
Identifies the principle components (PC) to use by using elbow knee method
}
