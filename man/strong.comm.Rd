% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/funs_top_clustering_bffs_v7.R
\name{strong.comm}
\alias{strong.comm}
\title{Best Friend Connector algorithm}
\usage{
strong.comm(m, thW = 0.6, thGroup = 0.6)
}
\arguments{
\item{m}{adjacency matrix of a patient similarity network}

\item{thW}{Default 0.6, double value between 0 and 1 determining the percentage of best friends to considers for a root (1BFS are the 60 percentage of the patients most similar to the root)}

\item{thGroup}{Default 0.6, double value between 0 and 1 determining the percentage of best friends to keep in the end as final subclass}
}
\value{
name of the nodes/columns which compose the best group with the most similar best friends and the root
}
\description{
Undefined
}
