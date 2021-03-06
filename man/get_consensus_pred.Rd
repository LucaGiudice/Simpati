% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/funs_tune.R
\name{get_consensus_pred}
\alias{get_consensus_pred}
\title{Determine the consensus class of a testing patient from a run of classification}
\usage{
get_consensus_pred(df)
}
\arguments{
\item{df}{Matrix which list the signature pathways with their corresponding predicted class per testing patient}
}
\value{
Produce a new matrix s.t. each gene has the rank of its original value in the sample
}
\description{
Given the classes assigned to a testing patient by all the signature PSN during the classification phase, it
determines the final class based on the frequency
}
