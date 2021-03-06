% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/funs_tools.R
\name{get_genes_desc}
\alias{get_genes_desc}
\title{Find features which separate the patient classes}
\usage{
get_genes_desc(mp, Aletter, Bletter)
}
\arguments{
\item{Aletter}{The letter signature of the labels of the first class samples}

\item{Bletter}{The letter signature of the labels of the second class samples}

\item{geno_p}{Matrix of profiles (genes x samples)}
}
\value{
Produce a dataframe s.t.  each row correspond to a gene situation
the columns are: rank_dist (higher and more it differs between classes)
diff_sd (higher and more its variance between classes is different)
}
\description{
Given a matrix of patient's features, this function ranks the features based on
how much their distribution values are different between classes
}
