% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/funs_tools.R
\name{get_B_prob}
\alias{get_B_prob}
\title{Determines testing patient probability of belonging to the second patient class based on BFC algorithm}
\usage{
get_B_prob(tepat, m_sim, m_sim_cl, vars_list, best_direction)
}
\arguments{
\item{tepat}{String character of the name of the testing profile}

\item{m_sim}{Matrix: Original PSN with all the profiles both training and testing}

\item{m_sim_cl}{Matrix: Clustering PSN of only training profiles}

\item{vars_list}{Simpati object: List of variables}

\item{best_direction}{Power info result}
}
\value{
A dataframe of information for the profile class prediction, precisely:
Btop: the top percentage in which the profile is included
BBe: the median of the profile's BB edges
DIFFe: the difference beteen the medians of profile's BB and AB edges
assumption: the assumption made with the function: BSign: the profile is fake of second
and signature class
}
\description{
Assumes that the testing patient belongs to the second class, it iterates the Best Friend Connector
algorithm to find how much it fits in the representative members of the class
}
