% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{get_p_matrix}
\alias{get_p_matrix}
\title{calculate p matrix}
\usage{
get_p_matrix(branch_length, eig = phangorn::edQt(), rate = 1)
}
\arguments{
\item{branch_length}{branch length}

\item{eig}{eigen object}

\item{rate}{rate}
}
\value{
p matrix
}
\description{
calculates the p matrix
}
