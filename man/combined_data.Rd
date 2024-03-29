% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_constructor.R
\name{combined_data}
\alias{combined_data}
\title{Retrieve both sample information and features}
\usage{
combined_data(object)
}
\arguments{
\item{object}{a MetaboSet object}
}
\description{
Returns a data frame with sample information plus all features as columns.
The data frame thus has one row per sample.
}
\examples{
combined_data(example_set)

}
