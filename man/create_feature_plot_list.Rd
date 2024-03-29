% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/effect_plots.R
\name{create_feature_plot_list}
\alias{create_feature_plot_list}
\title{Generate a list of plots}
\usage{
create_feature_plot_list(object, plot_fun)
}
\arguments{
\item{object}{a MetaboSet object, should contain only features to be plotted}

\item{plot_fun}{function, a notame plot function}
}
\value{
a list of ggplot objects
}
\description{
Helper function for generating a list of feature-wise plots given a plot function
}
