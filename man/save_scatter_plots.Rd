% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/effect_plots.R
\name{save_scatter_plots}
\alias{save_scatter_plots}
\title{Save scatter plots of each feature against a set variable}
\usage{
save_scatter_plots(
  object,
  x = "Injection_order",
  save = TRUE,
  file_path = NULL,
  format = "emf",
  all_features = FALSE,
  color = NULL,
  color_scale = NA,
  shape = NULL,
  title = "Feature_ID",
  subtitle = NULL,
  shape_scale = getOption("notame.shape_scale"),
  text_base_size = 14,
  point_size = 2,
  title_line_length = 40,
  theme = theme_bw(base_size = text_base_size),
  ...
)
}
\arguments{
\item{object}{a MetaboSet object}

\item{x}{character, name of the column to be used as x-axis}

\item{file_path}{character, a file path for PDF or prefix added to the file paths for other formats}

\item{format}{character, format in which the plots should be saved}

\item{all_features}{logical, should all features be used? If FALSE
(the default), flagged features are removed before visualization.}

\item{color}{character, name of the column to be used for coloring}

\item{color_scale}{the color scale as returned by a ggplot function.
Set to NA to choose the appropriate scale based on the class of the coloring variable.}

\item{shape}{character, name of the column used for shape}

\item{title, subtitle}{column names from fData to use as plot title/filename and subtitle.
Set to NULL for no title/subtitle, this creates running numbered filenames}

\item{shape_scale}{the shape scale as returned by a ggplot function}

\item{text_base_size}{integer, base size for text in figures}

\item{point_size}{numeric, size of the points}

\item{title_line_length}{integer, maximum length of the title line in characters, passed to stringr::str_wrap}

\item{theme}{a ggplot theme to be added to the plot}

\item{...}{other arguments to graphic device functions, like width and height}
}
\description{
Draws a scatterplots with a feature on y-axis and another variable on x-axis.
A separate plot is drawn and saved for each feature.
}
\examples{
\dontrun{
# Against injection order, colored by group
save_scatter_plots(
  object = merged_sample[1:10],
  x = "Injection_order",
  color = "Group",
  file_path = "./scatter_plots.pdf",
  format = "pdf"
)
}
# Plot one feature
save_scatter_plots(merged_sample[5, ], save = FALSE, color = "Group")
}
\seealso{
\code{\link[notame]{save_plot}}
}
