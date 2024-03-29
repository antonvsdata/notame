% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/batch_correction.R
\name{normalize_batches}
\alias{normalize_batches}
\title{Normalize batches}
\usage{
normalize_batches(object, batch, group, ref_label, population = "all", ...)
}
\arguments{
\item{object}{a MetaboSet object}

\item{batch, group}{character, column names of pData with batch labels and group labels}

\item{ref_label}{the label of the reference group i.e. the group that is constant through batches}

\item{population}{Identifier of population samples in group column
(all (default) or any type of samples present in group)}

\item{...}{additional parameters passed to batchCorr::normalizeBatches}
}
\value{
list, the object with normalized features and information on which
features were corrected by ref samples in each batch.
}
\description{
Normalize bathces by either reference samples of population median.
Uses normalizeBatches function from the batchCorr package
}
\examples{
# Batch correction
\dontrun{
batch_corrected <- normalize_batches(merged_sample, batch = "Batch", group = "QC", ref_label = "QC")
# Evaluate batch correction
pca_bhattacharyya_dist(merged_sample, batch = "Batch")
pca_bhattacharyya_dist(batch_corrected, batch = "Batch")
}
}
