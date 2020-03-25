## notame - Workflow for non-targeted LC-MS metabolic profiling 

This package can be used to analyze preprocessed LC-MS data in non-targeted metabolomics (notame, see?). Notame was developed at the research group of nutritional metabolomics at University of Eastern Finland. We use notame as a way to bundle together all the preprocessing methods we use for our non-targeted LC-MS metabolomics data, so it mainly consists of methods found in other packages, and a bunch of visualizations we have found useful.

### What does notame do acutallly?

Before we go into the list of features, it is good for you to know hot the workflow in our lab works. The first step is to take raw data files created by the LC-MS instrument and create a peak table using a peak picking software (we use [MS-DIAL](http://prime.psc.riken.jp/Metabolomics_Software/MS-DIAL/)). After peak picking with the dedicated software, we use R for data preprocessing, quality control, statistical analysis and visualization. We then use the obtained results in identification of the actual metabolites. During the years, we ended up with various scripts that were hard to handle and update, so we decided to make this package to keep things under control. For more information about our workflow, read the associated protocol paper ["NoTaMe": Workflow for Non-Targeted LC-MS Metabolic Profiling](https://www.preprints.org/manuscript/202002.0019/v1)

Here is a list of the current main functionalities of the package:

- Reading data from Excel spreadsheets created with MS-DIAL
- Data is stored in a custom object that holds all the information about the features and samples along with the feature abundance matrix. This allows for a simple interface for all of the functions in the package, as there is no need to juggle with different matrices/data frames.
- Drift correction: correcting for systematic drift in the intensity of molecular features using cubic spline correction (see [Kirwan & Broadhurst et al.](https://doi.org/10.1007/s00216-013-6856-7))
- Identifying and flagging (or removing) low-quality molecular features using quality metrics defined by [Broadhurst et al.](https://doi.org/10.1007/s11306-018-1367-3)
- Imputing missing values, multiple strategies available. Random forest imputation recommended, see [Kokla et al.](https://doi.org/10.1186/s12859-019-3110-0)
- Batch effect correction: correcting for systematic variation between batches. Multiple strategies available.
- A novel method for clustering similar molecular features
- A bunch of statistical analyses, both feature-wise tests and multivariate models
- A rather nice set of visualizations for use in quality control, explorative analysis and interpretation of results from statistical tests


### Installation and getting started

PACKAGE requires R version 3.5.0 or greater.

#### Installation

notame functions depend on a ton of other R packages. The packages you need to install depend on what you're using notame for: some packages are only needed for specific visualizations, others for batch effect correction, and others for common preprocessing tasks. This is why it's recommended to only install the packages that are really needed to make notame work. To do this, run:

```
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("antonvsdata/notame", build_vignettes = TRUE)
```

After installing the package, you can install rest of the packages you need on the fly OR use a handy function called ```install_dependencies```, which lets you install packages for core preprocessing, batch correction, specific visualizations etc.


### Documentation

For instructions on how to use the package, run browseVignettes("notame"). The Introduction vignette should get you started pretty well!

##### Troubleshooting installation

If ```devtools::install_github``` gives you a weird error, run the following line and try again. Read more from [the original issue](https://github.com/r-lib/devtools/issues/1900)  
```
devtools::install_github("r-lib/remotes", ref = "e56a41e1d0cad55cbe7d60b274b99ab7b7a76b5c")
```

If you are having problems with new package versions (Something like: Error: (converted from warning) package 'ggplot2' was built under R version 3.6.1), try running the following line to prevent the error:

```
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = TRUE)
```

Read more [on the issue of remotes package](https://github.com/r-lib/remotes/issues/403) and [the environment variables tutorial](https://github.com/r-lib/remotes#environment-variables)


### Credits and license

The notame package is written by Anton Klåvus for his master's thesis in Bioinformatics at Aalto university (published under former name Anton Mattsson). Notame is inspired by analysis scripts written by Jussi Paananen, Oskari Timonen and Anton Klåvus (formerly Mattsson) at University of Eastern Finland. The algorithm for clustering molecular features originating from the same compound is based on MATLAB code written by David Broadhurst, Professor of Data Science & Biostatistics in the School of Science, and director of the Centre for Integrative Metabolomics & Computational Biology at the Edith Covan University.

If you find any bugs or other things to fix, please submit an issue on GitHub! All contributions to the package are always welcome!

notame is published under an MIT license (tl;dr: it's really permissive!)


