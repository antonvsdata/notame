## amp - Automated Metabolomics Pipeline

This package can be used to analyze preprocessed LC-MS data. The starting point of the analyses conducted by this package is a peak table file, output from e.g. MS-DIAL.

The package contains functions for visualizing and further preprocessing LC-MS data, as well as handy ways of conducting simple statistical analyses.


### Installation and getting started

PACKAGE requires R version 3.5.0 or greater.

#### Installation

To install the package along with all the possible dependencies, run:

```
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("antonvsdata/amp", dependencies = c("Depends", "Imports", "Suggests"))
```
If you are having trouble installing some packages, you can install the bare minimum requirements by running the following code. Note that if you choose this type of installation, many of the functions in the package will not work as expected, as the pacakge has many dependencies.


```
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("antonvsdata/amp")
```

### Documentation

For instructions on how to use the package, run browseVignettes("amp"). The Introduction vignette should get you started pretty well!

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

