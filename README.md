## amp - Automated Metabolomics Pipeline

This package can be used to analyze preprocessed LC-MS data. The starting point of the analyses conducted by this package is a peak table file, output from e.g. MSDIAL.

The package contains functions for visualizing and further preprocessing LC-MS data, as well as handy ways of conducting statistical analyses, mainly based on univariate linear models.



### Prerequisites and installation

PACKAGE requires R version 3.5.0 or greater. All dependencies are installed automatically.

##### Install PACKAGE

To install the package, run:

```
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("antonvsdata/amp")
```

*NOTE:* If this gives you a weird error, run the following line and try again (ref: https://github.com/r-lib/devtools/issues/1900)  
```
devtools::install_github("r-lib/remotes", ref = "e56a41e1d0cad55cbe7d60b274b99ab7b7a76b5c")
```

### Documentation

For instructions on how to use the package, run browseVignettes("amp")