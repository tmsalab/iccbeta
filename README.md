
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Build
Status](https://travis-ci.org/tmsalab/iccbeta.svg)](https://travis-ci.org/tmsalab/iccbeta)
[![Package-License](http://img.shields.io/badge/license-GPL%20\(%3E=2\)-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html)
[![CRAN Version
Badge](http://www.r-pkg.org/badges/version/iccbeta)](https://cran.r-project.org/package=iccbeta)
[![CRAN
Status](https://cranchecks.info/badges/worst/iccbeta)](https://cran.r-project.org/web/checks/check_results_iccbeta.html)
[![RStudio CRAN Mirror’s Monthly
Downloads](http://cranlogs.r-pkg.org/badges/iccbeta?color=brightgreen)](http://www.r-pkg.org/pkg/iccbeta)
[![RStudio CRAN Mirror’s Total
Downloads](http://cranlogs.r-pkg.org/badges/grand-total/iccbeta?color=brightgreen)](http://www.r-pkg.org/pkg/iccbeta)
[![Coverage
status](https://codecov.io/gh/tmsalab/iccbeta/branch/master/graph/badge.svg)](https://codecov.io/github/tmsalab/iccbeta?branch=master)

# `iccbeta` R package

A function and vignettes for computing an intraclass correlation
described in Aguinis & Culpepper (in press). iccbeta quantifies the
share of variance in a dependent variable that is attributed to group
heterogeneity in slopes.

## Installation

You can install `iccbeta` from CRAN using:

``` r
install.packages("iccbeta")
```

Or, you can be on the cutting-edge development version on GitHub using:

``` r
if(!requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("tmsalab/iccbeta")
```

## Usage

To use the `iccbeta` package, load it into *R* using:

``` r
library("iccbeta")
```

From there, the `icc_beta()` function can be called using either
`lmer()` model object or individual components to compute the intraclass
correlation:

``` r
# Automatically calculate icc from model
results_model = icc_beta(<lmer-model>)

# Calculate icc from individual terms.
results_manual = icc_beta(X, l2id, T, vy)
```

## Authors

Steven Andrew Culpepper and Herman Aguinis

## Citing the `iccbeta` package

To ensure future development of the package, please cite `iccbeta`
package if used during an analysis or simulation studies. Citation
information for the package may be acquired by using in *R*:

``` r
citation("iccbeta")
```

## License

GPL (\>= 2)
