# iccbeta: Intraclass Correlation [![A hexagonal logo for iccbeta: a bundle of eight group regression lines meets at one dot and fans out at different slopes into an orange distribution of slopes labelled β₁ⱼ, above a solid bar with the word iccbeta cut out of it.](reference/figures/iccbeta-logo-animated-light.svg)](https://tmsalab.github.io/iccbeta/)

A function and vignettes for computing an intraclass correlation
described in Aguinis & Culpepper (in press). `iccbeta` quantifies the
share of variance in a dependent variable that is attributed to group
heterogeneity in slopes.

## Installation

You can install `iccbeta` from CRAN using:

``` r

install.packages("iccbeta")
```

Or, you can be on the cutting-edge development version on GitHub using:

``` r

if(!requireNamespace("remotes")) install.packages("remotes")
remotes::install_github("tmsalab/iccbeta")
```

## Usage

To use the `iccbeta` package, load it into *R* using:

``` r

library("iccbeta")
```

From there, calling the
[`icc_beta()`](https://tmsalab.github.io/iccbeta/reference/icc_beta.md)
function with either a `lmer()` model object or the desired individual
components will compute the intraclass correlation:

``` r
# Automatically calculate icc from model
results_model = icc_beta(<lmer-model>)

# Calculate icc from individual terms.
results_component = icc_beta(X, l2id, T, vy)
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
