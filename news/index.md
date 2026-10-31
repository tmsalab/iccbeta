# Changelog

## iccbeta 1.2.1

CRAN release: 2026-03-06

### Changes

- Added explicit dependencies on R (\>= 4.3.0), Rcpp (\>= 1.1.0), and
  RcppArmadillo (\>= 15.0.2-2)
- Removed CXX11 from `src/Makevars` and `src/Makevars.win` to avoid
  potential compilation issues with newer versions of Armadillo through
  RcppArmadillo.
- Switched README.Rmd to README.qmd to use Quarto for rendering.
- Fixed CITATION file to use [`c()`](https://rdrr.io/r/base/c.html)
  instead of [`personList()`](https://rdrr.io/r/utils/personList.html)
  and [`bibentry()`](https://rdrr.io/r/utils/bibentry.html) to avoid
  CRAN check notes.
- Switched from Travis-CI to GitHub Action workflows.

### Documentation

- Added a `pkgdown` website that deploys to
  <https://tmsalab.github.io/iccbeta/>.

## iccbeta 1.2.0

CRAN release: 2019-01-28

### Feature

- Added automatic calculation of icc on `lmer` model objects.

### Changes

- Added improved error messaging related to presence of missing values
  and dimension mis-match.
- Added ORCiDs to `DESCRIPTION`
- Increased version dependencies.
- Enabled use of OpenMP and C++11.
- Added `CITATION` information for the *R* package.

### Documentation

- Switched documentation to use Markdown.
- Improved contents of documentation.

### Deployment

- Enable the default TMSALab Travis-CI configuration.
- Enable code coverage checks.

## iccbeta 1.1.0

CRAN release: 2017-08-22

- Adds instructions on how to create `simICCdata2` to relevant helpdocs
- Switched over to using Rcpp to handle function registration instead of
  `src/init.c`.
- Added namespace requirement check for `lme4` and `RLRsim` to
  documentation.  
- Modified Travis-CI to use a later gcc

## iccbeta 1.0.1

CRAN release: 2017-04-26

- Added `src/init.c` to comply with R 3.4 requirements
- Added GitHub repository location information.
- Improved documentation flow.

## iccbeta 1.0.0

- Released initial version of iccbeta
