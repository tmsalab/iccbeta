## Test environments
* local OS X install, R 3.4.0
* ubuntu 12.04 (on travis-ci), R 3.4.0
* win-builder (devel and release)

## Feedback

The previously submitted version of this package had a CRAN maintainer request
preprint information regarding Aguinis & Culpepper (in press). This has been 
updated in the description and relevant documentation to point to the paper
published in 2015 on the authors respective website. Furthermore, we have
changed the sentence in the `DESCRIPTION` file from:

> iccbeta quantifies the share ...

To:

> This package quantifies the share ...

## R CMD check results

0 errors | 0 warnings | 1 note

We have one note related to spelling:

Possibly mis-spelled words in DESCRIPTION:
  Aguinis (11:18)
  Culpepper (11:28)
  Intraclass (3:25)
  iccbeta (11:50)
  intraclass (10:56)
  
The first two are the authors' last names. The last three are correct spellings,
with intraclass being the method and iccbeta being the package name.

## Reverse dependencies

There are no reverse dependencies for this packge. 
