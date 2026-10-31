# Simulated data example from Aguinis and Culpepper (2015).

A simulated data example from Aguinis and Culpepper (2015) to
demonstrate the
[`icc_beta`](https://tmsalab.github.io/iccbeta/reference/icc_beta.md)
function for computing the proportion of variance in the outcome
variable that is attributed to heterogeneity in slopes due to
higher-order processes/units.

## Usage

``` r
simICCdata
```

## Format

A data frame with 900 observations (i.e., 30 observations nested within
30 groups) on the following 6 variables.

- `l1id`:

  A within group ID variable.

- `l2id`:

  A group ID variable.

- `one`:

  A column of 1's for the intercept.

- `X1`:

  A simulated level 1 predictor.

- `X2`:

  A simulated level 1 predictor.

- `Y`:

  A simulated outcome variable.

## Source

Aguinis, H., & Culpepper, S.A. (2015). An expanded decision making
procedure for examining cross-level interaction effects with multilevel
modeling. *Organizational Research Methods*. Available at:
[http://www.hermanaguinis.com/pubs.html](http://www.hermanaguinis.com/pubs.md)

## Details

See Aguinis and Culpepper (2015) for the model used to simulate the
dataset.

## See also

[`lmer`](https://rdrr.io/pkg/lme4/man/lmer.html),
[`model.matrix`](https://rdrr.io/r/stats/model.matrix.html),
[`VarCorr`](https://rdrr.io/pkg/nlme/man/VarCorr.html),
[`LRTSim`](https://rdrr.io/pkg/RLRsim/man/LRTSim.html),
[`Hofmann`](https://tmsalab.github.io/iccbeta/reference/Hofmann.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(simICCdata)
if(requireNamespace("lme4")){ 
library("lme4")

# computing icca
vy <- var(simICCdata$Y)
lmm0 <- lmer(Y ~ (1|l2id), data = simICCdata, REML = FALSE)
VarCorr(lmm0)$l2id[1,1]/vy

# Create simICCdata2
grp_means = aggregate(simICCdata[c('X1','X2')], simICCdata['l2id'],mean)
colnames(grp_means)[2:3] = c('m_X1','m_X2')
simICCdata2 = merge(simICCdata, grp_means, by='l2id')

# Estimating random slopes model
lmm1  <- lmer(Y ~ I(X1-m_X1) + I(X2-m_X2) + (I(X1-m_X1) + I(X2-m_X2) | l2id),
              data = simICCdata2, REML = FALSE)
X <- model.matrix(lmm1)
p <- ncol(X)
T1 <- VarCorr(lmm1) $l2id[1:p,1:p]
# computing iccb
# Notice '+1' because icc_beta assumes l2ids are from 1 to 30.
icc_beta(X, simICCdata2$l2id+1, T1, vy)$rho_beta
} else {
 stop("Please install `lme4` to run the above example.")
}
} # }
```
