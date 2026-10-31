# A multilevel dataset from Hofmann, Griffin, and Gavin (2000).

A multilevel dataset from Hofmann, Griffin, and Gavin (2000).

## Usage

``` r
Hofmann
```

## Format

A data frame with 1,000 observations and 7 variables.

- `id`:

  a numeric vector of group ids.

- `helping`:

  a numeric vector of the helping outcome variable construct.

- `mood`:

  a level 1 mood predictor.

- `mood_grp_mn`:

  a level 2 variable of the group mean of mood.

- `cohesion`:

  a level 2 covariate measuring cohesion.

- `mood_grp_cent`:

  group-mean centered mood predictor.

- `mood_grd_cent`:

  grand-mean centered mood predictor.

## Source

Hofmann, D.A., Griffin, M.A., & Gavin, M.B. (2000). The application of
hierarchical linear modeling to management research. In K.J. Klein, &
S.W.J. Kozlowski (Eds.), Multilevel theory, research, and methods in
organizations: Foundations, extensions, and new directions (pp.
467-511). Hoboken, NJ: Jossey-Bass.

## References

Aguinis, H., & Culpepper, S.A. (2015). An expanded decision making
procedure for examining cross-level interaction effects with multilevel
modeling. *Organizational Research Methods*. Available at:
[http://hermanaguinis.com/pubs.html](http://hermanaguinis.com/pubs.md)

## See also

[`lmer`](https://rdrr.io/pkg/lme4/man/lmer.html),
[`model.matrix`](https://rdrr.io/r/stats/model.matrix.html),
[`VarCorr`](https://rdrr.io/pkg/nlme/man/VarCorr.html),
[`LRTSim`](https://rdrr.io/pkg/RLRsim/man/LRTSim.html),
[`simICCdata`](https://tmsalab.github.io/iccbeta/reference/simICCdata.md)

## Examples

``` r
if (FALSE) { # \dontrun{

if(requireNamespace("lme4") && requireNamespace("RLRsim")){ 
data(Hofmann)
library("lme4")

# Random-Intercepts Model
lmmHofmann0 = lmer(helping ~ (1|id), data = Hofmann)
vy_Hofmann = var(Hofmann[,'helping'])

# Computing icca
VarCorr(lmmHofmann0)$id[1,1]/vy_Hofmann

# Estimating Group-Mean Centered Random Slopes Model, no level 2 variables
lmmHofmann1  <- lmer(helping ~ mood_grp_cent + (mood_grp_cent |id),
                     data = Hofmann, REML = FALSE)
X_Hofmann = model.matrix(lmmHofmann1)
P = ncol(X_Hofmann)
T1_Hofmann  = VarCorr(lmmHofmann1)$id[1:P,1:P]

# Computing iccb
icc_beta(X_Hofmann, Hofmann[,'id'], T1_Hofmann, vy_Hofmann)$rho_beta

# Performing LR test
# Need to install 'RLRsim' package
library("RLRsim")
lmmHofmann1a  <- lmer(helping ~ mood_grp_cent + (1 | id),
                      data = Hofmann, REML = FALSE)
obs.LRT <- 2*(logLik(lmmHofmann1) - logLik(lmmHofmann1a))[1]
X <- getME(lmmHofmann1,"X")
Z <- t(as.matrix(getME(lmmHofmann1,"Zt")))
sim.LRT <- LRTSim(X, Z, 0, diag(ncol(Z)))
(pval <- mean(sim.LRT > obs.LRT))
} else {
 stop("Please install packages `RLRsim` and `lme4` to run the above example.")
}
} # }
```
