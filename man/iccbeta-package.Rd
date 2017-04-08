\name{iccbeta-package}
\alias{iccbeta-package}
\alias{iccbeta}
\docType{package}
\title{
Intraclass correlation used to assess variability of lower-order relationships across higher-order processes/units.
}
\description{
A function and vignettes for computing the intraclass correlation described in Aguinis & Culpepper (in press). iccbeta quantifies the share of variance in an outcome variable that is attributed to heterogeneity in slopes due to higher-order processes/units. 
}
\details{
\tabular{ll}{
Package: \tab iccbeta\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-11-21\cr
License: \tab GPL (>= 2)\cr
}
}
\author{
Steven Andrew Culpepper,
Herman Aguinis

Maintainer: Steven Andrew Culpepper <sculpepp@illinois.edu>
}
\references{
Aguinis, H., & Culpepper, S.A. (in press). An expanded decision making procedure for examining cross-level interaction effects with multilevel modeling. \emph{Organizational Research Methods}. Available at: \url{http://mypage.iu.edu/~haguinis/pubs.html}
}

\keyword{ package }

\examples{
\dontrun{
#Simulated Data Example
data(simICCdata)
  require(lme4)
  
  #computing icca
  vy = var(simICCdata$Y)
  lmm0 <- lmer(Y ~ (1|l2id),data=simICCdata,REML=F)
  VarCorr(lmm0)$l2id[1,1]/vy
    
  #Estimating random slopes model
  lmm1  <- lmer(Y~I(X1-m_X1)+I(X2-m_X2) +(I(X1-m_X1)+I(X2-m_X2)|l2id),data=simICCdata2,REML=F)
  X = model.matrix(lmm1)
  p=ncol(X)
  T1  = VarCorr(lmm1) $l2id[1:p,1:p]
  #computing iccb
  #Notice '+1' because icc_beta assumes l2ids are from 1 to 30.
  icc_beta(X,simICCdata2$l2id+1,T1,vy)$rho_beta
  
#Hofmann 2000 Example
data(Hofmann)
  require(lme4)

  #Random-Intercepts Model
  lmmHofmann0 = lmer(helping ~ (1|id),data=Hofmann)
  vy_Hofmann = var(Hofmann[,'helping'])
  #computing icca
  VarCorr(lmmHofmann0)$id[1,1]/vy_Hofmann

  #Estimating Group-Mean Centered Random Slopes Model, no level 2 variables
  lmmHofmann1  <- lmer(helping ~ mood_grp_cent + (mood_grp_cent |id),data=Hofmann,REML=F)
  X_Hofmann = model.matrix(lmmHofmann1)
  P = ncol(X_Hofmann)
  T1_Hofmann  = VarCorr(lmmHofmann1)$id[1:P,1:P]
  #computing iccb
  icc_beta(X_Hofmann,Hofmann[,'id'],T1_Hofmann,vy_Hofmann)$rho_beta
  
  #Performing LR test
  library('RLRsim')
  lmmHofmann1a  <- lmer(helping ~ mood_grp_cent + (1 |id),data=Hofmann,REML=F)
  obs.LRT <- 2*(logLik(lmmHofmann1)-logLik(lmmHofmann1a))[1]
  X <- getME(lmmHofmann1,"X")
  Z <- t(as.matrix(getME(lmmHofmann1,"Zt")))
  sim.LRT <- LRTSim(X, Z, 0, diag(ncol(Z)))
  (pval <- mean(sim.LRT > obs.LRT))
}
}
