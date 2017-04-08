\name{simICCdata}
\alias{simICCdata}
\docType{data}
\title{
Simulated data example from Aguinis and Culpepper (in press).
}
\description{
A simulated data example from Aguinis and Culpepper (in press) to demonstrate the icc_beta function for computing the proportion of variance in the outcome variable that is attributed to heterogeneity in slopes due to higher-order processes/units. 
}
\usage{data(simICCdata)}
\format{
  A data frame with 900 observations (i.e., 30 observations nested within 30 groups) on the following 6 variables.
  \describe{
    \item{\code{l1id}}{A within group ID variable.}
    \item{\code{l2id}}{A group ID variable.}
    \item{\code{one}}{A column of 1's for the intercept.}
    \item{\code{X1}}{A simulated level 1 predictor.}
    \item{\code{X2}}{A simulated level 1 predictor.}
    \item{\code{Y}}{A simulated outcome variable.}
  }
}
\details{
See Aguinis and Culpepper (in press) for the model used to simulate the dataset. 
}
\author{
Steven Andrew Culpepper,
Herman Aguinis

Maintainer: Steven Andrew Culpepper <sculpepp@illinois.edu>
}
\source{
Aguinis, H., & Culpepper, S.A. (in press). An expanded decision making procedure for examining cross-level interaction effects with multilevel modeling. \emph{Organizational Research Methods}. Available at: \url{http://mypage.iu.edu/~haguinis/pubs.html}
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\seealso{
\code{\link[lme4]{lmer}}, \code{\link{model.matrix}}, \code{\link[lme4]{VarCorr}}, \code{\link[RLRsim]{LRTSim}}, \code{\link{Hofmann}}
}
\examples{
\dontrun{
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
}
}
\keyword{datasets}
