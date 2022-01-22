\name{fit_bdmulti}
\alias{fit_bdmulti}
\title{
Fitting a birth-death diversification model across multiple trees
}
\description{
Fits a single birth-death model with potentially time-varying rates and potentially missing extant species to a list of multiple phylogenies. The notation is similar to one used in fit_bd in RPANDA package.
}
\usage{
fit_bdmulti(phylolist, tot_timelist, f.lamb, f.mu, lamb_par, mu_par, flist, meth = "Nelder-Mead", cst.lamb = FALSE, cst.mu = FALSE, expo.lamb = FALSE, expo.mu = FALSE, fix.mu = FALSE, dt = 0, cond = "crown")
}
\arguments{
  \item{phylolist}{
a list of objects of type 'phylo'
}
  \item{tot_timelist}{
a list of ages of phylogenies (crown age, or stem age if known)
}
  \item{f.lamb}{
a function specifying the hypothesized functional form (e.g. constant, linear, exponential, etc.) of the variation of the speciation rate \eqn{\lambda} with time. Any functional form may be used.
This function has two arguments: the first argument is time; the second argument is a numeric vector of the parameters of the time-variation (to be estimated).
}
  \item{f.mu}{
a function specifying the hypothesized functional form (e.g. constant, linear, exponential, etc.) of the variation of the extinction rate \eqn{\mu} with time. Any functional form may be used.
This function has two arguments: the first argument is time; the second argument is a numeric vector of the parameters of the time-variation (to be estimated).
}
  \item{lamb_par}{
a numeric vector of initial values for the parameters of f.lamb to be estimated (these values are used by the optimization algorithm). The length of this vector is used to compute the total number of parameters in the model, so to fit a model with constant speciation rate (for example), lamb_par should be a vector of length 1. Otherwise aic values will be wrong.
}
  \item{mu_par}{
a numeric vector of initial values for the parameters of f.mu to be estimated (these values are used by the optimization algorithm). The length of this vector is used to compute the total number of parameters in the model, so to fit a model without extinction (for example), mu_par should be empty (vector of length 0). Otherwise aic values will be wrong.
}
  \item{flist}{
a list of fractions of extant species included in each phylogeny
}
  \item{meth}{
optimization to use to maximize the likelihood function, see \emph{optim} for more details.The optimization setup generally requires more attention that in fit_bd, due to more complex likelihood surfaces.
}
  \item{cst.lamb}{
logical: should be set to TRUE only if f.lamb is constant (i.e. does not depend on time) to use analytical instead of numerical computation in order to reduce computation time.
}
  \item{cst.mu}{
logical: should be set to TRUE only if f.mu is constant (i.e. does not depend on time) to use analytical instead of numerical computation in order to reduce computation time.
}
  \item{expo.lamb}{
logical: should be set to TRUE only if f.lamb is exponential to use analytical instead of numerical computation in order to reduce computation time.
}
  \item{expo.mu}{
logical: should be set to TRUE only if f.mu is exponential to use analytical instead of numerical computation in order to reduce computation time.
}
  \item{fix.mu}{
logical: if set to TRUE, the extinction rate \eqn{\mu} is fixed and will not be optimized.
}
  \item{dt}{
the default value is 0. In this case, integrals in the likelihood are computed using R "integrate" function, which can be quite slow. If a positive dt is given as argument, integrals are computed using a piece-wise contant approximation, and dt represents the length of the intervals on which functions are assumed to be constant.
}
  \item{cond}{
conditioning to use to fit the model:
\itemize{ \item FALSE: no conditioning (not recommended);
\item "stem": conditioning on the survival of the stem lineage (use when the stem age is known, in this case tot_time should be the stem age);
\item "crown" (default): conditioning on a speciation event at the crown age and survival of the 2 daugther lineages (use when the stem age is not known, in this case tot_time should be the crown age).
}
}
}

\value{
	a list with the following components
  \item{model}{the name of the fitted model}
  \item{LH}{the maximum log-likelihood value}
  \item{aicc}{the second order Akaike's Information Criterion}
  \item{lamb_par}{a numeric vector of estimated f.lamb parameters, in the same order as defined in f.lamb}
  \item{mu_par}{a numeric vector of estimated f.mu parameters, in the same order as defined in f.mu (if fix.mu is FALSE)}
}

\details{
The lengths of lamb_par and mu_par are used to compute the total number of parameters in the model, so to fit a model with constant speciation rate (for example), lamb_par should be a vector of length 1. Otherwise aic values will be wrong. In the f.lamb and f.mu functions, the first argument (time) runs from the present to the past. Hence, if the parameter controlling the variation of \eqn{\lambda} with time is estimated to be positive (for example), this means that the speciation rate decreases from past to present.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
J Smyčka
}
\seealso{
\code{\link{RPANDA::fit_bd}}}
\examples{
#The model fit takes a long time (around 30 seconds on standard desktop) due to complex likelihood surface to be explored. Be patient or use rougher optimization setup.

#generate a list of 3 phylogenies with lambda=0.1 running over 30 time units
tlist=sim_bdmulti(3,0.1,0,30)


#get total time for each tree
tot_timefn=function(phylo){max(node.age(phylo)$ages)}
tot_timelist=lapply(tlist,tot_timefn)

#list of sampling proportion, complete sampling in our case
flist=as.list(rep(1,3))

#starting points for ML search
lamb_par<-c(1)
mu_par<-c()

#functions of time dependence (constant and zero)
f.lamb <-function(t,y){y[1]}
f.mu<-function(t,y){0}

#fit
multi=fit_bdmulti(tlist, tot_timelist, f.lamb, f.mu, lamb_par, mu_par,flist,fix.mu=T)

multi
multi$lamb_par
multi$LH
}


