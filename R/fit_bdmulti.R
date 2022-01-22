fit_bdmulti=function (phylolist, tot_timelist, f.lamb, f.mu, lamb_par, mu_par, flist, 
          meth = "Nelder-Mead", cst.lamb = FALSE, cst.mu = FALSE, expo.lamb = FALSE, 
          expo.mu = FALSE, fix.mu = FALSE, dt = 0, cond = "crown") 
{
  #the calculation of total n for aicc
  ntiplist=lapply(phylolist,Ntip)
  nobs=Reduce("+",ntiplist)
  
  if (fix.mu == FALSE) {
    init <- c(lamb_par, mu_par)
    p <- length(init)
    optimLH <- function(init) {
      lamb_par <- init[1:length(lamb_par)]
      mu_par <- init[(1 + length(lamb_par)):length(init)]
      f.lamb.par <- function(t) {
        abs(f.lamb(t, lamb_par))
      }
      f.mu.par <- function(t) {
        abs(f.mu(t, mu_par))
      }
      #the LH is now a sum of likelihoods for different trees, each tree can have a different total age (tot_time) 
      #and sampling proportion (f)
      lhsinpar=function(phylo, tot_time, f){likelihood_bd(phylo, tot_time, f.lamb.par, 
                                            f.mu.par, f, cst.lamb = cst.lamb, cst.mu = cst.mu, 
                                            expo.lamb = expo.lamb, expo.mu = expo.mu, dt = dt, 
                                            cond = cond)}
      likelihoodlist=mapply(lhsinpar, phylolist, tot_timelist, flist)
      LH=Reduce("+",likelihoodlist)
      
      return(-LH)
    }
    temp <- suppressWarnings(optim(init, optimLH, method = meth))
    lamb.par <- temp$par[1:length(lamb_par)]
    mu.par <- temp$par[(1 + length(lamb_par)):length(init)]
    f.lamb.par <- function(t) {
      f.lamb(t, lamb.par)
    }
    f.mu.par <- function(t) {
      f.mu(t, mu.par)
    }
    res <- list(model = "birth death", LH = -temp$value, 
                aicc = 2 * temp$value + 2 * p + (2 * p * (p + 1))/(nobs - 
                                                                     p - 1), lamb_par = lamb.par, mu_par = mu.par, 
                f.lamb = Vectorize(f.lamb.par), f.mu = Vectorize(f.mu.par))
  }
  else {
    init <- c(lamb_par)
    p <- length(init)
    optimLH <- function(init) {
      lamb_par <- init[1:length(lamb_par)]
      f.lamb.par <- function(t) {
        abs(f.lamb(t, lamb_par))
      }
      f.mu.par <- function(t) {
        abs(f.mu(t, mu_par))
      }
      #the LH is now a sum of likelihoods for different trees, each tree can have a different total age (tot_time) 
      #and sampling proportion (f)
      lhsinpar=function(phylo, tot_time, f){likelihood_bd(phylo, tot_time, f.lamb.par, 
                                                       f.mu.par, f, cst.lamb = cst.lamb, cst.mu = cst.mu, 
                                                       expo.lamb = expo.lamb, expo.mu = expo.mu, dt = dt, 
                                                       cond = cond)}
      likelihoodlist=mapply(lhsinpar, phylolist, tot_timelist, flist)
      LH=Reduce("+",likelihoodlist)
      
      return(-LH)
    }
    temp <- suppressWarnings(optim(init, optimLH, method = meth))
    lamb.par <- temp$par[1:length(lamb_par)]
    f.lamb.par <- function(t) {
      f.lamb(t, lamb.par)
    }
    f.mu.par <- function(t) {
      f.mu(t, mu_par)
    }
    res <- list(model = "birth.death", LH = -temp$value, 
                aicc = 2 * temp$value + 2 * p + (2 * p * (p + 1))/(nobs - 
                                                                     p - 1), lamb_par = lamb.par, f.lamb = Vectorize(f.lamb.par))
  }
  class(res) <- "fit.bd"
  return(res)
}