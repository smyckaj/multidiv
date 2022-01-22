sim_bdmulti=function(n, lambda, mu, time, n.cores=1){
  #helper function shrinking problem to one argument
  rbdtreesim=function(lam){rbdtree(lam,mu,time)}
  #use parallel apply if it is loaded
  if ("parallel" %in% (.packages())){
    tlist=mclapply(rep(lambda,n), rbdtreesim,mc.cores=n.cores)
  } else {
    tlist=lapply(rep(lambda,n), rbdtreesim)
  }
  return(tlist)
}
