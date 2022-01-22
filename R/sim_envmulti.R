sim_envmulti=function(n, env_data, lambda, mu, alpha, time, n.cores=1){
  #helper functions shrinking problem to one argument
  f.lamb=function(t,x,y){y[1]*exp(x*y[2])}
  f.mu=function(t,x,y){y[1]}
  renvtreesim=function(lam){
    #the control that trees with extant species are produced, the original return.all.extinct control does not work properly
    repeat{
      t=sim_env_bd(env_data,f.lamb,f.mu,lamb_par=c(lambda, alpha), mu_par=mu, time.stop=time,return.all.extinct=F, prune.extinct=T)$tree
      if (inherits(t, "phylo")) break
    }
  t
  }
  #use parallel apply if it is loaded
  if ("parallel" %in% (.packages())){
    tlist=mclapply(rep(lambda,n), renvtreesim,mc.cores=n.cores)
  } else {
    tlist=lapply(rep(lambda,n), renvtreesim)
  }

  return(tlist)
}
