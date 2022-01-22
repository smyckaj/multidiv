fit_envmulti=function (phylolist, env_data, tot_timelist, f.lamb, f.mu, lamb_par, 
          mu_par, df = NULL, flist, meth = "Nelder-Mead", cst.lamb = FALSE, 
          cst.mu = FALSE, expo.lamb = FALSE, expo.mu = FALSE, fix.mu = FALSE, 
          dt = 0, cond = "crown") 
{
  if (is.null(df)) {
    df <- smooth.spline(x = env_data[, 1], env_data[, 2])$df
  }
  spline_result <- sm.spline(env_data[, 1], env_data[, 2], 
                             df = df)
  env_func <- function(t) {
    predict(spline_result, t)
  }
  lower_bound_control <- 0.1
  upper_bound_control <- 0.1
  lower_bound <- min(env_data[, 1])
  upper_bound <- max(env_data[, 1])
  time_tabulated <- seq(from = lower_bound * (1 - lower_bound_control), 
                        to = upper_bound * (1 + upper_bound_control), length.out = 1 + 
                          1e+06)
  env_tabulated <- env_func(time_tabulated)
  env_func_tab <- function(t) {
    b <- upper_bound * (1 + upper_bound_control)
    a <- lower_bound * (1 - lower_bound_control)
    n <- length(env_tabulated) - 1
    index <- 1 + as.integer((t - a) * n/(b - a))
    return(env_tabulated[index])
  }
  f.lamb.env <- function(t, y) {
    f.lamb(t, env_func_tab(t), y)
  }
  f.mu.env <- function(t, y) {
    f.mu(t, env_func_tab(t), y)
  }
  #here we use fit_bdmulti instead of fit_bd
  res <- fit_bdmulti(phylolist, tot_timelist, f.lamb.env, f.mu.env, lamb_par, 
                mu_par, flist, meth=meth, cst.lamb, cst.mu, expo.lamb, expo.mu, 
                fix.mu, dt, cond)
  res$model <- "environmental birth death"
  res$f.lamb <- function(t) {
    f.lamb(t, env_func_tab(t), res$lamb_par)
  }
  if (fix.mu == FALSE) {
    res$f.mu <- function(t) {
      f.mu(t, env_func_tab(t), res$mu_par)
    }
  }
  class(res) <- "fit.env"
  return(res)
}