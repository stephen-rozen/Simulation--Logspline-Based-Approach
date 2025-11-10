## ----setup, include=FALSE--------------------------------------------------------------------------------------------------------------------------------
# library(tidyverse)


## ----Universal_Param_Ranges------------------------------------------------------------------------------------------------------------------------------
GM_seq <- c(.05, .5, 1, 2, 10)
GSD_seq <- c(seq(1.1, 2, .3), seq(2.5, 7, .5))
censored1_seq <- c(.025, .05, .1, seq(.15, .45, .1))
censored2_seq <- c(.025, .1, .25, .45)
censored3_seq <- censored2_seq


## ----method_specs----------------------------------------------------------------------------------------------------------------------------------------


LOD_estimation <- function(censored_sim, censored_bool, LOD, spline_obj, distro = "lognormal", ...){
  
  
  # core args
  core_args <- list(
  censored_sim     = censored_sim,
  censored_bool    = censored_bool)

# full method specs:
method_specs <- list(
  ros = list(
    fn   = ROS_estimate,
    args = core_args
  ),
   esros = list(
    fn   = ES_ROS_estimate,
    args = core_args
  ),
  mle = list(
    fn   = MLE_estimate,
    args = c(core_args,
             list(method     = "mle")
             )
  ),
   rmle = list(
    fn   = MLE_estimate,
    args = c(core_args,
             list(method     = "rmle")
             )
  ),
  bcmle = list(
    fn   = MLE_estimate,
    args = c(core_args,
             list(method     = "bcmle")
             )
  ),
  km = list(
    fn   = KM_estimate,
    args = c(core_args,
             list(...))
  ),
  root2 = list(
    fn   = substitution_estimate,
    args = c(core_args,
             list(method     = "root2"))
  ),
  half = list(
    fn   = substitution_estimate,
    args = c(core_args,
             list(method     = "half"))
  ),
  LOD = list(
    fn   = substitution_estimate,
    args = c(core_args,
             list(method     = "LOD"))
  ),
  zero = list(
    fn   = substitution_estimate,
    args = c(core_args,
             list(method     = "zero"))
  ),
  spline = list(
    fn   = spline_estimate,
    args = c(core_args, 
             list(spline_obj))
  ),
  rspline = list(
    fn   = robust_spline,
    args = c(core_args, 
             list(spline_obj))
  ),
  wlod = list(
    fn   = wlod2_estimate,
    args = c(core_args,
             list(LOD = LOD))
  ),
  beta = list(
    fn   = beta_sub,
    args = c(core_args,
             list(LOD = LOD))
  )
)

  # # Warning
  # if (!LOD_method %in% names(method_specs)) {
  #   stop("Unknown LOD_method: ", LOD_method)
  # }

  
  
  # create df
  df <- data.frame(method = names(method_specs)) %>%
    rowwise %>%
    mutate(est = list(do.call(method_specs[[method]]$fn, method_specs[[method]]$args)))
  df %>%
    unnest_longer(est)
  
}

#spl <- spline_fit(q$censored_sim, q$censored_bool)
# q2 <- LOD_estimation(q$censored_sim, q$censored_bool, q$LOD, spl)
# q2



## ----param_error-----------------------------------------------------------------------------------------------------------------------------------------
param_error <- function(trial_result_df){
  trial_result_df %>%
    unnest(result) %>%        # expand the 'result' list-column            
    mutate(
      mu    = log(GM),
      sigma = log(GSD),
      true_param = case_when(
        est_id == "mean" ~ exp(mu + sigma^2/2),
        est_id == "sd"   ~ sqrt((exp(sigma^2) - 1) * exp(2*mu + sigma^2)),
        est_id == "gm"   ~ GM,
        est_id == "gsd"  ~ GSD,
        TRUE             ~ NA_real_
      ),
      param_error = est - true_param
    ) %>%
    #flatten unnecessary columns
    # has columns:
    # samp_size, GSD, censored(maybe multiple), method, est_id,         
    # complete_case_error, param_error, -true_param, -est, 
    # -true_stat, -true_param, -mu, -sigma, -GM
    select(-true_param, -est, -true_stat, -true_param, -mu, -sigma, -GM)
}

