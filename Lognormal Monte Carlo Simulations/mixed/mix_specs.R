## --------------------------------------------------------------------------------------------------------------------------------------------------------
# library(tidyverse)



## ----Universal_Param_Ranges------------------------------------------------------------------------------------------------------------------------------
GM1_seq <- c(.05, .5)
GM2_seq <- c(1, 2, 10)
GSD1_seq <- c(seq(1.1, 2, .3), seq(2.5, 7, .5))
GSD2_seq <- c(1.2, 1.6, 2, 3, 5, 7)
p.mix_seq <- seq.int(3L, 7L) / 10
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



## ----param_error_mix-------------------------------------------------------------------------------------------------------------------------------------
param_error_mix <- function(trial_result_df){
  trial_result_df %>%
    unnest(result) %>%        # expand the 'result' list-column            
    mutate(
      p.mix = p.mix,
      mu1    = log(GM1),
      sigma1 = log(GSD1),
      mu2 = log(GM2),
      sigma2 = log(GSD2),
      mean1 = exp(mu1 + sigma1^2/2),
      mean2 = exp(mu2 + sigma2^2/2),
      sd1 = sqrt((exp(sigma1^2) - 1) * exp(2*mu1 + sigma1^2)),
      sd2 = sqrt((exp(sigma2^2) - 1) * exp(2*mu2 + sigma2^2)),
      GM = GM1^(1-p.mix)*GM2^(p.mix),
      GSD = exp(sqrt(
                            (1 - p.mix) * sigma1^2 +
                              p.mix       * sigma2^2 +
                            (1 - p.mix) * p.mix * (mu1 - mu2)^2
                            )),
      true_param = case_when(
        est_id == "mean" ~ (1-p.mix)*mean1 + p.mix*mean2,
        
        est_id == "sd" ~ sqrt(
        # mixture second moment
        (1 - p.mix)*(sd1^2 + mean1^2) +
        p.mix    *(sd2^2 + mean2^2)
        # minus square of the mixture mean
        - ((1 - p.mix)*mean1 + p.mix*mean2)^2),
        
        est_id == "gm"   ~ GM,
        est_id == "gsd"  ~  GSD,
        TRUE             ~ NA_real_
      ),
      param_error = est - true_param
    ) %>%
    # #flatten unnecessary columns
    # select(samp_size, GSD, censored, method, est_id, complete_case_error, param_error) 
    
    #flatten unnecessary columns
    # has columns:
    # samp_size, GSD, censored(maybe multiple), method, est_id,         
    # complete_case_error, param_error, 
    #  -est, -true_stat, -true_param, -GM, -GM1, 
    # -GM2, -p.mix, -mu1, -mu2, -sigma1, -sigma2, -mean1, -mean2, 
    # -sd1, -sd2, -GSD1, -GSD2
  
    select(-est, -true_stat, -true_param, -GM, -GM1, 
           -GM2, -p.mix, -mu1, -mu2, -sigma1, -sigma2, -mean1, -mean2, 
           -sd1, -sd2, -GSD1, -GSD2)
  
}


## ----constrainGSD----------------------------------------------------------------------------------------------------------------------------------------
constrainGSD <- function(param_df){
  param_df %>%
  mutate(
      p.mix = p.mix,
      mu1    = log(GM1),
      sigma1 = log(GSD1),
      mu2 = log(GM2),
      sigma2 = log(GSD2),
      mean1 = exp(mu1 + sigma1^2/2),
      mean2 = exp(mu2 + sigma2^2/2),
      sd1 = sqrt((exp(sigma1^2) - 1) * exp(2*mu1 + sigma1^2)),
      sd2 = sqrt((exp(sigma2^2) - 1) * exp(2*mu2 + sigma2^2)),
      GM = GM1^(1-p.mix)*GM2^(p.mix),
      GSD = exp(sqrt(
                            (1 - p.mix) * sigma1^2 +
                              p.mix       * sigma2^2 +
                            (1 - p.mix) * p.mix * (mu1 - mu2)^2
      ))
  ) %>%
  filter(GSD<=7) %>%
    select(-mu1, -mu2, -sigma1, -sigma2, -mean1, -mean2, 
           -sd1, -sd2, -GM, -GSD)
}

