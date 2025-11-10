## ----setup, include=FALSE----------------------------------------------------------------------------------------------------------------
source("../../basic_functions.R")
source("../one_dist_specs.R")


## ----------------------------------------------------------------------------------------------------------------------------------------
single_censor <- function(sim, fraction_censored) {
  # 1) LOD
  cutoff <- quantile(sim, probs = fraction_censored, na.rm = TRUE)
  
  # 2) Which values are below it?
  censored_bool <- sim < cutoff
  
  # 3) Create the “censored” version
  censored_sim <- sim
  censored_sim[censored_bool] <- cutoff
  
  # 4) Return
  list(
    censored_sim = censored_sim,
    censored_bool = censored_bool,
    cutoff = cutoff
  )
}


## ----------------------------------------------------------------------------------------------------------------------------------------
generate_ln_Ilod <- function(sample_size, GM, GSD, fraction_censored_1){
   # 1. Run core sim + "truth"
  sim <- rlnorm(n = sample_size, meanlog = log(GM), sdlog = log(GSD))
    #only need to apply log() on calls to rlnorm
  true_stats <- calculate_summary_stats(sim)

  # 2. Censor
  cen <- single_censor(sim, fraction_censored_1)


  # 3. Return "dataset object", list
  dataset_object <- list(
  censored_sim     = cen$censored_sim,
  censored_bool    = cen$censored_bool,
  LOD  = cen$cutoff,
  true_list = true_stats
  )

  return(dataset_object)
}
q <- generate_ln_Ilod(15, 1, 2, .1)
q

qbig <- generate_ln_Ilod(30, 1, 2, .1)

#substitution_estimate(q$censored_sim, q$censored_bool, "root2")



## ----------------------------------------------------------------------------------------------------------------------------------------

one_run_ln_ILOD <- function(sample_size, GM, GSD, fraction_censored, distro = "lognormal", ..., .max_tries = 30){
  
  #generate dataset
  dataset_object <- generate_ln_Ilod(sample_size, GM, GSD, fraction_censored)
  censored_sim <- dataset_object$censored_sim
  censored_bool <- dataset_object$censored_bool
  true_list <- dataset_object$true_list
  LOD <- dataset_object$LOD
  
  #fit the spline
  spline_obj<-spline_fit(censored_sim, censored_bool)
  
  # if NULL *and* we still have tries left, retry
  
  if (is.null(spline_obj) && .max_tries > 1 && sample_size > 13) {
    nooverwrite(censored_sim, 
                path = "converge_failures",
                name = paste(sample_size, GM, GSD, fraction_censored))
    message(paste("!!!!=======> Faliure", sample_size, GM, GSD, fraction_censored, distro, .max_tries))
    return(Recall(sample_size, GM, GSD, fraction_censored,
                  distro = distro, ..., .max_tries = .max_tries - 1))
  }
  # if NULL but no tries left, just continue with spline_obj = NULL


  
  #create estimate dataframe
  lod_est_df <- LOD_estimation(censored_sim = censored_sim, 
                               censored_bool = censored_bool, 
                               LOD = LOD, 
                               spline_obj = spline_obj, 
                               distro = distro, 
                               ...)
  
  final_df <- error_calculation(lod_est_df = lod_est_df, true_list = true_list)
  
  final_df

}

dee <- one_run_ln_ILOD(15, 1, 2, .1) 
dee



# LOD_estimation(censored_sim = q$censored_sim, censored_bool = q$censored_bool, LOD = q$LOD, forwardT = "log", reverseT = "exp", distro = "lognormal") %>%
#   error_calculation(q$true_list)
# 
# ROS_estimate(q$censored_sim, q$censored_bool)



