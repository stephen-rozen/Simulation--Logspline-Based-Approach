## ----setup, include=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------
source("../../basic_functions.R")
source("../mix_specs.R")


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
generate_IIln_ILOD <- function(sample_size, GM1, GSD1, GM2, GSD2, 
                                    fraction_censored_1, p.mix){
   # 1. Run core sim + "truth" 
  sim <- rlnormMix(sample_size, 
                   meanlog1 = log(GM1), 
                   sdlog1 = log(GSD1), 
                   meanlog2 = log(GM2), 
                   sdlog2 = log(GSD2),
                   p.mix = p.mix)
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
q <- generate_IIln_ILOD(30, 1, 1.2, 3, 1.1, .1, .5)
q


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

one_run_IIln_ILOD <- function(sample_size, GM1, GSD1, GM2, GSD2, 
                                    fraction_censored_1, p.mix,
                              distro = "lognormal", ..., .max_tries = 30){
  
  #generate dataset
  dataset_object <- generate_IIln_ILOD(sample_size, GM1, GSD1, GM2, GSD2, 
                                    fraction_censored_1, p.mix)
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
                name = paste(sample_size, GM1, GSD1, GM2, GSD2,
                  fraction_censored_1))
    message(paste("!!!!=======> Faliure", sample_size, GM1, GSD1, GM2, GSD2,
                  fraction_censored_1, distro, .max_tries))
    return(Recall(sample_size, GM1, GSD1, GM2, GSD2,
                  fraction_censored_1, p.mix,
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

dee <- one_run_IIln_ILOD(30, 1, 1.2, 3, 1.1, .1, .5)
dee



# LOD_estimation(censored_sim = q$censored_sim, censored_bool = q$censored_bool, LOD = q$LOD, forwardT = "log", reverseT = "exp", distro = "lognormal") %>%
#   error_calculation(q$true_list)
# 
# ROS_estimate(q$censored_sim, q$censored_bool)



