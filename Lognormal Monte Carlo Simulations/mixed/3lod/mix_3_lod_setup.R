## ----setup, include=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------
source("../../basic_functions.R")
source("../mix_specs.R")


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
triple_censor <- function(sim, fraction_censored_1,fraction_censored_2,fraction_censored_3) {
  
  # Split the simulated set into 3, representing data from 3 labs
  n <- length(sim)
  indices <- sample(n)
  split_indices <- split(indices, cut(seq_along(indices), 3, labels = FALSE))
  
  lab1 <- sim[split_indices[[1]]]
  lab2 <- sim[split_indices[[2]]]
  lab3 <- sim[split_indices[[3]]]
  
  # Assign each lab to their LOD
  labs_df <- tibble(lab = list(lab1, lab2, lab3), 
                        fraction_censored =  c(fraction_censored_1,fraction_censored_2,fraction_censored_3)
                        )
  
  # Created censored data sets for each lab
  labs_df <- labs_df %>%
    rowwise %>%
    mutate(cutoff = quantile(lab, probs = fraction_censored)) %>%
    mutate(censored_bool = list(lab < cutoff)) %>%
    mutate(censored_sim = list(replace(lab, censored_bool, cutoff)))
  
  
  # Return combined datasets from all labs
  list(
    censored_sim = flatten_dbl(labs_df$censored_sim),
    censored_bool = flatten_lgl(labs_df$censored_bool),
    cutoff = labs_df$cutoff
  )
}





## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
generate_IIln_IIILOD <- function(sample_size, GM1, GSD1, GM2, GSD2, 
                                    fraction_censored_1, fraction_censored_2,
                                 fraction_censored_3, p.mix){
   # 1. Run core sim + "truth" 
  sim <- rlnormMix(sample_size, 
                   meanlog1 = log(GM1), 
                   sdlog1 = log(GSD1), 
                   meanlog2 = log(GM2), 
                   sdlog2 = log(GSD2),
                   p.mix = p.mix)
  true_stats <- calculate_summary_stats(sim)
  
  # 2. Censor
  cen <- triple_censor(sim, fraction_censored_1,fraction_censored_2,fraction_censored_3)
  
  
  # 3. Return "dataset object", list
  dataset_object <- list(
  censored_sim     = cen$censored_sim,
  censored_bool    = cen$censored_bool,
  LOD  = cen$cutoff,
  true_list = true_stats
  )

  return(dataset_object)
}
q <- generate_IIln_IIILOD(30, 1, 1.2, 100, 1.1, .1, .3, .2, .5)
q


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

one_run_IIln_IIILOD <- function(sample_size, GM1, GSD1, GM2, GSD2, 
                                    fraction_censored_1, fraction_censored_2,
                                fraction_censored_3, p.mix,
                              distro = "lognormal", ..., .max_tries = 30){
  
  #generate dataset
  dataset_object <- generate_IIln_IIILOD(sample_size, GM1, GSD1, GM2, GSD2, 
                                    fraction_censored_1, fraction_censored_2,
                                    fraction_censored_3, p.mix)
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
                  fraction_censored_1,fraction_censored_2,
                  fraction_censored_3))
    message(paste("!!!!=======> Faliure", sample_size, GM1, GSD1, GM2, GSD2,
                  fraction_censored_1,fraction_censored_2,
                  fraction_censored_3, distro, .max_tries))
    return(Recall(sample_size, GM1, GSD1, GM2, GSD2,
                  fraction_censored_1, fraction_censored_2,fraction_censored_3, p.mix,
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

dee <- one_run_IIln_IIILOD(30, 1, 1.2, 100, 1.1, .1, .3, .2, .5)
dee



# LOD_estimation(censored_sim = q$censored_sim, censored_bool = q$censored_bool, LOD = q$LOD, forwardT = "log", reverseT = "exp", distro = "lognormal") %>%
#   error_calculation(q$true_list)
# 
# ROS_estimate(q$censored_sim, q$censored_bool)



