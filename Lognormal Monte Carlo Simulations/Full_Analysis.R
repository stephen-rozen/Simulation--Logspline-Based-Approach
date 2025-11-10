## --------------------------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)


## --------------------------------------------------------------------------------------------------------------------------------------------------

# "param_error" or 'complete_case_error'
# 'complete_case_error' refers to error relative to full uncensored dataset
# applies your selection to the whole markdown file

err_type <- 'param_error'


## --------------------------------------------------------------------------------------------------------------------------------------------------


filtration <- function(results_df, n_range, grouping, stat_to_estimate, 
                       num_slice = 1, metric = 'rmse', e6_filter = F){
  # stat_to_estimate takes a vector with any of the following: 
    # c('mean', 'sd', 'gm', 'gsd')
  # grouping takes a vector with any of the following:
    # for single lod:
        # c( 'censored','samp_size','GSD')
    # for multiple lod:
        # c( 'censored1', 'censored2','censored3', 'samp_size','GSD')
    # more esoteric use:
        # grouping by average percent censored:
          # c('avg_censored','samp_size','GSD')
        # grouping combos of censoring rather than evaluating each individually:
          # list(c( 'censored1', 'censored2','censored3'),'samp_size','GSD')
  # n_range takes a 1 element list with a vector of desired sample size range. The following is a format example:
    # n =  c(14:19) ---> select spline-inclusive iterates
 
  
  if('avg_censored' %in% grouping){
    results_df <- results_df %>%
      mutate(avg_censored = rowMeans(across(contains("censored")), na.rm = TRUE))
  } 
  
  if(e6_filter){
    #filter if error on sample > 10^6
    results_df <- results_df %>%
      filter(complete_case_error < 10^6)
  }
  
  
  # --- bin GSD into 14 equal-width bins  ---
  if (any(grouping == "GSD") && dplyr::n_distinct(results_df$GSD, na.rm = TRUE) > 14) {
    g_min <- suppressWarnings(min(results_df$GSD, na.rm = TRUE))
    g_max <- suppressWarnings(max(results_df$GSD, na.rm = TRUE))
    # protect against degenerate ranges
    if (is.finite(g_min) && is.finite(g_max) && g_min < g_max) {
      brks <- seq(g_min, g_max, length.out = 15)  # 14 bins
    } else {
      # fall back to a single bin if range is degenerate
      brks <- c(g_min, g_max)
    }
    results_df <- results_df %>%
      dplyr::mutate(GSD_bin = cut(GSD, breaks = brks, include.lowest = TRUE, right = TRUE))

    # replace "GSD" with "GSD_bin" in the grouping vector
    grouping <- replace(grouping, grouping == "GSD", "GSD_bin")
  }
  # ----------------------------------------------------------
  
  results_df %>%
     filter(!is.na(.data[[err_type]]), 
            samp_size %in% n_range,
            est_id == stat_to_estimate
            ) %>%
    group_by(method, est_id, across(all_of(grouping))) %>%
    summarize(rmse =
    sqrt(mean(.data[[err_type]]**2)),
    bias =                        
      mean(.data[[err_type]])) %>%
    group_by(across(all_of(grouping)), est_id) %>%
  slice_min(abs(.data[[metric]]), n = num_slice, with_ties = T) %>%
    ungroup()
}

auto_analysis <- function(results_df, analysis_grid, num_slice = 1, metric = 'rmse', e6_filter = F){
   # num_slice specifies how many methods to return for each margin. 
    # Typical args are 1 (the best method), 3 (top 3), or 14 (all the methods)
  # metric takes 'rmse' or 'bias', chooses which to order by
    pmap(
      analysis_grid,
         function(n, grouper, stat_to_estimate){
            filtration(
                      results_df, 
                       n_range = n, 
                       grouping = grouper,  
                       stat_to_estimate = stat_to_estimate,
                       num_slice = num_slice,
                       metric = metric,
                       e6_filter = e6_filter
                      )
         }
    )
}


## --------------------------------------------------------------------------------------------------------------------------------------------------
comparison_filtration <- function(results_df, 
                                  n_range, grouping,
                                  stat_to_estimate,
                                  methods_to_compare,
                                  metric = 'rmse',
                                  e6_filter = F){ 
  # stat_to_estimate takes a vector with any of the following: 
    # c('mean', 'sd', 'gm', 'gsd')
 # grouping takes a vector with any of the following:
    # for single lod:
        # c( 'censored','samp_size','GSD')
    # for multiple lod:
        # c( 'censored1', 'censored2','censored3', 'samp_size','GSD')
    # more esoteric use:
        # grouping by average percent censored:
          # c('avg_censored','samp_size','GSD')
        # grouping combos of censoring rather than evaluating each individually:
          # list(c( 'censored1', 'censored2','censored3'),'samp_size','GSD')
  # n_range takes a list containing vectors of desired sample size range. The following is a format example:
    # n = list(
    #   # c(5:13),  # select non-spline iterates
    #     c(14:19)) # select spline-inclusive iterates
  # methods_to_compare takes a vector of your desired included methods. e.g:
    # list(c('km', 'spline', 'rspline'))
  if('avg_censored' %in% grouping){
    results_df <- results_df %>%
      mutate(avg_censored = rowMeans(across(contains("censored")), na.rm = TRUE))
  }  
  if(e6_filter){
    results_df <- results_df %>%
      filter(complete_case_error < 10^6)
  }
  
  
  # --- bin GSD into 14 equal-width bins if requested ---
  if (any(grouping == "GSD") && dplyr::n_distinct(results_df$GSD, na.rm = TRUE) > 14) {
    g_min <- suppressWarnings(min(results_df$GSD, na.rm = TRUE))
    g_max <- suppressWarnings(max(results_df$GSD, na.rm = TRUE))
    # protect against degenerate ranges
    if (is.finite(g_min) && is.finite(g_max) && g_min < g_max) {
      brks <- seq(g_min, g_max, length.out = 15)  # 14 bins
    } else {
      # fall back to a single bin if range is degenerate
      brks <- c(g_min, g_max)
    }
    results_df <- results_df %>%
      dplyr::mutate(GSD_bin = cut(GSD, breaks = brks, include.lowest = TRUE, right = TRUE))

    # replace "GSD" with "GSD_bin" in the grouping vector
    grouping <- replace(grouping, grouping == "GSD", "GSD_bin")
  }
  # ----------------------------------------------------------
  
  results_df %>%
     filter(
       !is.na(.data[[err_type]]), 
            .data[[err_type]] != Inf,
            samp_size %in% n_range,
            method %in% methods_to_compare,
            est_id == stat_to_estimate
       ) %>%
    group_by(method, est_id, across(all_of(grouping))) %>%
    summarize(rmse =
    sqrt(mean(.data[[err_type]]**2)),
    bias =                        
      mean(.data[[err_type]])) %>%
    group_by(across(all_of(grouping)), est_id) %>%
  slice_min(abs(.data[[metric]]),
            n = length(methods_to_compare), 
            with_ties = T
            ) %>%
    ungroup()
}

comparison_auto_analysis <- function(results_df, comparison_grid, metric = 'rmse', e6_filter = F){
  # metric takes 'rmse' or 'bias', chooses which to order by
    pmap(
      comparison_grid,
      function(n, grouper, stat_to_estimate, methods_to_compare){
            comparison_filtration(
                                  results_df, 
                                  n_range = n,
                                  grouping = grouper,  
                                  stat_to_estimate = stat_to_estimate, 
                                  methods_to_compare = methods_to_compare,
                                  metric = metric, 
                                   e6_filter = e6_filter
                                  )
      }
    )
}


## --------------------------------------------------------------------------------------------------------------------------------------------------
best_overall <- function(results_df, 
                         metric = 'rmse',
                         num_slice = 3,
                         max_samp_size = max(results_df$samp_size), 
                         min_samp_size = min(results_df$samp_size),
                         e6_filter = F){
  # metric takes 'rmse' or 'bias'
  # min/max_samp_size is the INCLUSIVE bounds of sample size
  # VALUE: Returns best method for each stat (mean, sd, gm, gsd)
   if(e6_filter){
    results_df <- results_df %>%
      filter(complete_case_error < 10^6)
   }
  
  results_df %>%
  filter(!is.na(.data[[err_type]]), 
         samp_size >= min_samp_size,
         samp_size <= max_samp_size) %>%
  group_by(method, est_id) %>%
  summarize(
    bias = mean(.data[[err_type]]),
    rmse = sqrt(mean(.data[[err_type]]**2))
    ) %>%
  group_by(est_id) %>%
  slice_min(abs(.data[[metric]]), n = num_slice) %>%
  select(method,est_id, rmse, bias)  %>%
    ungroup()
}


## --------------------------------------------------------------------------------------------------------------------------------------------------
extremes <- function(results_df, methods = NULL){
  # methods takes vector of methods, e.g. c('spline', 'rspline').
  # If no arg is passed, it returns extremes of all methods
  # to exclude, follow this: extremes(df, methods = setdiff(unique(df$method), "spline"))
  results_df %>%
  dplyr::filter(!is.na(.data[[err_type]])) %>%
    dplyr::filter(if (is.null(methods)) TRUE else .data$method %in% methods) %>%
    dplyr::arrange(dplyr::desc(abs(.data[[err_type]]))) %>%
    slice_head(n = 1e3)
}


## --------------------------------------------------------------------------------------------------------------------------------------------------
spline_fails <- function(results_df){
   results_df %>%
  filter(
    is.na(.data[[err_type]]) | is.infinite(.data[[err_type]]),
         method %in% c('spline', 'rspline')
    ) %>%
   filter(samp_size > 14)
}


## --------------------------------------------------------------------------------------------------------------------------------------------------
basic_analysis <- function(results_df, analysis_grid, metric = 'rmse', 
                           e6_filter = F, n_slice_auto = 20, n_slice_bests = 20,
                          best_max_n = max(results_df$samp_size),
                          best_min_n = min(results_df$samp_size), name){
  #specify best_max_n

  comparison_grid_km <- expand_grid(analysis_grid, 
                                  methods_to_compare = list(c('km', 'spline', 'rspline')))
  comparison_grid_ros <-  expand_grid(analysis_grid, 
                                  methods_to_compare = list(c('ros', 'esros')))
    
  
  analysis_lst <- list(
    
    overall_best = best_overall(results_df, 
                         metric,
                         num_slice = n_slice_bests,
                         max_samp_size = best_max_n, 
                         min_samp_size = best_min_n,
                         e6_filter = e6_filter),
    
    marginal = auto_analysis(results_df, 
                             analysis_grid, 
                             num_slice = n_slice_auto,
                             metric = metric, 
                             e6_filter = e6_filter) %>%
      list_rbind(),
    
    np_comparison = comparison_auto_analysis(results_df, comparison_grid_km, metric = metric, e6_filter = e6_filter)%>%
      list_rbind(),
    
    ros_comparison = comparison_auto_analysis(results_df, comparison_grid_ros, metric = metric, e6_filter = e6_filter)%>%
      list_rbind(),
    
    spline_faliures = spline_fails(results_df),
    
    extreme_spline = extremes(results_df, c('spline', 'rspline')),
    
    extreme_all = extremes(results_df)
    
    
  )
  
  # var_nm  <- deparse(substitute(name))
  base_nm <- if (grepl("_results$", name)) sub("_results$", "_basic_analysis", name)
             else paste0(name, "_basic_analysis")
  file_nm <- paste0(base_nm, ".rds")
  
  saveRDS(analysis_lst, file = file.path("analysis", file_nm))
  invisible(analysis_lst)
}




## --------------------------------------------------------------------------------------------------------------------------------------------------

print_all <- function(lst) {
  walk(lst, function(df) {
      print(head(df, 100))
  }
  )
}


## --------------------------------------------------------------------------------------------------------------------------------------------------

all_bests <- function(results_df, metric = 'rmse', e6_filter = F, 
                      best_max_n = max(results_df$samp_size),
                      best_min_n = min(results_df$samp_size), name){
  overall_best <- best_overall(results_df, 
               metric,
               max_samp_size = best_max_n, 
               min_samp_size = best_min_n,
               slice_min = 1000,
               e6_filter = e6_filter)
  
  
  base_nm <- if (grepl("_results$", name)) sub("_results$", "_full_bests", name)
  else paste0(name, "_full_bests")
  file_nm <- paste0(base_nm, ".rds")
  
  saveRDS(analysis_lst, file = file.path("analysis", file_nm))
  invisible(analysis_lst)
}
  
  