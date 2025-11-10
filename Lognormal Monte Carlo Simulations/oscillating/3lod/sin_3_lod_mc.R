## --------------------------------------------------------------------------------------------------------------------------------------------------------
source("sin_3_lod_setup.R")
library(future)
library(furrr)
library(purrr)
library(rlang)
# library(tidyverse)

if (future::supportsMulticore()) {
  # Linux: use forking
  future::plan(future::multicore, workers = parallel::detectCores(logical = FALSE))
} else {
  # Windows fallback
  future::plan(future::multisession, workers = max(1, future::availableCores() - 1))
}
# future::plan(future::multicore, workers = parallel::detectCores(logical = FALSE))
# current_plan <- future::plan()
# message("Using future plan: ", class(current_plan)[1], 
#         " with workers = ", future::nbrOfWorkers())

#error handling
safe_run <- possibly(one_run_sin_IIILOD, otherwise = NULL)

# compute_one_rep <- function(rep_id, params_df, outdir) {
#   set.seed(sample.int(.Machine$integer.max, 1))  # reseed so retries differ
#   out <- params_df %>%
#     mutate(result = furrr::future_pmap(
#       list(.data$samp_size, .data$GM, .data$GSD, .data$censored1, .data$censored2, .data$censored3),
#       safe_run,
#       .options = furrr::furrr_options(seed = TRUE)
#     )) %>%
#     param_error()
#   saveRDS(out, file.path(outdir, sprintf("rep_%03d.rds", rep_id)))
#   invisible(TRUE)
# }

compute_one_rep <- function(rep_id, params_df, outdir) {
  set.seed(sample.int(.Machine$integer.max, 1))  # reseed so retries differ
  args <- with(params_df, list(samp_size, GM, GSD, censored1, censored2, censored3))
  res  <- furrr::future_pmap(
    args,
    safe_run,
    .progress = T,
    .options = furrr::furrr_options(seed = TRUE, stdout = TRUE, scheduling = 1L)
  )
  out <- params_df
  out$result <- res
  out <- param_error(out)
  
  saveRDS(out, file.path(outdir, sprintf("rep_%03d.rds", rep_id)))
  invisible(TRUE)
}





## --------------------------------------------------------------------------------------------------------------------------------------------------------
sin_IIILOD_small_params <- expand_grid(
  samp_size = seq(7, 19, 2),
  GM        = GM_seq,
  GSD       = GSD_seq,
  censored1  = censored1_seq,
  censored2 = censored2_seq,
  censored3 = censored3_seq
  ) %>% 
  filter(censored1 < censored2) %>%
  filter(censored2 <= censored3)
sin_IIILOD_small_params


## --------------------------------------------------------------------------------------------------------------------------------------------------------

execute_loop(50, sin_IIILOD_small_params)



## --------------------------------------------------------------------------------------------------------------------------------------------------------
sin_IIILOD_medium_params <- expand_grid(
  samp_size = seq(20, 90, 10),
  GM        = GM_seq,
  GSD       = GSD_seq,
  censored1  = censored1_seq,
  censored2 = censored2_seq,
  censored3 = censored3_seq
  ) %>% 
  filter(censored1 < censored2) %>%
  filter(censored2 <= censored3)
sin_IIILOD_medium_params


## --------------------------------------------------------------------------------------------------------------------------------------------------------

execute_loop(50, sin_IIILOD_medium_params)

# sin_IIILOD_medium_20 <- sin_IIILOD_medium_params %>%
#   filter(samp_size == 20)
# try(execute_loop(50, sin_IIILOD_medium_20))
# 
# sin_IIILOD_medium_30 <- sin_IIILOD_medium_params %>%
#   filter(samp_size == 30)
# try(execute_loop(50, sin_IIILOD_medium_30))
# 
# sin_IIILOD_medium_40 <- sin_IIILOD_medium_params %>%
#   filter(samp_size == 40)
# try(execute_loop(50, sin_IIILOD_medium_40))
# 
# sin_IIILOD_medium_50 <- sin_IIILOD_medium_params %>%
#   filter(samp_size == 50)
# try(execute_loop(50, sin_IIILOD_medium_50))
# 
# sin_IIILOD_medium_60 <- sin_IIILOD_medium_params %>%
#   filter(samp_size == 60)
# try(execute_loop(50, sin_IIILOD_medium_60))
# 
# sin_IIILOD_medium_70 <- sin_IIILOD_medium_params %>%
#   filter(samp_size == 70)
# try(execute_loop(50, sin_IIILOD_medium_70))
# 
# sin_IIILOD_medium_80 <- sin_IIILOD_medium_params %>%
#   filter(samp_size == 80)
# try(execute_loop(50, sin_IIILOD_medium_80))
# 
# sin_IIILOD_medium_90 <- sin_IIILOD_medium_params %>%
#   filter(samp_size == 90)
# try(execute_loop(50, sin_IIILOD_medium_90))


## --------------------------------------------------------------------------------------------------------------------------------------------------------
sin_IIILOD_large_params <- expand_grid(
  samp_size = seq(100, 500, 100),
  GM        = GM_seq,
  GSD       = GSD_seq,
  censored1  = censored1_seq,
  censored2 = censored2_seq,
  censored3 = censored3_seq
  ) %>% 
  filter(censored1 < censored2) %>%
  filter(censored2 <= censored3)


## --------------------------------------------------------------------------------------------------------------------------------------------------------

execute_loop(50, sin_IIILOD_large_params)
# 
# sin_IIILOD_large_100 <- sin_IIILOD_large_params %>%
#   filter(samp_size == 100)
# try(execute_loop(50, sin_IIILOD_large_100))
# 
# sin_IIILOD_large_200 <- sin_IIILOD_large_params %>%
#   filter(samp_size == 200)
# try(execute_loop(50, sin_IIILOD_large_200))
# 
# sin_IIILOD_large_300 <- sin_IIILOD_large_params %>%
#   filter(samp_size == 300)
# try(execute_loop(50, sin_IIILOD_large_300))
# 
# sin_IIILOD_large_400 <- sin_IIILOD_large_params %>%
#   filter(samp_size == 400)
# try(execute_loop(50, sin_IIILOD_large_400))
# 
# sin_IIILOD_large_500 <- sin_IIILOD_large_params %>%
#   filter(samp_size == 500)
# try(execute_loop(50, sin_IIILOD_large_500))