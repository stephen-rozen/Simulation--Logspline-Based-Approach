## ----setup, include=FALSE--------------------------------------------------------------------------------------------------------------------------------
source("sin_2_lod_setup.R")
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
#error handling
safe_run <- possibly(one_run_sin_IILOD, otherwise = NULL)



# compute_one_rep <- function(rep_id, params_df, outdir) {
#   set.seed(sample.int(.Machine$integer.max, 1))  # reseed so retries differ
#   out <- params_df %>%
#     mutate(result = furrr::future_pmap(
#       list(.data$samp_size, .data$GM, .data$GSD, .data$censored1, .data$censored2),
#       safe_run,
#       .options = furrr::furrr_options(seed = TRUE)
#     )) %>%
#     param_error()
#   saveRDS(out, file.path(outdir, sprintf("rep_%03d.rds", rep_id)))
#   invisible(TRUE)
# }


compute_one_rep <- function(rep_id, params_df, outdir) {
  set.seed(sample.int(.Machine$integer.max, 1))  # reseed so retries differ
  args <- with(params_df, list(samp_size, GM, GSD, censored1, censored2))
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
sin_IILOD_small_params <- expand_grid(
  samp_size = seq(7, 19, 2),
  GM        = GM_seq,
  GSD       = GSD_seq,
  censored1  = censored1_seq,
  censored2 = censored2_seq
  ) %>% 
  filter(censored1 < censored2)


## --------------------------------------------------------------------------------------------------------------------------------------------------------

execute_loop(60, sin_IILOD_small_params)



## --------------------------------------------------------------------------------------------------------------------------------------------------------
sin_IILOD_medium_params <- expand_grid(
  samp_size = seq(20, 90, 10),
  GM        = GM_seq,
  GSD       = GSD_seq,
  censored1  = censored1_seq,
  censored2 = censored2_seq
  ) %>% 
  filter(censored1 < censored2)


## --------------------------------------------------------------------------------------------------------------------------------------------------------

execute_loop(60, sin_IILOD_medium_params)



## --------------------------------------------------------------------------------------------------------------------------------------------------------
sin_IILOD_large_params <- expand_grid(
  samp_size = seq(100, 500, 100),
  GM        = GM_seq,
  GSD       = GSD_seq,
  censored1  = censored1_seq,
  censored2 = censored2_seq
  ) %>% 
  filter(censored1 < censored2)


## --------------------------------------------------------------------------------------------------------------------------------------------------------
execute_loop(60, sin_IILOD_large_params)

