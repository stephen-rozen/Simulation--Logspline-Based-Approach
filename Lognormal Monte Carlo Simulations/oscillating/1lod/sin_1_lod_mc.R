## ----setup, include=FALSE----------------------------------------------------------------------------------------
source("sin_1_lod_setup.R")
library(future)
library(furrr)
library(purrr)
library(rlang)
# library(tidyverse)

future::plan(future::multisession, workers = max(1, future::availableCores() - 1))
# future::plan(future::multicore, workers = parallel::detectCores(logical = FALSE))
#error handling
safe_run <- possibly(one_run_sin_ILOD, otherwise = NULL)

compute_one_rep <- function(rep_id, params_df, outdir) {
  set.seed(sample.int(.Machine$integer.max, 1))  # reseed so retries differ
  out <- params_df %>%
    mutate(result = furrr::future_pmap(
      list(.data$samp_size, .data$GM, .data$GSD, .data$censored),
      safe_run,
      .progress = T,
      .options = furrr::furrr_options(seed = TRUE)
      
    )) %>%
    param_error()
  saveRDS(out, file.path(outdir, sprintf("rep_%03d.rds", rep_id)))
  invisible(TRUE)
}





## ----------------------------------------------------------------------------------------------------------------
sin_ILOD_small_params <- expand_grid(
  samp_size = seq(7, 19, 2),
  GM        = GM_seq,
  GSD       = GSD_seq,
  censored  = censored1_seq
)

sin_ILOD_small_params


## ----------------------------------------------------------------------------------------------------------------
execute_loop(60, sin_ILOD_small_params)



## ----------------------------------------------------------------------------------------------------------------
sin_ILOD_medium_params <- expand_grid(
  samp_size = seq(20, 90, 10),
  GM        = GM_seq,
  GSD       = GSD_seq,
  censored  = censored1_seq
)

sin_ILOD_medium_params


## ----------------------------------------------------------------------------------------------------------------
execute_loop(60, sin_ILOD_medium_params)


## ----------------------------------------------------------------------------------------------------------------
sin_ILOD_large_params <- expand_grid(
  samp_size = seq(100, 500, 100),
  GM        = GM_seq,
  GSD       = GSD_seq,
  censored  = censored1_seq
)

sin_ILOD_large_params


## ----------------------------------------------------------------------------------------------------------------
execute_loop(60, sin_ILOD_large_params)

