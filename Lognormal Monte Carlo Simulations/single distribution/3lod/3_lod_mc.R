## --------------------------------------------------------------------------------------------------------------------------------------------------------
source("3_lod_setup.R")
library(future)
library(furrr)
library(purrr)
library(rlang)
library(tidyverse)

future::plan(future::multisession, workers = max(1, future::availableCores() - 1))
#error handling
safe_run <- possibly(one_run_ln_IIILOD, otherwise = NULL)

compute_one_rep <- function(rep_id, params_df, outdir) {
  set.seed(sample.int(.Machine$integer.max, 1))  # reseed so retries differ
  out <- params_df %>%
    mutate(result = furrr::future_pmap(
      list(.data$samp_size, .data$GM, .data$GSD, .data$censored1, .data$censored2, .data$censored3),
      safe_run,
      .progress = T,
      .options = furrr::furrr_options(seed = TRUE)
    )) %>%
    param_error()
  saveRDS(out, file.path(outdir, sprintf("rep_%03d.rds", rep_id)))
  invisible(TRUE)
}





## --------------------------------------------------------------------------------------------------------------------------------------------------------
ln_IIILOD_small_params <- expand_grid(
  samp_size = seq(5, 19, 2),
  GM        = GM_seq,
  GSD       = GSD_seq,
  censored1  = censored1_seq,
  censored2 = censored2_seq,
  censored3 = censored3_seq
  ) %>% 
  filter(censored1 < censored2) %>%
  filter(censored2 <= censored3)
ln_IIILOD_small_params


## --------------------------------------------------------------------------------------------------------------------------------------------------------

execute_loop(50, ln_IIILOD_small_params)



## --------------------------------------------------------------------------------------------------------------------------------------------------------
ln_IIILOD_medium_params <- expand_grid(
  samp_size = seq(20, 90, 10),
  GM        = GM_seq,
  GSD       = GSD_seq,
  censored1  = censored1_seq,
  censored2 = censored2_seq,
  censored3 = censored3_seq
  ) %>% 
  filter(censored1 < censored2) %>%
  filter(censored2 <= censored3)
ln_IIILOD_medium_params


## --------------------------------------------------------------------------------------------------------------------------------------------------------

execute_loop(50, ln_IIILOD_medium_params)



## --------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(123)
ln_IIILOD_large_params <- expand_grid(
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

# ln_IIILOD_large_100 <- ln_IIILOD_large_params %>%
#   filter(samp_size == 100)
# try(execute_loop(50, ln_IIILOD_large_100))
# 
# ln_IIILOD_large_200 <- ln_IIILOD_large_params %>%
#   filter(samp_size == 200)
# try(execute_loop(50, ln_IIILOD_large_200))
# 
# ln_IIILOD_large_300 <- ln_IIILOD_large_params %>%
#   filter(samp_size == 300)
# try(execute_loop(50, ln_IIILOD_large_300))
# 
# ln_IIILOD_large_400 <- ln_IIILOD_large_params %>%
#   filter(samp_size == 400)
# try(execute_loop(50, ln_IIILOD_large_400))
# 
# ln_IIILOD_large_500 <- ln_IIILOD_large_params %>%
#   filter(samp_size == 500)
# try(execute_loop(50, ln_IIILOD_large_500))
# 
# 
