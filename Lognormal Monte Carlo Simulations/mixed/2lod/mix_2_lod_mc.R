## ----setup, include=FALSE--------------------------------------------------------------------------------------------------------------------------------
source("mix_2_lod_setup.R")
library(future)
library(furrr)
library(purrr)
library(rlang)
library(tidyverse)

future::plan(future::multisession, workers = max(1, future::availableCores() - 1))

#error handling
safe_run <- possibly(one_run_IIln_IILOD, otherwise = NULL)
# 
compute_one_rep <- function(rep_id, params_df, outdir) {
  set.seed(sample.int(.Machine$integer.max, 1))  # reseed so retries differ
  args <- with(params_df, list(samp_size, GM1, GSD1, GM2, GSD2, censored1, censored2, p.mix))
  res  <- furrr::future_pmap(
    args,
    safe_run,
    .progress = T,
    .options = furrr::furrr_options(seed = TRUE, stdout = TRUE, scheduling = 1L)
  )
  out <- params_df
  out$result <- res
  out <- param_error_mix(out)

  saveRDS(out, file.path(outdir, sprintf("rep_%03d.rds", rep_id)))
  invisible(TRUE)
}

#error handling
# safe_run <- possibly(one_run_IIln_IILOD, otherwise = NULL)
# 
# compute_one_rep <- function(rep_id, params_df, outdir) {
#   set.seed(sample.int(.Machine$integer.max, 1))  # reseed so retries differ
#   args <- with(params_df, list(samp_size, GM1, GSD1, GM2, GSD2, censored1, censored2, p.mix))
#   res  <- pmap(
#     args, 
#     safe_run
#   )
#   out <- params_df
#   out$result <- res
#   out <- param_error_mix(out)
#   
#   saveRDS(out, file.path(outdir, sprintf("rep_%03d.rds", rep_id)))
#   invisible(TRUE)
# }


## --------------------------------------------------------------------------------------------------------------------------------------------------------
param_args <- list(
  GM1 = GM1_seq,
  GM2 = GM2_seq,
  GSD1 = GSD1_seq,
  GSD2 = GSD2_seq,
  censored1  = censored1_seq,
  censored2 = censored2_seq
  )


## --------------------------------------------------------------------------------------------------------------------------------------------------------
# files<-list.files("IIln_IILOD_medium_results", full.names = TRUE, pattern = "\\.rds$")
#  results <- purrr::map(files, readRDS) %>% purrr::list_rbind()
#  saveRDS(results, "IIln_IILOD_medium_results_fin.rds")
#  unlink("mc_checkpoints", recursive = TRUE)


## --------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(123)
IIln_IILOD_small_params <- expand_grid(
  samp_size = seq(5, 19, 2),
  !!!param_args
) %>% 
  filter(censored1 < censored2) %>%
  mutate(p.mix = sample(p.mix_seq, size = n(), replace = T)) %>%
  constrainGSD()

IIln_IILOD_small_params


## --------------------------------------------------------------------------------------------------------------------------------------------------------

execute_loop(15, IIln_IILOD_small_params)




## --------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(123)
IIln_IILOD_medium_params <- expand_grid(
  samp_size = seq(20, 90, 10),
  !!!param_args
) %>% 
  filter(censored1 < censored2) %>%
  mutate(p.mix = sample(p.mix_seq, size = n(), replace = T))%>%
  constrainGSD()

IIln_IILOD_medium_params


## --------------------------------------------------------------------------------------------------------------------------------------------------------

execute_loop(15, IIln_IILOD_medium_params)



## --------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(123)
IIln_IILOD_large_params <- expand_grid(
  samp_size = seq(100, 500, 100),
  !!!param_args
) %>% 
  filter(censored1 < censored2) %>%
  mutate(p.mix = sample(p.mix_seq, size = n(), replace = T))%>%
  constrainGSD()

IIln_IILOD_large_params



## --------------------------------------------------------------------------------------------------------------------------------------------------------
# execute_loop(15, IIln_IILOD_large_params)

IIln_IILOD_large_100_.3 <- IIln_IILOD_large_params %>%
  filter(samp_size == 100, p.mix == .3)
try(execute_loop(15, IIln_IILOD_large_100_.3))

IIln_IILOD_large_100_.4 <- IIln_IILOD_large_params %>%
  filter(samp_size == 100, p.mix == .4)
try(execute_loop(15, IIln_IILOD_large_100_.4))

IIln_IILOD_large_100_.5 <- IIln_IILOD_large_params %>%
  filter(samp_size == 100, p.mix == .5)
try(execute_loop(15, IIln_IILOD_large_100_.5))

IIln_IILOD_large_100_.6 <- IIln_IILOD_large_params %>%
  filter(samp_size == 100, p.mix == .6)
try(execute_loop(15, IIln_IILOD_large_100_.6))

IIln_IILOD_large_100_.7 <- IIln_IILOD_large_params %>%
  filter(samp_size == 100, p.mix == .7)
try(execute_loop(15, IIln_IILOD_large_100_.7))


IIln_IILOD_large_200_.3 <- IIln_IILOD_large_params %>%
  filter(samp_size == 200, p.mix == .3)
try(execute_loop(15, IIln_IILOD_large_200_.3))

IIln_IILOD_large_200_.4 <- IIln_IILOD_large_params %>%
  filter(samp_size == 200, p.mix == .4)
try(execute_loop(15, IIln_IILOD_large_200_.4))

IIln_IILOD_large_200_.5 <- IIln_IILOD_large_params %>%
  filter(samp_size == 200, p.mix == .5)
try(execute_loop(15, IIln_IILOD_large_200_.5))

IIln_IILOD_large_200_.6 <- IIln_IILOD_large_params %>%
  filter(samp_size == 200, p.mix == .6)
try(execute_loop(15, IIln_IILOD_large_200_.6))

IIln_IILOD_large_200_.7 <- IIln_IILOD_large_params %>%
  filter(samp_size == 200, p.mix == .7)
try(execute_loop(15, IIln_IILOD_large_200_.7))


IIln_IILOD_large_300_.3 <- IIln_IILOD_large_params %>%
  filter(samp_size == 300, p.mix == .3)
try(execute_loop(15, IIln_IILOD_large_300_.3))

IIln_IILOD_large_300_.4 <- IIln_IILOD_large_params %>%
  filter(samp_size == 300, p.mix == .4)
try(execute_loop(15, IIln_IILOD_large_300_.4))

IIln_IILOD_large_300_.5 <- IIln_IILOD_large_params %>%
  filter(samp_size == 300, p.mix == .5)
try(execute_loop(15, IIln_IILOD_large_300_.5))

IIln_IILOD_large_300_.6 <- IIln_IILOD_large_params %>%
  filter(samp_size == 300, p.mix == .6)
try(execute_loop(15, IIln_IILOD_large_300_.6))

IIln_IILOD_large_300_.7 <- IIln_IILOD_large_params %>%
  filter(samp_size == 300, p.mix == .7)
try(execute_loop(15, IIln_IILOD_large_300_.7))


IIln_IILOD_large_400_.3 <- IIln_IILOD_large_params %>%
  filter(samp_size == 400, p.mix == .3)
try(execute_loop(15, IIln_IILOD_large_400_.3))

IIln_IILOD_large_400_.4 <- IIln_IILOD_large_params %>%
  filter(samp_size == 400, p.mix == .4)
try(execute_loop(15, IIln_IILOD_large_400_.4))

IIln_IILOD_large_400_.5 <- IIln_IILOD_large_params %>%
  filter(samp_size == 400, p.mix == .5)
try(execute_loop(15, IIln_IILOD_large_400_.5))

IIln_IILOD_large_400_.6 <- IIln_IILOD_large_params %>%
  filter(samp_size == 400, p.mix == .6)
try(execute_loop(15, IIln_IILOD_large_400_.6))

IIln_IILOD_large_400_.7 <- IIln_IILOD_large_params %>%
  filter(samp_size == 400, p.mix == .7)
try(execute_loop(15, IIln_IILOD_large_400_.7))


IIln_IILOD_large_500_.3 <- IIln_IILOD_large_params %>%
  filter(samp_size == 500, p.mix == .3)
try(execute_loop(15, IIln_IILOD_large_500_.3))

IIln_IILOD_large_500_.4 <- IIln_IILOD_large_params %>%
  filter(samp_size == 500, p.mix == .4)
try(execute_loop(15, IIln_IILOD_large_500_.4))

IIln_IILOD_large_500_.5 <- IIln_IILOD_large_params %>%
  filter(samp_size == 500, p.mix == .5)
try(execute_loop(15, IIln_IILOD_large_500_.5))

IIln_IILOD_large_500_.6 <- IIln_IILOD_large_params %>%
  filter(samp_size == 500, p.mix == .6)
try(execute_loop(15, IIln_IILOD_large_500_.6))

IIln_IILOD_large_500_.7 <- IIln_IILOD_large_params %>%
  filter(samp_size == 500, p.mix == .7)
try(execute_loop(15, IIln_IILOD_large_500_.7))



