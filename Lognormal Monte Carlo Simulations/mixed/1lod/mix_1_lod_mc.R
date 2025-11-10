## ----setup, include=FALSE--------------------------------------------------------------------------------------------------------------------------------
source("mix_1_lod_setup.R")
library(future)
library(furrr)
library(purrr)
library(rlang)
library(tidyverse)

future::plan(future::multisession, workers = max(1, future::availableCores() - 1))

#error handling
safe_run <- possibly(one_run_IIln_ILOD, otherwise = NULL)

compute_one_rep <- function(rep_id, params_df, outdir) {
  set.seed(sample.int(.Machine$integer.max, 1))  # reseed so retries differ
  out <- params_df %>%
    mutate(result = furrr::future_pmap(
      list(.data$samp_size, .data$GM1, .data$GSD1, .data$GM2, .data$GSD2, .data$censored, .data$p.mix),
      safe_run,
      .progress = T,
      .options = furrr::furrr_options(seed = TRUE)
    )) %>%
    param_error_mix()
  saveRDS(out, file.path(outdir, sprintf("rep_%03d.rds", rep_id)))
  invisible(TRUE)
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
param_args <- list(
  GM1 = GM1_seq,
  GM2 = GM2_seq,
  GSD1 = GSD1_seq,
  GSD2 = GSD2_seq,
  censored  = censored1_seq
  )


## --------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(123)
IIln_ILOD_small_params <- expand_grid(
  samp_size = seq(5, 19, 2),
  !!!param_args
)  %>%
  mutate(p.mix = sample(p.mix_seq, size = n(), replace = T))%>%
  constrainGSD()

IIln_ILOD_small_params


## --------------------------------------------------------------------------------------------------------------------------------------------------------

execute_loop(25, IIln_ILOD_small_params)




## --------------------------------------------------------------------------------------------------------------------------------------------------------
# files<-list.files("mc_checkpoints", full.names = TRUE, pattern = "\\.rds$")
#  results <- purrr::map(files, readRDS) %>% purrr::list_rbind()
#  saveRDS(results, "IIln_ILOD_small_results.rds")
#  unlink("mc_checkpoints", recursive = TRUE)


## --------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(123)
IIln_ILOD_medium_params <- expand_grid(
  samp_size = seq(20, 90, 10),
  !!!param_args
)  %>%
  mutate(p.mix = sample(p.mix_seq, size = n(), replace = T))%>%
  constrainGSD()

IIln_ILOD_medium_params


## --------------------------------------------------------------------------------------------------------------------------------------------------------
IIln_ILOD_medium_20 <- IIln_ILOD_medium_params %>%
  filter(samp_size == 20)
try(execute_loop(25, IIln_ILOD_medium_20))

IIln_ILOD_medium_30 <- IIln_ILOD_medium_params %>%
  filter(samp_size == 30)
try(execute_loop(25, IIln_ILOD_medium_30))

IIln_ILOD_medium_40 <- IIln_ILOD_medium_params %>%
  filter(samp_size == 40)
try(execute_loop(25, IIln_ILOD_medium_40))

IIln_ILOD_medium_50 <- IIln_ILOD_medium_params %>%
  filter(samp_size == 50)
try(execute_loop(25, IIln_ILOD_medium_50))

IIln_ILOD_medium_60 <- IIln_ILOD_medium_params %>%
  filter(samp_size == 60)
try(execute_loop(25, IIln_ILOD_medium_60))

IIln_ILOD_medium_70 <- IIln_ILOD_medium_params %>%
  filter(samp_size == 70)
try(execute_loop(25, IIln_ILOD_medium_70))

IIln_ILOD_medium_80 <- IIln_ILOD_medium_params %>%
  filter(samp_size == 80)
try(execute_loop(25, IIln_ILOD_medium_80))

IIln_ILOD_medium_90 <- IIln_ILOD_medium_params %>%
  filter(samp_size == 90)
try(execute_loop(25, IIln_ILOD_medium_90))


## --------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(123)
IIln_ILOD_large_params <- expand_grid(
  samp_size = seq(100, 500, 100),
  !!!param_args
) %>%
  mutate(p.mix = sample(p.mix_seq, size = n(), replace = T))%>%
  constrainGSD()

IIln_ILOD_large_params


## --------------------------------------------------------------------------------------------------------------------------------------------------------
IIln_ILOD_large_100_.3 <- IIln_ILOD_large_params %>%
  filter(samp_size == 100, p.mix == .3)
try(execute_loop(25, IIln_ILOD_large_100_.3))

IIln_ILOD_large_100_.4 <- IIln_ILOD_large_params %>%
  filter(samp_size == 100, p.mix == .4)
try(execute_loop(25, IIln_ILOD_large_100_.4))

IIln_ILOD_large_100_.5 <- IIln_ILOD_large_params %>%
  filter(samp_size == 100, p.mix == .5)
try(execute_loop(25, IIln_ILOD_large_100_.5))

IIln_ILOD_large_100_.6 <- IIln_ILOD_large_params %>%
  filter(samp_size == 100, p.mix == .6)
try(execute_loop(25, IIln_ILOD_large_100_.6))

IIln_ILOD_large_100_.7 <- IIln_ILOD_large_params %>%
  filter(samp_size == 100, p.mix == .7)
try(execute_loop(25, IIln_ILOD_large_100_.7))


IIln_ILOD_large_200_.3 <- IIln_ILOD_large_params %>%
  filter(samp_size == 200, p.mix == .3)
try(execute_loop(25, IIln_ILOD_large_200_.3))

IIln_ILOD_large_200_.4 <- IIln_ILOD_large_params %>%
  filter(samp_size == 200, p.mix == .4)
try(execute_loop(25, IIln_ILOD_large_200_.4))

IIln_ILOD_large_200_.5 <- IIln_ILOD_large_params %>%
  filter(samp_size == 200, p.mix == .5)
try(execute_loop(25, IIln_ILOD_large_200_.5))

IIln_ILOD_large_200_.6 <- IIln_ILOD_large_params %>%
  filter(samp_size == 200, p.mix == .6)
try(execute_loop(25, IIln_ILOD_large_200_.6))

IIln_ILOD_large_200_.7 <- IIln_ILOD_large_params %>%
  filter(samp_size == 200, p.mix == .7)
try(execute_loop(25, IIln_ILOD_large_200_.7))


IIln_ILOD_large_300_.3 <- IIln_ILOD_large_params %>%
  filter(samp_size == 300, p.mix == .3)
try(execute_loop(25, IIln_ILOD_large_300_.3))

IIln_ILOD_large_300_.4 <- IIln_ILOD_large_params %>%
  filter(samp_size == 300, p.mix == .4)
try(execute_loop(25, IIln_ILOD_large_300_.4))

IIln_ILOD_large_300_.5 <- IIln_ILOD_large_params %>%
  filter(samp_size == 300, p.mix == .5)
try(execute_loop(25, IIln_ILOD_large_300_.5))

IIln_ILOD_large_300_.6 <- IIln_ILOD_large_params %>%
  filter(samp_size == 300, p.mix == .6)
try(execute_loop(25, IIln_ILOD_large_300_.6))

IIln_ILOD_large_300_.7 <- IIln_ILOD_large_params %>%
  filter(samp_size == 300, p.mix == .7)
try(execute_loop(25, IIln_ILOD_large_300_.7))


IIln_ILOD_large_400_.3 <- IIln_ILOD_large_params %>%
  filter(samp_size == 400, p.mix == .3)
try(execute_loop(25, IIln_ILOD_large_400_.3))

IIln_ILOD_large_400_.4 <- IIln_ILOD_large_params %>%
  filter(samp_size == 400, p.mix == .4)
try(execute_loop(25, IIln_ILOD_large_400_.4))

IIln_ILOD_large_400_.5 <- IIln_ILOD_large_params %>%
  filter(samp_size == 400, p.mix == .5)
try(execute_loop(25, IIln_ILOD_large_400_.5))

IIln_ILOD_large_400_.6 <- IIln_ILOD_large_params %>%
  filter(samp_size == 400, p.mix == .6)
try(execute_loop(25, IIln_ILOD_large_400_.6))

IIln_ILOD_large_400_.7 <- IIln_ILOD_large_params %>%
  filter(samp_size == 400, p.mix == .7)
try(execute_loop(25, IIln_ILOD_large_400_.7))


IIln_ILOD_large_500_.3 <- IIln_ILOD_large_params %>%
  filter(samp_size == 500, p.mix == .3)
try(execute_loop(25, IIln_ILOD_large_500_.3))

IIln_ILOD_large_500_.4 <- IIln_ILOD_large_params %>%
  filter(samp_size == 500, p.mix == .4)
try(execute_loop(25, IIln_ILOD_large_500_.4))

IIln_ILOD_large_500_.5 <- IIln_ILOD_large_params %>%
  filter(samp_size == 500, p.mix == .5)
try(execute_loop(25, IIln_ILOD_large_500_.5))

IIln_ILOD_large_500_.6 <- IIln_ILOD_large_params %>%
  filter(samp_size == 500, p.mix == .6)
try(execute_loop(25, IIln_ILOD_large_500_.6))

IIln_ILOD_large_500_.7 <- IIln_ILOD_large_params %>%
  filter(samp_size == 500, p.mix == .7)
try(execute_loop(25, IIln_ILOD_large_500_.7))


