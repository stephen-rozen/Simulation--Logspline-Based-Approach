## ----setup, include=FALSE-------------------------------------------------------------------------------------------------------------------------------
#devtools::install_github("cran/NADA")
library(NADA)
library(logspline)
library(EnvStats)
library(cubature)
library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(tibble)


## ----calculate_summary_stats_on_uncensored_data---------------------------------------------------------------------------------------------------------
calculate_summary_stats <- function(some_data) {
  # will use statistics on the sim, instead of the param used to generate the distro, as the benchmark
  true_mean <- mean(some_data) 
  true_sd <- sd(some_data)
  true_gm <- exp(mean(log(some_data)))
  true_gsd <- exp(sd(log(some_data)))
  
  list(
    mean  = true_mean,
    sd    = true_sd, 
    gm = true_gm,
    gsd    = true_gsd
  )
  
}


## ----Complete_Case_aka_Uncensored_Error_Calculation-----------------------------------------------------------------------------------------------------
#error refers to "complete case error", the error relative to the statistic from the uncensored distribution
error_calculation <- function(lod_est_df, true_list){
  true_stat <- unlist(true_list)
  
  lod_est_df %>%
    mutate(true_stat = rep(unlist(true_list), length.out = n())) %>%
    mutate(complete_case_error = est-true_stat)
  # error compared to stat of full uncen data
  #could implement case_when matching
  
}


## ----Setup_Helpers--------------------------------------------------------------------------------------------------------------------------------------
saveRDS_named <- function(object, file = NULL, envir = parent.frame(), ...) {
  name <- deparse(substitute(object))
  if (is.null(file)) {
    file <- paste0(name, ".rds")
  }
  saveRDS(object, file = file, ...)
}

#create dir for checkpoints
create_dir <- function(outdir) {
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  } else {
    message("Directory already exists; skipping dir.create().")
  }
}


## ----Looping_helpers------------------------------------------------------------------------------------------------------------------------------------
# retry wrapper
run_rep_with_retry <- function(rep_id, params_df, outdir, maxtries = 6, verbose = TRUE) {
  for (a in seq_len(maxtries)) {
    ok <- tryCatch({
      compute_one_rep(rep_id, params_df, outdir)
      TRUE
    }, error = function(e) {
      if (verbose) {
        msg <- substr(conditionMessage(e), 1L, 200L)  # avoid huge prints
        message(sprintf("[%s] rep %03d attempt %d/%d failed: %s",
                        format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                        rep_id, a, maxtries, msg))
      }
      FALSE
    })
    if (ok) return(invisible(TRUE))
    if (a < maxtries) Sys.sleep(min(30, 2^(a - 1)))  # backoff, skip after last
  }
  stop(sprintf("rep %03d failed after %d attempts", rep_id, maxtries))
}


# one-button runner:
# call like: execute_loop(50, ln_IILOD_large_params)
execute_loop <- function(n_reps, params) {
  

  params_sym <- rlang::ensym(params)
  params_df  <- rlang::eval_tidy(params_sym)
  
  base_name  <- sub("_params$", "", rlang::as_string(params_sym))
  results_sym <- rlang::sym(paste0(base_name, "_results"))
  outdir     <- rlang::as_string(results_sym)    # checkpoint folder = results name
  
  create_dir(outdir)

  # cross-platform progress bar ----------------------------------------------
  make_pb <- function(total, title = "Monte Carlo Progress") {
    if (.Platform$OS.type == "windows" && interactive()) {
      pb <- utils::winProgressBar(title = title, label = "Starting…", min = 0, max = total, width = 400)
      attr(pb, "pb_kind") <- "win"
    } else {
      pb <- utils::txtProgressBar(min = 0, max = total, style = 3)
      attr(pb, "pb_kind") <- "txt"
    }
    pb
  }
  
  set_pb <- function(pb, value, label = NULL) {
    kind <- attr(pb, "pb_kind")
    if (identical(kind, "win")) {
      utils::setWinProgressBar(pb, value, label = label)
    } else if (identical(kind, "txt")) {
      utils::setTxtProgressBar(pb, value)
    }
  }
  
  close_pb <- function(pb) {
    if (!is.null(pb)) close(pb)
  }
  # ------
  pb <- make_pb(n_reps)
  on.exit(close_pb(pb), add = TRUE)
  

  repeat {
    done_count   <- length(list.files(outdir, pattern = "^rep_\\d+\\.rds$"))
    if (done_count >= n_reps) break

    remaining_ids <- seq.int(done_count + 1, n_reps)

    # update per rep
    for (rep_id in remaining_ids) {
      run_rep_with_retry(rep_id, params_df, outdir)  # saves rep_XXX.rds
      done_now <- length(list.files(outdir, pattern = "^rep_\\d+\\.rds$"))
      set_pb(pb, done_now, label = sprintf("Completed %d/%d", done_now, n_reps))
      # utils::flush.console()  # uncomment if using Rgui and it feels laggy
    }
  }
}



## ----Modified_Oldlogspline_function---------------------------------------------------------------------------------------------------------------------------
oldlogspline <- function (uncensored, right, left, interval, lbound, ubound, 
    nknots, knots, penalty, delete = TRUE) 
{
  #   annoying console printouts converted to messages, permits
  #   use of suppressMessages
    nsample <- rep(0, 6)
    if (!missing(uncensored)) 
        uncensored <- unstrip(uncensored)
    if (!missing(right)) 
        right <- unstrip(right)
    if (!missing(left)) 
        left <- unstrip(left)
    if (!missing(interval)) 
        interval <- unstrip(interval)
    if (!missing(knots)) 
        knots <- unstrip(knots)
    if (!missing(interval)) {
        if (length(interval[1, ]) != 2) 
            stop("interval must have two columns")
        if (min(abs(interval[, 1] - interval[, 2])) < 0) 
            stop("not all lower bounds smaller than upper bounds")
        nsample[3] <- length(interval)/2
        nsample[1] <- length(interval)/2
        if (!missing(lbound)) 
            interval[interval[, 1] < lbound, 1] <- lbound
        if (!missing(ubound)) 
            interval[interval[, 2] > ubound, 2] <- ubound
        sample <- as.vector(t(interval))
        ror <- order(interval[, 1], interval[, 2])
        if (nsample[3] > 1) {
            ro1 <- interval[ror[(1:(nsample[3] - 1))], 1] == 
                interval[ror[2:nsample[3]], 1]
            ro2 <- interval[ror[(1:(nsample[3] - 1))], 2] == 
                interval[ror[2:nsample[3]], 2]
            nsample[6] <- nsample[3] - sum(ro1 + ro2 == 2)
        }
        else nsample[6] <- 1
    }
    if (!missing(uncensored)) {
        uncensored2 <- uncensored[!is.na(uncensored)]
        u2 <- length(uncensored) - length(uncensored2)
        if (u2 > 0) 
            print(paste("***", u2, " NAs ignored in uncensored"))
        uncensored <- uncensored2
        if (nsample[1] > 0) 
            sample <- c(uncensored, sample)
        if (nsample[1] == 0) 
            sample <- uncensored
        nsample[1] <- length(uncensored) + nsample[1]
        nsample[2] <- length(uncensored)
        uncensored <- sort(uncensored)
        if (nsample[2] > 1) 
            nsample[6] <- sum(uncensored[2:nsample[2]] != uncensored[1:(nsample[2] - 
                1)]) + 1 + nsample[6]
        else nsample[6] <- nsample[6] + 1
    }
    if (nsample[1] == 0) 
        stop("you either need uncensored or interval censored data")
    if (!missing(right)) {
        if (nsample[1] > 0) 
            sample <- c(sample, right)
        if (nsample[1] == 0) 
            sample <- right
        nsample[1] <- length(right) + nsample[1]
        nsample[4] <- length(right)
        right <- sort(right)
        if (nsample[4] > 1) {
            nsample[6] <- sum(right[2:nsample[4]] != right[1:(nsample[4] - 
                1)]) + 1 + nsample[6]
        }
        else nsample[6] <- nsample[6] + 1
    }
    if (!missing(left)) {
        if (nsample[1] > 0) 
            sample <- c(sample, left)
        if (nsample[1] == 0) 
            sample <- left
        nsample[1] <- length(left) + nsample[1]
        nsample[5] <- length(left)
        left <- sort(left)
        if (nsample[5] > 1) {
            nsample[6] <- sum(left[2:nsample[5]] != left[1:(nsample[5] - 
                1)]) + 1 + nsample[6]
        }
        else nsample[6] <- nsample[6] + 1
    }
    if (missing(penalty)) 
        penalty <- log(nsample[1])
    n1 <- 4 * nsample[1]^0.2 + 1
    if (!missing(nknots)) 
        n1 <- nknots + 1
    if (!missing(knots)) 
        n1 <- length(knots) + 1
    if (!missing(knots)) {
        nknots <- length(knots)
        knots <- sort(knots)
        iautoknot <- 0
        if (knots[1] > min(sample)) 
            stop("first knot must be smaller than smallest sample")
        if (knots[nknots] < max(sample)) 
            stop("last knot should be larger than largest sample")
    }
    else {
        if (missing(nknots)) 
            nknots <- 0
        knots <- vector(mode = "double", length = max(nknots, 
            50))
        iautoknot <- 1
    }
    xbound <- c(1, 0, 0, 0, 0)
    if (!missing(lbound)) {
        xbound[2] <- 1
        xbound[3] <- lbound
        if (lbound > min(sample)) 
            stop("lbound should be smaller than smallest sample")
    }
    if (!missing(ubound)) {
        xbound[4] <- 1
        xbound[5] <- ubound
        if (ubound < max(sample)) 
            stop("ubound should be larger than largest sample")
    }
    SorC <- vector(mode = "integer", length = 35)
    SorC[1] <- 1
    SorC[17] <- 0
    nsample[6] <- nsample[6] - 1
    if (length(table(sample)) < 3) 
        stop("Not enough unique values")
    z <- .C("logcensor", as.integer(c(delete, 0, 0, 0, 0)), as.integer(c(iautoknot, 
        0, 0, 0, 0)), as.double(c(sample, 0, 0, 0, 0)), as.integer(c(nsample, 
        0, 0, 0, 0)), bd = as.double(c(xbound, 0, 0, 0, 0)), 
        SorC = as.integer(c(SorC, 0, 0, 0, 0)), nk = as.integer(nknots), 
        kt = as.double(c(knots, 0, 0, 0, 0)), cf = as.double(c(knots, 
            0, 0, 0, 0)), as.double(c(penalty, 0, 0, 0, 0)), 
        as.double(c(sample, 0, 0, 0, 0)), as.double(c(sample, 
            0, 0, 0, 0)), logl = as.double(rep(0, n1 + 1 + 10)), 
        PACKAGE = "logspline")
    SorC <- z$SorC
    if (SorC[1] == -1 && SorC[28] == 0 && nsample[1] != nsample[2] && 
        nsample[2] > 15) {
        SorC <- vector(mode = "integer", length = 35)
        SorC[1] <- 1
        SorC[17] <- 1
        z <- .C("logcensor", as.integer(c(delete, 0, 0, 0, 0)), 
            as.integer(c(iautoknot, 0, 0, 0, 0)), as.double(c(sample, 
                0, 0, 0, 0)), as.integer(c(nsample, 0, 0, 0, 
                0)), bd = as.double(c(xbound, 0, 0, 0, 0)), SorC = as.integer(c(SorC, 
                0, 0, 0, 0)), nk = as.integer(nknots), kt = as.double(c(knots, 
                0, 0, 0, 0)), cf = as.double(c(knots, 0, 0, 0, 
                0)), as.double(c(penalty, 0, 0, 0, 0)), as.double(c(sample, 
                0, 0, 0, 0)), as.double(c(sample, 0, 0, 0, 0)), 
            logl = as.double(rep(0, n1 + 1 + 10)), PACKAGE = "logspline")
    }
    bound <- c(z$bd[2], z$bd[3], z$bd[4], z$bd[5])
    SorC <- z$SorC
    if (abs(SorC[1]) > 2) {
        for (i in 3:abs(SorC[1])) message(paste("===> warning: knot ", 
            SorC[i - 1], " removed - double knot"))
        if (SorC[1] < 0) 
            SorC[1] <- -1
        if (SorC[1] == 23) 
            SorC[1] <- -3
    }
    if (abs(SorC[1]) > 3) {
        cat("* several double knots suggests that your data is *\n")
        cat("* strongly rounded, attention might be required   *\n")
        SorC[1] <- 1
    }
    if (SorC[1] == -3) 
        stop("* too many double knots")
    if (SorC[1] == -1 && SorC[28] == 0) 
        stop("* no convergence")
    if (SorC[28] > 0) 
        message(paste("* convergence problems, smallest number of knots", 
            " tried is ", SorC[28] + 1, " *\n"))
    if (SorC[1] == 2) 
        stop("* sample is too small")
    if (SorC[1] == -2) 
        stop(paste("* too many knots, at most ", SorC[2], "knots possible"))
    if (SorC[22] == 1) {
        cat("possible discontinuity at lower end\n")
        cat(paste("consider rerunning with lbound=", z$kt[1], 
            "\n"))
    }
    if (SorC[22] == 3) {
        message("possible infinite density at lower end")
        message("running program with fewer knots")
    }
    if (SorC[21] == 1) 
        cat("running with maximum degrees of freedom\n")
    if (SorC[25] > 0) 
        cat("* problems are possibly due to a very heavy right tail *\n")
    if (SorC[24] > 0) 
        cat("* problems are possibly due to a very heavy left tail *\n")
    if (SorC[23] == 3) {
        cat("possible infinite density at upper end\n")
        cat("running program with fewer knots\n")
    }
    if (SorC[23] == 1) {
        cat("possible discontinuity at upper end\n")
        cat(paste("consider rerunning with ubound=", z$kt[z$nk], 
            "\n"))
    }
    if (delete && SorC[28] > 0) 
        delete <- 3
    coef <- z$cf[1:(z$nk + 2)]
    uu <- 3:z$nk
    if (delete == FALSE) 
        uu <- 1
    fit <- list(coef = coef, knots = z$kt[1:z$nk], bound = bound, 
        logl = z$logl[uu], penalty = penalty, sample = nsample[1], 
        delete = delete)
    class(fit) <- "oldlogspline"
    fit
}


## ----LOD_handling_methods-------------------------------------------------------------------------------------------------------------------------------
estimate <- function(object){
  
  list(
    "mean"  = mean(object),
    "sd"    = sd(object),
    "gm" = exp(mean(log(object))),
    "gsd"   = exp(sd(log(object)))
  )
  
}

ES_ROS_estimate <- function(censored_sim, censored_bool){
  fit<- elnormCensored(censored_sim, censored_bool, method = "impute.w.qq.reg")
  fit_alt <- elnormAltCensored(censored_sim, censored_bool, method = "impute.w.qq.reg")
  
  meanlog <- unname(fit$parameters["meanlog"])
  sdlog <- unname(fit$parameters["sdlog"])
  mean <- unname(fit_alt$parameters["mean"])
  sd <- mean * unname(fit_alt$parameters["cv"])

  list(
    "mean"  = mean,
    "sd"    = sd,
    "gm"    = exp(meanlog),
    "gsd"   = exp(sdlog)
  )

}

ROS_estimate <- function(censored_sim, censored_bool){
  myros <- ros(obs = censored_sim, censored = censored_bool)
  
  ros_fit <- as.data.frame(myros)$modeled
  
  estimate(ros_fit)
}

MLE_estimate <- function(censored_sim, censored_bool, method){
  #accepts methods "rmle", "bcmle", "mle"
  fit<- elnormCensored(censored_sim, censored_bool, "mle")
  meanlog <- unname(fit$parameters["meanlog"])
  sdlog <- unname(fit$parameters["sdlog"])   
  
  if(method == 'mle'){
    return(
      list(
    "mean"  = exp(meanlog + sdlog^2/2),
    "sd"    = sqrt((exp(sdlog^2) - 1) * exp(2*meanlog + sdlog^2)),
    "gm"    = exp(meanlog),
    "gsd"   = exp(sdlog)
  )
    )
  }
  
  
  else if(method == 'bcmle'){
    # ingredients for ψ(g)
  n <- length(censored_sim)
  g <- sdlog^2 / 2
  
  # ψ(g) with K = 4 (so total terms = 1 + 4 = 5)
  psi5 <- {
    res <- 1
    for(k in 1:4) {
      num      <- (n - 1)^(2*k - 1) * g^k
      den_prod <- if(k > 1) prod(n + seq(1, by = 2, length.out = k - 1)) else 1
      den      <- n^k * den_prod * factorial(k)
      res      <- res + num / den
    }
    res
  }
  
  # bias‐corrected arithmetic mean
  mean_bc <- exp(meanlog) * psi5
  
  # plug‐in arithmetic SD, plus geometric mean & SD
  sd_bc <- sqrt((exp(sdlog^2) - 1) * exp(2*meanlog + sdlog^2))
  gm    <- exp(meanlog)
  gsd   <- exp(sdlog)
  
  list(
    mean = mean_bc,
    sd   = sd_bc,
    gm   = gm,
    gsd  = gsd
  )
  }
  
  else if(method == 'rmle'){
    # compute plotting positions
  pp <- hc.ppoints(censored_sim, censored_bool)

  # impute each censored at its fitted quantile
  imputed_values <- qlnorm(pp[censored_bool], 
                            meanlog, sdlog)

  # assemble full vector
  all_vals <- censored_sim
  all_vals[censored_bool] <- imputed_values

  return( estimate(all_vals) )
  }
  
  
}



calc_cenfit_negatives <- function(cf) {
  #helper for KM method
  if (!inherits(cf, "cenfit"))
    stop("Input must be a cenfit object.")
  s <- cf@survfit
  if (!is.null(s$strata))
    stop("Stratified cenfit objects not supported.")

  # 1) Find minimum event time
  ev_times <- s$time[s$n.event > 0]
  if (length(ev_times) == 0)
    stop("No detected events to estimate from.")
  min_time <- min(ev_times, na.rm = TRUE)

  # 2) Compute shift (only if min_time < 0)
  shift_amt <- if (min_time < 0) -min_time else 0

  # 3) Shift times if needed
  cf2 <- if (shift_amt > 0) {
    s2 <- s
    s2$time <- s2$time + shift_amt
    new("cenfit", survfit = s2)
  } else {
    cf
  }

  # 4) Extract mean and sd
  m    <- mean(cf2)["mean"]
  sd_v <- sd(cf2)

  # 5) Shift mean back (variance unaffected)
  if (shift_amt > 0) {
    m <- m - shift_amt
  }

  list(
    mean = unname(m),
    sd   = unname(sd_v)
  )
}






KM_estimate <- function(censored_sim, censored_bool, ...) {
  #cenmle and cenfit can take groups (like different labs or collection sites)
  
  # if any non‐positive values, abort and return NAs
  if (any(censored_sim <= 0, na.rm = TRUE)) {
    return(list(
      mean = NA_real_,
      sd   = NA_real_,
      gm   = NA_real_,
      gsd  = NA_real_
    ))
  }

  # fit KM to the raw data
  mycenfit  <- cenfit(obs = censored_sim, censored = censored_bool, ...)
  
  # fit KM to the log‐transformed data
  logcenfit <- cenfit(obs = log(censored_sim), censored = censored_bool, ...)
  accurate_log_vals <- calc_cenfit_negatives(logcenfit)
  
  
  list(
    mean = unname(mean(mycenfit)[1]),
    sd   = sd(mycenfit),
    gm   = exp((accurate_log_vals)[[1]]),
    gsd  = exp((accurate_log_vals)[[2]])
  )
}

substitution_estimate <- function(censored_sim, censored_bool, method){
  # methods: root2, half, LOD, zero
  # returns substitution object
  # only works if censored vector has LOD at censored indices (our sim does)
  
  # grab the LOD‐values
  subs_vals <- censored_sim[censored_bool]
  
  # 
  censored_sim[censored_bool] <- switch(
    method,
    root2 = subs_vals / sqrt(2),
    half  = subs_vals / 2,
    LOD   = subs_vals,
    zero  = 0
  )
  #print(censored_sim)
  estimate(censored_sim)
}

spline_fit <- function(censored_sim, censored_bool){
  # attempt logspline fit
  
  
  # fitting with lbound causes instability
  # Kooperberg and Stone note: "The user should be wary about extrapolation 
  # beyond the range of the data. In particular when all observations beyond a 
  # certain point are censored, as in Type I or Type II censoring, the 
  # reliability of conclusions about the right tail of the density may be 
  # severely limited.)" ... However, using the min of the data seemed to create 
  # instability in fitting oldlogspline
  
  lbound <- 0
  ubound <- max(censored_sim, na.rm = TRUE)

  fit <- tryCatch(
    suppressMessages({
      oldlogspline(
                 uncensored = censored_sim[!censored_bool],
                 left       = censored_sim[censored_bool],
                 lbound     = lbound,
                 ubound = ubound
                 )
      }),
    error = function(e) {
      message("Logspline fit failed: ", e$message)
      NULL
    }
  )

  return(fit)
}



spline_estimate <- function(censored_sim, censored_bool, fit) {
  # require vector length <14
  if (length(censored_sim) < 14) {
    return(list(
      mean = NA_real_,
      sd   = NA_real_,
      gm   = NA_real_,
      gsd  = NA_real_
    ))
  }

  # no zeros/negatives
  if (any(censored_sim <= 0, na.rm = TRUE)) {
    warning("Non-positive values detected; returning NAs for all estimates.")
    return(list(
      mean = NA_real_,
      sd   = NA_real_,
      gm   = NA_real_,
      gsd  = NA_real_
    ))
  }

  if (is.null(fit)) {
    return(list(
      mean = NA_real_,
      sd   = NA_real_,
      gm   = NA_real_,
      gsd  = NA_real_
    ))
  }

  lower_bound <- 0
  upper_bound <- max(censored_sim, na.rm = TRUE)
  samp_first_moment <- mean(censored_sim)
  samp_second_moment <- mean(censored_sim^2)

  # helper for fitted density
  pdf_x <- function(x) doldlogspline(x, fit)
  
  # vectorized integrand
  f_vec <- function(x) {
    if (is.matrix(x)) {
      xs    <- as.numeric(x[1, ])    # extract the single dimension
      n_pts <- ncol(x)
    } else {
      xs    <- as.numeric(x)
      n_pts <- length(xs)
    }
    vals <- xs * pdf_x(xs)           # length‑N numeric
    # return as 1×N matrix so cubature knows its vectorized
    matrix(vals, nrow = 1, ncol = n_pts)
  }
  
  # compute first moment
  first_moment <- tryCatch({
  # 1) run hcubature()
  res1 <- hcubature(
    f               = f_vec,
    lowerLimit      = lower_bound,
    upperLimit      = upper_bound,
    vectorInterface = TRUE,
    maxEval         = 1e4,
    tol             = 1e-6,
    absError        = 1e-8
  )
  
  # extract the integral
  integ1 <- res1$integral
  
  # threshold‐check against samp_first_moment
  if (!is.finite(integ1) || abs(integ1 - samp_first_moment) > 1e6) {
    warning("Integral invalid or differs by > 1e6; returning NA")
    NA_real_
  } else {
    integ1
  }
  

  
}, error = function(e) {
  message("Error computing first moment: ", e$message)
  NA_real_
})

  
    # vectorized integrand for second moment (x^2 * pdf(x))
  g_vec <- function(x) {
    # handle either a 1×N matrix or a plain vector
    if (is.matrix(x)) {
      xs    <- as.numeric(x[1, ])    # extract the single dimension
      n_pts <- ncol(x)
    } else {
      xs    <- as.numeric(x)
      n_pts <- length(xs)
    }
    vals <- xs^2 * pdf_x(xs)         # compute x^2 * pdf(x)
    matrix(vals, nrow = 1, ncol = n_pts)
  }
  
  # compute second moment
  second_moment <- tryCatch({
  # 1) run hcubature()
  res2 <- hcubature(
    f               = g_vec,
    lowerLimit      = lower_bound,
    upperLimit      = upper_bound,
    vectorInterface = TRUE,
    maxEval         = 1e4,
    tol             = 1e-6,
    absError        = 1e-8
  )
  
  # extract the integral
  integ <- res2$integral
  
  # threshold‐check
  if (!is.finite(integ) || abs(integ - samp_second_moment) > 1e6) {
    warning("Integral invalid or differs by > 1e6; returning NA")
    NA_real_
  } else {
    integ
  }
  
}, error = function(e) {
  message("Error computing second moment: ", e$message)
  NA_real_
})


  var_orig <- second_moment - first_moment^2

     # ---------- log-scale moments via hcubature (vectorized) ----------
  log_upper <- log(upper_bound)
  # finite surrogate for -Inf; exp(-700) ~ 9.86e-305 (safe, ≈0)
  log_lower <- -700

  # sample log moments (for safegaurds)
  samp_log_first_moment <- mean(log(censored_sim))
  samp_log_second_central_moment <- mean( (log(censored_sim) - samp_log_first_moment)^2 )


  # vectorized integrand for E[log X] = integrate(t f(e^t) e^t dt)
  l1_vec <- function(tmat) {
    if (is.matrix(tmat)) {
      ts   <- as.numeric(tmat[1, ])
      npts <- ncol(tmat)
    } else {
      ts   <- as.numeric(tmat)
      npts <- length(ts)
    }
    xs   <- exp(ts)                 # x = e^t
    dens <- numeric(npts)
    # avoid evaluating pdf at exactly 0
    keep <- xs > 0
    dens[keep] <- pdf_x(xs[keep])
    vals <- ts * dens * xs          # t f(e^t) e^t
    matrix(vals, nrow = 1, ncol = npts)
  }

  mu_log <- tryCatch({
    res_log1 <- hcubature(
      f               = l1_vec,
      lowerLimit      = log_lower,
      upperLimit      = log_upper,
      vectorInterface = TRUE,
      maxEval         = 1e4,
      tol             = 1e-6,
      absError        = 1e-8
    )
    integ_log1 <- res_log1$integral
    if (!is.finite(integ_log1) || abs(integ_log1 - samp_log_first_moment) > log(1e6)) {
      warning("Log first moment invalid or differs by > 1e6; returning NA")
      NA_real_
    } else {
      integ_log1
    }
  }, error = function(e) {
    message("Error computing mu_log (hcubature): ", e$message)
    NA_real_
  })

  # vectorized integrand for Var[log X] = integrate((t - mu_log)^2 f(e^t) e^t dt)
  l2c_vec <- function(tmat) {
    if (is.matrix(tmat)) {
      ts   <- as.numeric(tmat[1, ])
      npts <- ncol(tmat)
    } else {
      ts   <- as.numeric(tmat)
      npts <- length(ts)
    }
    xs   <- exp(ts)
    dens <- numeric(npts)
    keep <- xs > 0
    dens[keep] <- pdf_x(xs[keep])
    vals <- (ts - mu_log)^2 * dens * xs
    matrix(vals, nrow = 1, ncol = npts)
  }

  var_log <- tryCatch({
    res_log2 <- hcubature(
      f               = l2c_vec,
      lowerLimit      = log_lower,
      upperLimit      = log_upper,
      vectorInterface = TRUE,
      maxEval         = 1e4,
      tol             = 1e-6,
      absError        = 1e-8
    )
    integ_log2 <- res_log2$integral
    if (!is.finite(integ_log2) || abs(integ_log2 - samp_log_second_central_moment) > log(1e6)) {
      warning("Log variance invalid or differs by > 1e6; returning NA")
      NA_real_
    } else {
      integ_log2
    }
  }, error = function(e) {
    message("Error computing var_log (hcubature): ", e$message)
    NA_real_
  })

  # back-transform
  gm  <- if (is.na(mu_log)) NA_real_ else exp(mu_log)
  gsd <- if (is.na(var_log)) NA_real_ else exp(sqrt(var_log))


  list(
    mean = first_moment,
    sd   = sqrt(var_orig),
    gm   = gm,
    gsd  = gsd
  )
}



robust_spline <- function(censored_sim, censored_bool, fit) {
  # censored_sim: numeric vector of observations (with LODs at censored positions)
  # censored_bool: logical vector (TRUE = censored, FALSE = detected)
  if (length(censored_sim) != length(censored_bool)) {
    stop("'censored_sim' and 'censored_bool' must be the same length")
  }
  
  if (length(censored_sim) < 14) {
    # message("Need at least 14 uncensored points for logspline; returning NAs.")
    return(list(
      mean = NA_real_,
      sd   = NA_real_,
      gm   = NA_real_,
      gsd  = NA_real_
    ))
  }

   if (is.null(fit)) {
    return(list(
      mean = NA_real_,
      sd   = NA_real_,
      gm   = NA_real_,
      gsd  = NA_real_
    ))
  }

  # compute ROS plotting positions
  pp <- hc.ppoints(censored_sim, censored_bool)

  # impute each censored at its fitted quantile
  imputed_values <- qoldlogspline(pp[censored_bool], fit)

  # assemble full vector
  all_vals <- censored_sim
  all_vals[censored_bool] <- imputed_values

  estimate(all_vals)
}





beta_sub <- function(censored_sim, censored_bool, LOD){
  # LOD can be scalar or vector
  # only works if censored vector has LOD at censored indices, our sim does
  
  #calculate an average field LOD according to Hewett and Gasner
  if(length(LOD)> 1){
    average_limit_df<- tibble(dl = LOD) %>%
      rowwise() %>%
      mutate(mi = sum(censored_sim[censored_bool] == dl)) %>%
      mutate(summand = mi * log(dl)) %>%
      ungroup()
  
    LOD <- exp(
      1/sum(average_limit_df$mi) * sum(average_limit_df$summand)
      )
  }

  
  
  n <- length(censored_sim)
  k <- length(censored_sim[censored_bool == T])
  
  # handle edge cases k=1 or k=0
  if (k == 0) {
    # no censored values: use direct estimates
    mean_est <- mean(censored_sim)
    sd_est   <- sd(censored_sim)
    gm_est   <- exp(mean(log(censored_sim)))
    gsd_est  <- exp(sd(log(censored_sim)))
    return(list(
      mean = mean_est,
      sd   = sd_est,
      gm   = gm_est,
      gsd  = gsd_est
    ))
  }
  if (k == n) {
    # all values censored: cannot estimate
    return(list(
      mean = NA_real_,
      sd   = NA_real_,
      gm   = NA_real_,
      gsd  = NA_real_
    ))
  }
  
  ln_detects <- log(censored_sim[censored_bool == F])
  if (any(censored_sim <= 0)) stop("Values must be positive for log-transformation")
  y_hat <- (1 / (n-k)) * sum(ln_detects)
  z <- qnorm((k/n))
  f_z <- dnorm(z, 0,1) / (1-pnorm(z, 0, 1))
  s_hat_y <- (y_hat - log(LOD))/ (f_z - z)
  f_sy_z <- (1 - pnorm(z - (s_hat_y/n), 0 ,1)) /
        (1 - pnorm(z, 0, 1))

  beta_mean <- (n/k) * pnorm(z-s_hat_y, 0 ,1) *
    exp(-s_hat_y * z + (s_hat_y**2)/2)
  sim_tmp <- censored_sim ## avoid overwriting censored_sim
  sim_tmp[censored_bool == T] <- beta_mean*LOD
  mean_est <- mean(sim_tmp)

  beta_GM <- exp(-(n-k)*n/k*log(f_sy_z)-s_hat_y*z-(n-k)/(2*k*n)*s_hat_y**2)
  sim_tmp <- censored_sim
  sim_tmp[censored_bool == T] <- beta_GM*LOD
  gm_est <- exp(mean(log(sim_tmp)))

  ratio <- mean_est / gm_est

if (is.na(ratio) || gm_est <= 0) {
  # if can't compute
  s_y <- NA_real_
} else if (ratio <= 1) {
  ## Hewett and ganser note this in methods
  s_y <- 0
} else {
  s_y <- sqrt(2 * n/(n - 1) * log(ratio))
}
  gsd_est <- exp(s_y)

#   q05 <- exp(
#     log(gm_est) - s_y**2/(2*n) + qnorm(.05)*s_y
#   )

  list(
    "mean" = mean_est,
    "sd" = s_y,
    "gm" = gm_est,
    "gsd" = gsd_est
  )
}

wlod2_estimate <- function(censored_sim, censored_bool, LOD) {
  
  #calculate an average field LOD according to Hewett and Gasner
  if(length(LOD)> 1){
    average_limit_df<- tibble(dl = LOD) %>%
      rowwise() %>%
      mutate(mi = sum(censored_sim[censored_bool] == dl)) %>%
      mutate(summand = mi * log(dl)) %>%
      ungroup()
  
    LOD <- exp(
      1/sum(average_limit_df$mi) * sum(average_limit_df$summand)
      )
  }
  LOD <- unname(LOD)
  
  # ——— Helper: omega under lognormal assumption ———
  calc_omega <- function(sd_obs, mean_obs, pct_below) {
    # ω = a · t^b,    t = sd(X≥LOD) / mean(X≥LOD)
    t_ratio <- sd_obs / mean_obs
    a_const <- 1e-3 * (1482  - 9.166  * pct_below)
    b_const <- -1e-3 * (166   + 2.569  * pct_below)
    a_const * t_ratio^b_const
  }

  # ——— Split data  ———
  n_total     <- length(censored_sim)
  observed    <- censored_sim[!censored_bool]
  n_obs       <- length(observed)
  n_cens      <- sum(censored_bool)

  # ——— Compute omega and substitute ———
  mean_obs    <- mean(observed)
  sd_obs      <- sd(observed)
  pct_below   <- n_cens / n_total * 100
  omega       <- calc_omega(sd_obs, mean_obs, pct_below)

  substituted <- omega * (LOD / 2)

  # ——— Summary statistics per the paper ———
  # Arithmetic mean
  am <- (n_cens * substituted + n_obs * mean_obs) / n_total

  # Geometric mean
  mean_log_obs <- mean(log(observed))
  gm <- exp((n_cens * log(substituted) + n_obs * mean_log_obs) / n_total)

  # Pooled sample SD (denominator n–1)
  sum_sq_cens <- n_cens * (substituted - am)^2
  sum_sq_obs  <- sum((observed - am)^2)
  sd <- sqrt((sum_sq_cens + sum_sq_obs) / (n_total - 1))


  list(
    "mean" = am,
    "sd"   = sd,
    "gm"   = gm,
    "gsd"  = NA_real_
  )
}



## -------------------------------------------------------------------------------------------------------------------------------------------------------
nooverwrite <- function(object, path, name, ext = "rds", sep = "_", ...) {
  # ensure the directory exists
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  
  # build the base filepath (without numeric suffix)
  base <- file.path(path, name)
  file <- paste0(base, ".", ext)
  
  # if it exists, find the next available numbered filename
  if (file.exists(file)) {
    i <- 1L
    repeat {
      candidate <- paste0(base, sep, i, ".", ext)
      if (!file.exists(candidate)) {
        file <- candidate
        break
      }
      i <- i + 1L
    }
  }
  
  # save without notification
  saveRDS(object, file = file, ...)
}


