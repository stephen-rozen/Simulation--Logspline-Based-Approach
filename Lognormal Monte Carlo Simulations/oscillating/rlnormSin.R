library(nleqslv)


rlnorm_sin_fast <- function(
    n, meanlog = 0, sdlog = 1, period_hours = 24,
    calib_n = min(n, 2048),           # smaller grid for solving
    include_noise = TRUE, noise_sdlog = 0.10,
    seed = NULL, max_tries = 20,
    ftol = 1e-7, xtol = 1e-7, maxit = 200, eps = 1e-8
){
  if (!is.null(seed)) set.seed(seed)

  # --- Calibration grid (small) ---
  t_cal <- seq(0, 48, length.out = calib_n)
  s_cal <- sin(2*pi*t_cal/period_hours)
  # noise_cal <- if (include_noise) stats::rlnorm(calib_n, 0, noise_sdlog) else rep(1, calib_n) ##remove for better stability https://chatgpt.com/s/t_68a84ae9fa888191a09e28cfa9f86401
  noise_cal <- rep(1, calib_n)
  
  # positivity reparam: base > amp > 0
  unpack <- function(theta){
    amp  <- exp(theta[2])
    base <- exp(theta[1]) + amp + eps
    c(base = base, amp = amp)
  }

  resid_fun <- function(theta){
    p <- unpack(theta)
    core <- p["base"] + p["amp"] * s_cal
    if (!all(is.finite(core)) || any(core <= 0)) return(c(1e6, 1e6))
    lx <- log(core * noise_cal)  # deterministic within the solve
    c(mean(lx) - meanlog, stats::sd(lx) - sdlog)
  }

  # --- Multistart (no recursion) ---
  sol <- NULL
  for(i in seq_len(max_tries)){
    x0 <- if (i == 1) c(0,0) else stats::rnorm(2, 0, 1)
    try_sol <- try(
      nleqslv::nleqslv(x = x0, fn = resid_fun,
                       control = list(ftol = ftol, xtol = xtol, maxit = maxit)),
      silent = TRUE
    )
    if (!inherits(try_sol, "try-error") && try_sol$termcd == 1) { sol <- try_sol; break }
  }
  if (is.null(sol) || sol$termcd != 1) stop("Failed to converge after ", max_tries, " attempts.")

  # --- Final long series (single allocation) ---
  p <- unpack(sol$x)
  t <- seq(0, 48, length.out = n)
  core <- p["base"] + p["amp"] * sin(2*pi*t/period_hours)
  if (include_noise) core <- core * stats::rlnorm(n, 0, noise_sdlog)
  core
}

# Try + fallback
safe_rlnorm_sin <- function(
    n, meanlog = 0, sdlog = 1, ...,
    .retry_ftol = 1e-6, .retry_xtol = 1e-6, .retry_max_tries = 50,
    .on_fail = c("NA", "error", "fallback_lnorm")
){
  .on_fail <- match.arg(.on_fail)
  
  bad <- function(x) {
    inherits(x, "try-error") || !is.numeric(x) || length(x) != n || any(!is.finite(x))
  }
  
  # 1) fast attempt (your current defaults)
  ok1 <- try(rlnorm_sin_fast(n = n, meanlog = meanlog, sdlog = sdlog, ...), silent = TRUE)
  if (!bad(ok1)) return(ok1)
  
  # 2) laxer retry (more iterations, looser tolerances / more starts)
  ok2 <- try(
    rlnorm_sin_fast(
      n = n, meanlog = meanlog, sdlog = sdlog,
      max_tries = .retry_max_tries, ftol = .retry_ftol, xtol = .retry_xtol, ...),
    silent = TRUE
  )
  if (!bad(ok2)) return(ok2)
  
  # 3) final fallback policy
  if (.on_fail == "fallback_lnorm") {
    return(stats::rlnorm(n = n, meanlog = meanlog, sdlog = sdlog))
  } else if (.on_fail == "error") {
    stop("sin-modulated generator failed after fast + retry attempts")
  } else { # "NA"
    warning("sin-modulated generator failed; returning NA vector")
    return(rep(NA_real_, n))
  }
}
