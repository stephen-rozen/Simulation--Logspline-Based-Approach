# fix_path <- function(x) normalizePath(x, winslash = "/")
# # Example:
# 
# library(callr)
# run_clean <- function(path, echo = TRUE) {
#   callr::r(
#     function(p, echo) {
#       source(p, chdir = TRUE, echo = echo, print.eval = echo)
#     },
#     args = list(normalizePath(path), echo),
#     show = TRUE
#   )
# }

# ----------------------------------------------------------
# run_clean <- function(p){ source(p, chdir = TRUE, echo = T, print.eval = T) }

library(callr)
run_clean <- function(path, echo = TRUE) {
  callr::r(
    function(p, echo) {
      source(p, chdir = TRUE, echo = echo, print.eval = echo)
    },
    args = list(normalizePath(path), echo),
    show = TRUE
  )
}
# Uses the Rproj root at: LOD_method_assessment/

if (!requireNamespace("rprojroot", quietly = TRUE)) install.packages("rprojroot")

# Find the project root (dir that contains *.Rproj; falls back to git root if needed)
root <- rprojroot::find_root(
  rprojroot::has_file_pattern("\\.Rproj$") | rprojroot::is_git_root
)

# Base dir where this script's MC files live
mc_dir <- file.path(root, "Lognormal Monte Carlo Simulations")

# Relative script paths (from mc_dir).
scripts_rel <- c(
  file.path("mixed", "3lod", "mix_3_lod_mc_lin.R"),
  file.path("single distribution", "2lod", "2_lod_mc.R"),
  file.path("single distribution", "3lod", "3_lod_mc.R"),
  file.path("mixed", "1lod", "mix_1_lod_mc.R"),
  file.path("mixed", "2lod", "mix_2_lod_mc.R"),
  file.path("oscillating", "2lod", "sin_2_lod_mc.R"),
  file.path("oscillating", "3lod", "sin_3_lod_mc.R"),
  file.path("oscillating", "1lod", "sin_1_lod_mc.R"),
  #etc... run all the files
)

# Absolute, normalized paths (portable across Windows/Linux)
scripts <- normalizePath(file.path(mc_dir, scripts_rel), winslash = "/", mustWork = FALSE)

# Run them (continue on error)
for (p in scripts) {
  if (!file.exists(p)) {
    message("Not found: ", p)
    next
  }
  message("Running: ", p)
  try(run_clean(p), silent = FALSE)
}



