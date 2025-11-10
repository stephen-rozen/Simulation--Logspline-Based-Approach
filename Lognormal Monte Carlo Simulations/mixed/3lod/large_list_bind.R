library(purrr)
library(dplyr)


# list dirs with the right pattern
list_large_dirs <- function(path = ".") {
  dirs <- list.dirs(path, full.names = TRUE, recursive = FALSE)
  grep("IIln_IIILOD_large_[0-9]+_\\.[0-9]+$", dirs, value = TRUE)
}
getwd()
list_large_dirs()
# list files in each dir
list_files_in_dirs <- function(path = ".") {
  dirs <- list_large_dirs(path)
  files <- map(dirs, function(x) sort(list.files(x, full.names = TRUE)))
  names(files) <- basename(dirs)
  files
}

# align by index and bind across dirs
bind_files_by_index <- function(path = ".") {
  files <- list_files_in_dirs(path)
  n_files <- length(files[[1]])   # same length in every dir
  
  map(seq_len(n_files), function(i) {
    ith_files <- map(files, function(x) x[[i]])
    dfs <- map(ith_files, function(x) readRDS(x))
    bind_rows(dfs, .id = "source_dir")
  })
}

combined_list <- bind_files_by_index("path/to/results")
results <- purrr::map(combined_list, readRDS) %>% purrr::list_rbind()
saveRDS(results, "IIln_IIILOD_large_results.rds")