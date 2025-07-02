suppressPackageStartupMessages({
  library(parallel)
  library(benchmarkme)
  library(pryr)
  library(processx)
  library(hdf5r)
})

.this_file_path <- NULL
try({.this_file_path <- dirname(sys.frame(1)$ofile)}, silent = TRUE)
if (is.null(.this_file_path)) .this_file_path <- "."
.MONITOR_MEM_SH_PATH <- normalizePath(
  file.path(.this_file_path, "monitor_mem.sh"),
  mustWork = FALSE
)

TimerMemoryCollection = function(silent = TRUE) {
  env = environment()
  env$timings = list()
  env$delay = 0.15
  pid = Sys.getpid()
  with_timer = function(message, expr) {
    start = Sys.time()
    if (!silent) {
      cat(paste0(message, '...\n'))
    }
    result = NULL
    aborted = FALSE
    curr_process <- process$new(
        command = .MONITOR_MEM_SH_PATH,
        args = c("-p", as.character(pid)),
        stdout = "|"
      )
    process_pid <- curr_process$get_pid()
    Sys.sleep(env$delay)
    tryCatch({      
      result = invisible(eval(substitute(expr), parent.frame()))
    }, error = function(e) {
      aborted = TRUE
      stop(e)
    }, finally = {
      duration = as.numeric(difftime(Sys.time(), start, units = 'secs'))
      processx::run("kill", args = c(as.character(process_pid)))
      stdout_output <- curr_process$read_all_output()
      con <- textConnection(stdout_output) 
      df <- read.csv(con, header = FALSE, sep = ",", strip.white = TRUE, 
        col.names = c("Integer", "Percentage"))
      close(con)
      peak_mem <- max(df$Integer, na.rm = TRUE)
      percent <- max(df$Percentage, na.rm = TRUE)
      env$timings[[message]] = list(
        duration = duration,
        max_mem= peak_mem/1024/1024,
        mem_percent = percent,
        aborted = aborted
      )
      if (!silent) {
        time_str = format_time(duration)
        status = if (aborted) 'aborted after' else 'took'
        cat(sprintf('%s %s %s\n\n', message, status, time_str))
      }
    })
    gc()
    return(invisible(result))
  }
  
  format_time = function(duration, unit = NULL) {
    duration = as.numeric(duration)
    
    if (!is.null(unit)) {
      converted = switch(unit,
                        "s" = duration,
                        "ms" = duration * 1000,
                        "us" = duration * 1000000,
                        "µs" = duration * 1000000,
                        "ns" = duration * 1000000000,
                        "m" = duration / 60,
                        "h" = duration / 3600,
                        "d" = duration / 86400,
                        stop("Unsupported unit: ", unit))
      return(paste0(format(converted, scientific = FALSE), unit))
    }
    
    units = list(
      list(threshold = 86400, suffix = 'd'),
      list(threshold = 3600, suffix = 'h'),
      list(threshold = 60, suffix = 'm'),
      list(threshold = 1, suffix = 's'),
      list(threshold = 0.001, suffix = 'ms'),
      list(threshold = 0.000001, suffix = 'µs'),
      list(threshold = 0.000000001, suffix = 'ns')
    )
    
    parts = c()
    
    for (unit in units) {
      threshold = unit$threshold
      suffix = unit$suffix
      
      if (duration >= threshold || 
          (length(parts) == 0 && threshold == 0.000000001)) {
        if (threshold >= 1) {
          value = as.integer(duration %/% threshold)
          duration = duration %% threshold
        } else {
          value = as.integer((duration / threshold) %% 1000)
        }
        if (value > 0 || (length(parts) == 0 && threshold == 0.000000001)) {
          parts = c(parts, paste0(value, suffix))
        }
        if (length(parts) == 2) break
      }
    }
    
    if (length(parts) > 0) paste(parts, collapse = ' ') else 'less than 1ns'
  }
  
  print_summary = function(sort = TRUE, unit = NULL) {
    cat('\n--- Timing Summary ---\n')
    
    if (length(env$timings) == 0) {
      cat('no timings recorded.\n')
      return(invisible(NULL))
    }
    
    if (sort) {
      durations = sapply(env$timings, function(x) x$duration)
      items = names(env$timings)[order(durations, decreasing = TRUE)]
    } else {
      items = names(env$timings)
    }
    
    total_time = sum(sapply(env$timings, function(x) x$duration))
    
    for (msg in items) {
      info = env$timings[[msg]]
      duration = info$duration
      percentage = if (total_time > 0) (duration / total_time) * 100 else 0
      max_mem = info$max_mem
      mem_percent = info$mem_percent


      status = if (info$aborted) 'aborted after' else 'took'
      time_str = format_time(duration, unit)
      
      cat(sprintf('%s %s %s (%.1f%%) using %.2f GiB (%.1f%%)\n', msg, status, time_str, percentage,max_mem,mem_percent))
    }
    
    cat(sprintf('\nTotal time: %s\n', format_time(total_time, unit)))
  }
  
  to_dataframe = function(sort = TRUE, unit = NULL) {
    if (length(env$timings) == 0) {
      return(data.frame(
        operation = character(0),
        duration = numeric(0),
        duration_unit = character(0),
        aborted = logical(0),
        percentage = numeric(0),
        memory = numeric(0),
        memory_unit = character(0),
        percent_mem = numeric(0)
      ))
    }
    
    if (sort) {
      durations = sapply(env$timings, function(x) x$duration)
      items = names(env$timings)[order(durations, decreasing = TRUE)]
    } else {
      items = names(env$timings)
    }
    
    total = sum(sapply(env$timings, function(x) x$duration))
    
    ops = items
    durs = sapply(items, function(msg) env$timings[[msg]]$duration)

    
    if (!is.null(unit)) {
      durs = switch(unit,
                   "s" = durs,
                   "ms" = durs * 1000,
                   "us" = durs * 1000000,
                   "µs" = durs * 1000000,
                   "ns" = durs * 1000000000,
                   "m" = durs / 60,
                   "h" = durs / 3600,
                   "d" = durs / 86400,
                   stop("Unsupported unit: ", unit))
      duration_unit = unit
    } else {
      duration_unit = "s" 
    }
    
    aborts = sapply(items, function(msg) env$timings[[msg]]$aborted)
    pcts = sapply(items, function(msg) {
      if (total > 0) (env$timings[[msg]]$duration / total) * 100 else 0
    })
    memory = sapply(items, function(msg) env$timings[[msg]]$max_mem)
    memory_unit = sapply(items, function(msg) {"GiB"})
    percent_mem = sapply(items, function(msg) env$timings[[msg]]$mem_percent)
    data.frame(
      operation = ops,
      duration = durs,
      duration_unit = duration_unit,
      aborted = aborts,
      percentage = pcts,
      memory = memory,
      memory_unit = memory_unit,
      percent_mem = percent_mem
    )
  }
  structure(list(
    with_timer = with_timer,
    print_summary = print_summary,
    to_dataframe = to_dataframe
  ), class = "TimerCollection")
}

system_info = function() {
  cp = parallel::detectCores(logical = FALSE)
  cl = parallel::detectCores(logical = TRUE)
  
  mem_info = system('free -g', intern = TRUE)
  mem_line = strsplit(mem_info[2], "\\s+")[[1]]
  tm = as.numeric(mem_line[2])
  am = as.numeric(mem_line[4])
  
  cat('\n--- System Information ---\n')
  cat(sprintf('Node: %s\n', Sys.info()['nodename']))
  cat(sprintf('CPU: %s physical cores, %s logical cores\n', cp, cl))
  cat(sprintf('Memory: %.1f GB available / %.1f GB total\n', am, tm))
}


read_h5ad_obs <- function(file_path) {
  h5_file <- H5File$new(file_path, mode = "r")
  obs_group <- h5_file[["obs"]]
  
  obs_list <- list()
  for (col_name in list.datasets(obs_group)) {
    if (col_name != "_index") {
      obs_list[[col_name]] <- obs_group[[col_name]][]
    }
  }
  
  if ("__categories" %in% list.groups(obs_group)) {
    cat_group <- obs_group[["__categories"]]
    for (cat_col in list.datasets(cat_group)) {
      codes <- obs_list[[cat_col]]
      levels <- cat_group[[cat_col]][]
      obs_list[[cat_col]] <- factor(levels[codes + 1])
    }
  }
  
  obs_df <- as.data.frame(obs_list)
  rownames(obs_df) <- obs_group[["_index"]][]
  
  h5_file$close()
  return(obs_df)
}
