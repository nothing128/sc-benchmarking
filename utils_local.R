suppressPackageStartupMessages({
  library(parallel)
  library(benchmarkme)
  library(pryr)
})

TimerCollection = function(silent = TRUE) {
  env = environment()
  env$timings = list()
  
  with_timer = function(message, expr) {
    start = Sys.time()
    
    if (!silent) {
      cat(paste0(message, '...\n'))
    }
    
    result = NULL
    aborted = FALSE
    
    tryCatch({
      result = eval(substitute(expr), parent.frame())
    }, error = function(e) {
      aborted = TRUE
      stop(e)
    }, finally = {
      duration = as.numeric(difftime(Sys.time(), start, units = 'secs'))
      env$timings[[message]] = list(
        duration = duration,
        aborted = aborted
      )
    })
    
    return(result)
  }
  
  format_time = function(duration) {
    duration = as.numeric(duration)
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
  
  print_summary = function(sort = TRUE) {
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
      status = if (info$aborted) 'aborted after' else 'took'
      time_str = format_time(duration)
      
      cat(sprintf('%s %s %s (%.1f%%)\n', msg, status, time_str, percentage))
    }
    
    cat(sprintf('\nTotal time: %s\n', format_time(total_time)))
  }
  
  to_dataframe = function(sort = TRUE) {
    if (length(env$timings) == 0) {
      return(data.frame(
        operation = character(0),
        duration_seconds = numeric(0),
        aborted = logical(0),
        percentage = numeric(0)
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
    aborts = sapply(items, function(msg) env$timings[[msg]]$aborted)
    pcts = sapply(items, function(msg) {
      if (total > 0) (env$timings[[msg]]$duration / total) * 100 else 0
    })
    
    data.frame(
      operation = ops,
      duration_seconds = durs,
      aborted = aborts,
      percentage = pcts
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