library('futile.logger')

initialize_swarm <- function(par, fn, s) {
  n <- length(par)
  swarm <- list()
  swarm$particles <- replicate(n = s, simplify = FALSE, expr = {
    position <- best_position <- sample(x = c(FALSE, TRUE), size = n, replace = TRUE)
    velocity <- velocity_0 <- velocity_1 <- numeric(length = n)
    value <- best_value <- -Inf
    particle <- list('position' = position,
                     'velocity_0' = velocity_0,
                     'velocity_1' = velocity_1,
                     'velocity' = velocity,
                     'value' = value,
                     'best_position' = best_position,
                     'best_value' = best_value)
    
    particle$value <- particle$best_value <- fn(particle$position)
    particle
  })
  
  swarm$best_value <- -Inf
  swarm$best_position <- logical(n)
  swarm
}

update_swarm_ibpso <- function(fn, swarm,  w, c_p, c_g,ibpso =NA) {
  
  for (i in seq_along(swarm$particles)) {
    particle <- swarm$particles[[i]]
    particle$value <- fn(particle$position)
    n <- length(particle$position)
    
    if (is.na(particle$value))
      flog.error('The function value is NA!', name = 'console')
    
    ## Compare current performance to its best performance
    if (particle$value > particle$best_value) {
      particle$best_value <- particle$value
      particle$best_position <- particle$position
    }
    
    ## Compare current performance to globally best performance
    if (particle$value > swarm$best_value) {
      swarm$best_value <- particle$value
      swarm$best_position <- particle$position
    }
    
    ## Change velocities of the particle
    d11 <- d01 <- d12 <- d02 <- logical(length = n)
    c1r1 <- runif(n = n, min = 0, max = c_p)
    c2r2 <- runif(n = n, min = 0, max = c_g)
    pbp <- particle$best_position
    
    if (any(is.na(ibpso))){
      gbp <- swarm$best_position
    } else {
      gbp <- ibpso
    }
    
    d11[pbp] <- c1r1[pbp]
    d01[pbp] <- -c1r1[pbp]
    d01[!pbp] <- c1r1[!pbp]
    d11[!pbp] <- -c1r1[!pbp]
    d12[gbp] <- c2r2[gbp]
    d02[gbp] <- -c2r2[gbp]
    d02[!gbp] <- c2r2[!gbp]
    d12[!gbp] <- -c2r2[!gbp]
    
    particle$velocity_1 <- w * particle$velocity_1 + d11 + d12
    particle$velocity_0 <- w * particle$velocity_0 + d01 + d02
    
    ## Change velocities of bit-change
    pp <- particle$position
    particle$velocity[!pp] <- particle$velocity_1[!pp]
    particle$velocity[pp] <- particle$velocity_0[pp]
    
    ## Update position and value
    s <- 1.0 / (1.0 + exp(-particle$velocity))
    logsub <- runif(n = n) < s
    particle$position[logsub] <- !particle$position[logsub]
    ## leave unchanged otherwise
    swarm$particles[[i]] <- particle
  }
  swarm
}


ibpso <- function (par, fn, ..., control = list(), debug = FALSE) 
{
  flog.logger("console", INFO)
  logfile <- "bpso.log"
  if (file.exists(logfile)) 
    file.remove(logfile)
  flog.logger("logfile", DEBUG, appender = appender.file(logfile))
  if (!debug) {
    flog.threshold(0L, name = "debug")
  }
  
  n <- length(par)
  con <- list(trace = TRUE, REPORT = 10L, fnscale = 1, maxit = 1000L, 
              s = floor(10 + 2 * sqrt(n)), w = 1/(2 * log(2)), c.p = 1 + 
                log(2), c.g = 1 + log(2), v.max = NA_real_, maxit.stagnate = Inf, 
              trace.stats = FALSE)
  names_con <- names(con)
  con[(names_control <- names(control))] <- control
  if (length(no_names <- names_control[!(names_control %in% 
                                         names_con)])) 
    warning("unknown names in control: ", paste(no_names, 
                                                collapse = ", "))
  fn1 <- function(par) fn(par)/con$fnscale
  trace <- con$trace > 0L
  REPORT <- con$REPORT
  trace_stats <- con$trace.stats
  if (!is.na(con$v.max)) {
    .NotYetUsed("v.max", error = FALSE)
  }
  swarm <- initialize_swarm(par = par, fn = fn1, s = con$s)
  if (trace_stats) {
    trace_stats_ret <- list()
    tsc <- 1L
  }
  it <- 1L
  it_stagnate <- 0L
  
  while (TRUE) {
    best_value <- swarm$best_value
    best_position <- swarm$best_position
    
    swarm <- update_swarm(fn = fn1, swarm = swarm, w = con$w, 
                          c_p = con$c.p, c_g = con$c.g)
    
    #ibpso 
    cond_ibpso <- isTRUE(all.equal(target = best_position, current = swarm$best_position))
    if (cond_ibpso){
      g_best <- rowSums(sapply(swarm$particles,function(x) x$best_position))/con$s > 0.5
      swarm <- update_swarm_ibpso(fn = fn1, swarm = swarm, w = con$w, 
                            c_p = con$c.p, c_g = con$c.g,ibpso = g_best)
    } else {
      swarm <- update_swarm_ibpso(fn = fn1, swarm = swarm, w = con$w, 
                                  c_p = con$c.p, c_g = con$c.g)
    }
    
    
    cond <- isTRUE(all.equal(target = best_value, current = swarm$best_value))
    if (cond) {
      it_stagnate <- it_stagnate + 1L
    }
    else {
      it_stagnat <- 0L
    }
    
    if (trace) {
      if (it%%REPORT == 0L) {
        tmp <- lapply(X = swarm$particles, FUN = function(x) 1/(1 + 
                                                                  exp(-x$velocity)))
        prob_bit_change <- mean(Reduce(f = `+`, 
                                       x = tmp)/length(tmp))
        flog.info("iteration = %4d | value = %.3f | P(change) = %1.3f", 
                  it, swarm$best_value, prob_bit_change, name = "info")
      }
    }
    if (trace_stats) {
      if (it%%REPORT == 0L) {
        trace_stats_ret[[tsc]] <- list(iteration = it, 
                                       par = swarm$best_position, value = swarm$best_value, 
                                       f = vapply(swarm$particles, function(x) x$value, 
                                                  FUN.VALUE = numeric(1L)), x = do.call(rbind, 
                                                                                        lapply(swarm$particles, function(x) x$position)))
        tsc <- tsc + 1L
      }
    }
    if (it == con$maxit || it_stagnate == con$maxit.stagnate) {
      break
    }
    it <- it + 1L
  }
  
  list(par = swarm$best_position, value = swarm$best_value, 
       stats = if (trace_stats) trace_stats_ret else NULL)
}
