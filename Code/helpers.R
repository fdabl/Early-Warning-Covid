library('pomp')
library('grid')
library('trend')
library('dplyr')
library('tidyr')
library('spaero')
library('ggpubr')
library('ggrepel')
library('ggplot2')
library('stringr')
library('EpiNow2')
library('COVID19')
library('latex2exp')
library('covid19us')
library('doParallel')
library('countrycode')
library('RColorBrewer')
library('covidregionaldata')


# ggplot themes
ptheme <- theme(
  plot.title = element_text(hjust = .5),
  text = element_text(size = 14),
  legend.key = element_blank(),
  panel.background = element_blank(),
  legend.position = 'top',
  axis.line.x = element_line(),
  axis.line.y = element_line(),
  strip.background = element_blank()
)

ptheme2 <- theme(
  plot.title = element_text(hjust = .5),
  text = element_text(size = 14),
  legend.key = element_blank(),
  legend.title = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  panel.border = element_blank(),
  legend.position = 'top',
  axis.line = element_blank(),
  strip.background = element_blank()
)



#' Function that creates the pomp SEIR object
create_SEIR <- function(
  times,
  beta_par_t,
  eta_par_t,
  t0 = min(times),
  params = c(sigma = 1 / 5.2, gamma = 1 / 10, eta = 0, rho = 1),
  init = c(S0 = 1e5 - 100, E0 = 0, I0 = 100, R0 = 0, N = 1e5)
) {
  
  SEIR_template <- list(
    
    exposeS = list(
      'rate = ((beta_par_t * S0 / S) * I / N + eta_par_t) * S;',
      c(S = -1, E = 1, I = 0, R = 0, cases = 0)
    ),
    
    infectI = list(
      'rate = sigma * E;',
      c(S = 0, E = -1, I = 1, R = 0, cases = 0)
    ),
    
    recoverI = list(
      'rate = gamma * I;',
      c(S = 0, E = 0, I = -1, R = 1, cases = 1)
    )
  )
  
  rinit <- Csnippet(
    '
    S = S0;
    E = E0;
    I = I0;
    R = R0;
    cases = 0;
    '
  )
  
  rmeas <- Csnippet('reports = rbinom(cases, rho);')
  rprocess <- do.call(pomp::gillespie_hl, SEIR_template)
  
  covar <- data.frame(
    'time' = times,
    'beta_par_t' = beta_par_t,
    'eta_par_t' = eta_par_t
  )
  
  covar <- covariate_table(covar, times = 'time')
  
  params <- c(params, init)
  pomp::simulate(
    times = times, t0 = t0, params = params,
    paramnames = names(params), rprocess = rprocess,
    rmeasure = rmeas, rinit = rinit, covar = covar,
    statenames = c('S', 'E', 'I', 'R', 'cases'),
    accumvars = 'cases', obsnames = 'reports'
  )
}


#' Function that specifies the forcing of R(t)
get_forcing <- function(
  R0, lowest_point, nr_decrease, nr_constant, nr_increase, nr_total
) {
  
  end_first <- seq(R0, lowest_point, length.out = nr_decrease)
  quiet_period <- seq(lowest_point, lowest_point, length.out = nr_constant + 1)[-1]
  
  if (nr_increase == 0) {
    
    Rt <- c(end_first, quiet_period)
    Rt <- c(Rt, rep(lowest_point, nr_total - length(Rt)))
    
  } else {
    reaching_criticality <- seq(lowest_point, 1, length.out = nr_increase + 1)[-1]
    
    rest <- nr_total - nr_decrease - nr_constant - nr_increase
    
    slope <- reaching_criticality[2] - reaching_criticality[1]
    moving_further <- 1 + seq(rest) * slope
    
    Rt <- c(end_first, quiet_period, reaching_criticality, moving_further)
  }
  
  Rt[seq(nr_total)]
}


#' Function that computes the early warning indicators
get_ews <- function(
  y, bw = 5, ws = 5, type = 'gaussian',
  detrending = TRUE, backward_only = TRUE
) {
  
  stats <- get_stats(
    y,
    center_trend = ifelse(bw == 'simulation', 'assume_zero', 'local_constant'),
    center_kernel = type,
    center_bandwidth = bw,
    stat_trend = 'local_constant',
    stat_kernel = type,
    stat_bandwidth = ws,
    backward_only = backward_only
  )
  
  if (!detrending) {
    stats_ndt <- get_stats(
      y,
      center_trend = 'assume_zero',
      center_kernel = type,
      center_bandwidth = bw,
      stat_trend = 'local_constant',
      stat_kernel = type,
      stat_bandwidth = ws,
      backward_only = backward_only
    )
    
    # These indicators require an estimate of the mean!
    stats_ndt$stats$mean <- stats$stats$mean
    stats_ndt$stats$index_of_dispersion <- stats$stats$index_of_dispersion
    stats_ndt$stats$coefficient_of_variation <- stats$stats$coefficient_of_variation
    
    return(stats_ndt$stats)
  }
  
  stats$stats
}


#' Simulates from the SEIR model and computes early warning indicators
simulate_ews <- function(
  lowest_point, nr_decrease,
  nr_constant, nr_increase,
  R0 = 3, eta_start = 50, eta = 1e-5, type = 'gaussian',
  nsims = 250, bw = 5, ws = 5, detrending = TRUE, fit_ARMA = FALSE,
  nr_surrogates = 250, cut_first_window = FALSE
) {
  
  N <- 1e6
  I0 <- N * 0.01
  S0 <- N - I0
  
  rho <- 1
  gamma <- 1 / 10
  sigma <- 1 / 5.2
  
  grid <- expand.grid(
    lowest_point = lowest_point,
    nr_decrease = nr_decrease,
    nr_constant = nr_constant,
    nr_increase = nr_increase
  )
  
  create_simres <- function(
    R0, lowest_point, nr_decrease, nr_constant,
    nr_increase, times
  ) {
    
    Rt <- get_forcing(
      R0, lowest_point, nr_decrease, nr_constant,
      nr_increase, length(times)
    )
    
    beta_par_t <- Rt * gamma
    eta_par_t <- c(rep(0, eta_start), rep(eta, length(times) - eta_start))
    
    sim <- create_SEIR(
      times, beta_par_t = beta_par_t, eta_par_t = eta_par_t,
      params = c(sigma = sigma, gamma = gamma, eta = eta, rho = rho),
      init = c(S0 = S0, E0 = 0, I0 = I0, R0 = 0, N = N)
    )
    
    simres <- simulate(sim, nsim = nsims, format = 'data.frame')
    simres
  }
  
  
  get_early_warnings <- function(simres, nr_decrease, nr_constant, nr_increase) {
    
    kendalls <- c()
    ylengths <- c()
    pvals <- c()
    
    # Detrend using the mean of the simulations
    if (bw == 'simulation') {
      
      ymean <- simres %>% 
        group_by(time) %>% 
        summarize(ymean = mean(cases))
      
      ymean <- ymean$ymean
      
    } else {
      ymean <- rep(0, 2000)
    }
    
    for (i in seq(nsims)) {
      sim <- simres[simres$.id == i, ]
      y <- sim$cases
      
      end_first <- nr_decrease + nr_constant
      start_second <- end_first + abs(nr_increase)
      
      ymeansel <- ymean[seq(end_first, start_second)]
      ysel <- y[seq(end_first, start_second)]
      ylengths <- c(ylengths, length(ysel))
      
      ysel <- ysel - ymeansel
      
      if (bw == 'dynamic') {
        bws <- list('25' = 3, '50' = 5, '100' = 10, '200' = 20)
        bw <- bws[[as.character(nr_constant)]]
      }
      
      ews <- get_ews(ysel, bw = bw, ws = ws, detrending = detrending, type = type)
      
      if (bw == 'simulation') {
        ews$mean <- ymeansel
        ews$index_of_dispersion <- ews$variance / ews$mean
        ews$coefficient_of_variation <- sqrt(ews$variance) / ews$mean
      }
      
      if (fit_ARMA) {
        ys <- ysel - ews$mean
        len <- length(ys)
        
        # get best fitting ARMA(p, q) model
        arima_coefs <- select_arima(ys)
        
        surrogates <- tryCatch(
          replicate(nr_surrogates, arima.sim(model = arima_coefs, n = len)),
          error = function(e) {
            NULL
          }
        )
        
        if (!is.null(surrogates)) {
          stats_ews <- make_inference(surrogates, ews, ws, bw, cut_first_window)
        } else {
          stats_ews <- cbind(rep(NA, 10), rep(NA, 10), rep(NA, 10))
        }
        
        sews <- c(stats_ews[, 1], stats_ews[, 2], stats_ews[, 3])
        names(sews) <- c(
          paste0(names(ews), '_tau'),
          paste0(names(ews), '_pval'),
          paste0(names(ews), '_pvalabs')
        )
        
        kendalls <- rbind(kendalls, sews)
        
      } else {
        
        # AUC route
        get_tau <- function(x) {
          x <- na.omit(x) # Remove NAs and NaNs
          x <- x[x != -Inf & x != Inf] # Remove -Inf for decay time
          
          if (length(x) < 3) {
            return(NA)
          } else {
            mk.test(x)$estimates[3]
          }
        }
        
        shorten <- function(x) {
          to_remove <- ws - 1
          len <- length(x)
          
          # If there are enough data points left after
          # removing the first rolling window, do it
          # otherwise return a vector of NAs
          if ((len - to_remove) >= 2) {
            return(x[-seq(to_remove)])
            
          } else {
            return(rep(NA, len))
          }
        }
        
        if (cut_first_window) {
          ews <- lapply(ews, shorten)
        }
      
        ix <- seq(length(ysel))
        
        kendall <- sapply(ews, get_tau)
        kendalls <- rbind(kendalls, kendall)
        
      }
    }
    
    result <- cbind(
      kendalls, times = seq(nsims),
      lowest_point = lowest_point,
      nr_decrease = nr_decrease,
      nr_constant = nr_constant,
      nr_increase = nr_increase,
      ylen_mean = mean(ylengths),
      ylen_sd = sd(ylengths)
    )
    
    result
  }
  
  # Null result is the same for a fixed number of nr_decrease
  stopifnot(length(nr_decrease) == 1 && length(lowest_point) == 1)
  
  times_null <- seq(max(nr_decrease) + max(nr_constant) + max(nr_increase))
  simres_null <- create_simres(R0, lowest_point, nr_decrease, nr_constant[1], 0, times_null)
  
  grid_null <- expand.grid(
    lowest_point = lowest_point,
    nr_constant = nr_constant,
    nr_increase = nr_increase
  )
  
  null_results <- foreach(i = seq(nrow(grid_null)), .combine = 'rbind') %dopar% {
    
    config_null <- grid_null[i, ]
    lo <- config_null$lowest_point
    nr_constant <- config_null$nr_constant
    nr_increase <- config_null$nr_increase
    
    result <- get_early_warnings(simres_null, nr_decrease, nr_constant, nr_increase)
    result <- cbind(result, 'type' = 0)
    
    result
  }
  
  test_results <- foreach (i = seq(nrow(grid)), .combine = 'rbind') %dopar% {
    
    config <- grid[i, ]
    lo <- config$lowest_point
    nr_decrease <- config$nr_decrease
    nr_constant <- config$nr_constant
    nr_increase <- config$nr_increase
    
    times <- nr_decrease + nr_constant + nr_increase
    simres <- create_simres(R0, lo, nr_decrease, nr_constant, nr_increase, seq(times))
    result <- get_early_warnings(simres, nr_decrease, nr_constant, nr_increase)
    result <- cbind(result, 'type' = 1)
    
    result
  }
  
  
  d <- data.frame(cbind(rbind(test_results, null_results), eta = eta, ws = ws))
  d$bw <- bw
  d$smoothing <- type
  rownames(d) <- NULL
  d
}


#' Compute the AUC, using the simulation results from 'simulate_ews'
compute_AUC <- function(d, ews_names, nr_draws = 1e6, AUC = TRUE, alpha = 0.05) {
  
  combinations <- distinct(select(d, starts_with('nr')))
  # test_combinations <- combinations[combinations$nr_increase != 0, ]
  # null_combinations <- combinations[combinations$nr_increase == 0, ]
  
  res <- c()
  
  for (i in seq(nrow(combinations))) {
    
    comb <- as.numeric(combinations[i, ])
    nr_decrease <- comb[1]
    nr_constant <- comb[2]
    nr_increase <- comb[3]
    
    aucs <- c()
    
    for (ews in ews_names) {
      
      if (AUC) {
        
        x0 <- d[
          d$nr_decrease == nr_decrease &
            d$nr_constant == nr_constant &
            d$nr_increase == nr_increase & 
            d$type == 0, ews
          ]
        
        x1 <- d[
          d$nr_decrease == nr_decrease &
            d$nr_constant == nr_constant &
            d$nr_increase == nr_increase &
            d$type == 1, ews
          ]
        
          auc <- mean(
            sample(x1, nr_draws, replace = TRUE) >
            sample(x0, nr_draws, replace = TRUE), na.rm = TRUE
          )
          
        # Get true positive rate
        } else {
          x <- d[
            d$nr_decrease == nr_decrease &
              d$nr_constant == nr_constant &
              d$nr_increase == nr_increase & 
              d$type == 1, ews
          ]
          auc <- mean(x < alpha, na.rm = TRUE)
        }
      
      aucs <- c(aucs, auc)
    }
    
    ylens <- d[
      d$nr_decrease == nr_decrease &
        d$nr_constant == nr_constant &
        d$nr_increase == nr_increase, c('ylen_mean', 'ylen_sd')
      ] %>% distinct %>% as.data.frame
    
    res <- rbind(
      res, c(nr_decrease, nr_constant, nr_increase, aucs, ylens)
    )
  }
  
  res <- data.frame(apply(res, 2, unlist))
  colnames(res) <- c('nr_decrease', 'nr_constant', 'nr_increase', ews_names, 'ylen_mean', 'ylen_sd')
  res
}


#' Estimate R(t) using the R package EpiNow2 (Abbott et al., 2020)
estimate_rt <- function(dat, country, delay = TRUE) {
  
  reported_cases <- dat[dat$country == country, c('date', 'cases_new')]
  
  if (nrow(reported_cases) == 0) {
    stop('Country data frame has zero rows.')
  }
  
  colnames(reported_cases) <- c('date', 'confirm')
  
  # distributions from the literature
  generation_time <- get_generation_time(disease = 'SARS-CoV-2', source = 'ganyani')
  incubation_period <- get_incubation_period(disease = 'SARS-CoV-2', source = 'lauer')
  
  # define reporting delay as lognormal with mean of 4 days and sd of 1 day
  reporting_delay <- list(mean = convert_to_logmean(4, 1),
                          mean_sd = 0.1,
                          sd = convert_to_logsd(4, 1),
                          sd_sd = 0.1,
                          max = 15)
  
  if (delay) {
    delays <- delay_opts(incubation_period, reporting_delay)
  } else {
    delays <- delay_opts(incubation_period)
  }
  
  start <- Sys.time()
  out <- epinow(reported_cases = reported_cases, 
                generation_time = generation_time,
                delays = delays,
                rt = rt_opts(prior = list(mean = 1, sd = 1)),
                # here we define the quality of the gaussian process approximation
                # if the fit to data appears poor try increasing basis_prop and
                # potentially the boundary_scale (see ?gp_opts for details)
                # though this will likely increase runtimes.
                gp = gp_opts(basis_prop = 0.2),
                # in some instances stan chains can get stuck when estimating in 
                # these instances consider using the future fitting mode by passing 
                # `future = TRUE, max_execution_time = 60 * 30` to stan_opts and calling 
                # `future::plan("multiprocess")` prior to running epinow this will time out
                # chains after 30 minutes but still return results from completed chains
                stan = stan_opts(samples = 4000, cores = 4, chains = 4),
                horizon = 0,
                target_folder = "results",
                logs = file.path("logs", Sys.Date()),
                CrIs = c(0.50, 0.60, 0.70, 0.80, 0.90, 0.95),
                return_output = TRUE,
                # output = 'samples',
                verbose = TRUE)
  end <- Sys.time()
  print(end - start)
  out
}


#' Create AUC plot (Figure 6 & 7 in the paper)
create_AUC_plot <- function(
  aucs,
  ylab = 'Early warning indicators',
  xlab = expression('Number of days to reach' ~ R[t] ~ ' = 1'),
  title = 'Area under the Curve',
  subtitle = '', legend_name = 'AUC',
  ytext_size = 12, legend = FALSE, colpallete = rev(brewer.pal(8, 'RdYlBu'))
) {
  
  aucs <- arrange(aucs, nr_constant)
  aucs$nr_constant <- factor(aucs$nr_constant, labels = paste0(unique(aucs$nr_constant), ' days constant'))
  
  p <- ggplot(
    aucs, aes(y = factor(ews), x = factor(nr_increase), fill = value)
  ) +
    geom_tile() +
    scale_fill_gradientn(
      name = legend_name,
      colors = colpallete, limits = c(0, 1)
    ) +
    guides(
      fill = guide_colorbar(
        title.position = 'top', title.hjust = 0.5,
        title.size = 10, barheight = 10, barwidth = 0.75
      )
    ) +
    facet_grid(factor(nr_constant) ~ ., as.table = FALSE) +
    ylab(ylab) +
    xlab(xlab) +
    ggtitle(title) +
    ptheme2 +
    labs(subtitle = subtitle) +
    theme(
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = ytext_size),
      legend.position = ifelse(legend, 'right', c(2, 2)),
      plot.subtitle = element_text(size = 12, hjust = 0.5)
    )
  
  if (ytext_size == 0) {
    p <- p + theme(axis.ticks.y = element_blank())
  }
  
  p
}


normalize <- function(x, ...) (x - min(x, ...)) / (max(x, ...) - min(x, ...))


#' Selects best-fitting ARIMA(p, q) model
#' 
#' @param yprocessed filtered & detrended time-series
#' @returns coefficient of best fitting model
select_arima <- function(yprocessed) {
  models <- list()
  
  w <- 1
  
  for (p in seq(5)) {
    for (q in seq(5)) {
      
      models[[w]] <- tryCatch(
        arima(yprocessed, order = c(p, 0, q)),
        
        # Return NULL when model not stationary
        error = function(e) {
          NULL
        }
      )
      
      w <- w + 1
    }
  }
  
  # Select only those that are stationary
  models_ix <- sapply(models, function(x) !is.null(x))
  models <- models[which(models_ix)]
  
  aics <- sapply(models, AIC)
  m <- models[[which.min(aics)]]
  coefs <- coef(m)
  nnames <- names(coefs)
  
  # Get model coefficients for arima.sim
  arima_coefs <- list(
    intercept = coefs['intercept'],
    ar1 = coefs[grepl('^ar', nnames)],
    ma1 = coefs[grepl('^ma', nnames)]
  )
  
  arima_coefs
}



#' Get boostrapped p-value and observed Kendall's tau
#' 
#' @param surrogates surrogate time-series
#' @param ws rolling window size
#' @param bw bandwidth of centering
#' @param cut_first_window remove estimates for y_i where i < ws
#' 
#' @returns bootstrapped pvalue and observed Kendall's tau
make_inference <- function(surrogates, ews, ws, bw, cut_first_window, ...) {
  
  get_tau <- function(x) {
    x <- na.omit(x) # Remove NAs and NaNs
    x <- x[x != -Inf & x != Inf] # Remove -Inf for decay time
    
    if (length(x) < 3) {
      return(NA)
    } else {
      mk.test(x)$estimates[3]
    }
  }
  
  shorten <- function(x) {
    to_remove <- ws - 1
    len <- length(x)
    
    # If there are enought data points left after
    # removing the first rolling window, do it
    # otherwise return a vector of NAs
    if ((len - to_remove) >= 2) {
      return(x[-seq(to_remove)])
      
    } else {
      return(rep(NA, len))
    }
  }
  
  if (cut_first_window) {
    ews <- lapply(ews, shorten)
  }
  
  taus <- sapply(ews, get_tau)
  
  kendall_null <- sapply(seq(ncol(surrogates)), function(i) {
    ts <- surrogates[, i]
  
    surr_stats <- get_stats(
      ts,
      center_trend = 'local_constant',
      center_kernel = 'uniform',
      center_bandwidth = bw,
      stat_trend = 'local_constant',
      stat_kernel = 'uniform',
      stat_bandwidth = ws,
      backward_only = TRUE
    )
    
    surr_ews <- surr_stats$stats
    
    if (cut_first_window) {
      surr_ews <- lapply(surr_ews, shorten)
    }
    
    taus <- sapply(surr_ews, get_tau)
    taus
  })
  
  ## Observe a tau value at least as large as the null taus
  pvals <- rowMeans(kendall_null >= taus, na.rm = TRUE)
  
  ## Observe a tau value at least as large as the null taus in absolute value
  pvals_abs <- rowMeans(abs(kendall_null) >= abs(taus), na.rm = TRUE)
  
  round(cbind(taus, pvals, pvals_abs), 3)
}



#' Creates data frame from episamples (EpiNow2 object)
create_dat <- function(res, country, cutoff = 0.99, start_date = NULL, end_date = NULL) {
  episamples <- res$samples
  observations <- res$observations
  
  Rdat <- episamples %>% 
    filter(parameter == 'R') %>% 
    group_by(date) %>% 
    summarize(
      R_mean = mean(value),
      R_sd = sd(value),
      R_lo01 = quantile(value, 0.01),
      R_hi99 = quantile(value, 0.99),
      R_lo05 = quantile(value, 0.05),
      R_hi95 = quantile(value, 0.95),
      R_lo10 = quantile(value, 0.10),
      R_hi90 = quantile(value, 0.90),
      crossing_prob = mean(value >= 1)
    )
  
  infections <- episamples %>% 
    filter(parameter == 'infections') %>% 
    group_by(date) %>% 
    summarize(
      infections_mean = mean(value),
      infections_sd = sd(value),
      inf_lo01 = quantile(value, 0.01),
      inf_hi99 = quantile(value, 0.99),
      inf_lo05 = quantile(value, 0.05),
      inf_hi95 = quantile(value, 0.95),
      inf_lo10 = quantile(value, 0.10),
      inf_hi90 = quantile(value, 0.90)
    )
  
  # Compute end of first and start of second wave based on a cutoff
  n <- nrow(Rdat)
  start <- ifelse(is.null(start_date), 1, which(Rdat$date == start_date))
  end <- ifelse(is.null(end_date), n, which(Rdat$date == end_date))
  sel <- seq(start, end) 
  
  # Find those that are below crossing probability
  below_1_ix <- start - 1 + which((1 - Rdat$crossing_prob[sel]) >= cutoff)
  R_mean_min_ix <- which(min(Rdat$R_mean[below_1_ix]) == Rdat$R_mean)
  
  end_first_ix <- R_mean_min_ix
  start_second_ix <- head(which(Rdat$crossing_prob[seq(end_first_ix, end)] >= cutoff), 1)
  start_second_ix <- end_first_ix + start_second_ix
  
  Rdat <- Rdat %>% 
    mutate(
      end_first_ix = end_first_ix,
      end_first_date = Rdat$date[end_first_ix],
      start_second_ix = start_second_ix,
      start_second_date = Rdat$date[start_second_ix],
      country = country
    )
  
  d <- left_join(infections, Rdat, by = 'date')
  d <- left_join(d, observations, by = 'date')
  
  # Find the lowest point after the crossing
  Rbetween1 <- d[seq(end_first_ix, start_second_ix), 'R_mean', drop = TRUE]
  Rmin <- end_first_ix + which(diff(Rbetween1) > 0)[1] - 1
  
  Rdiffs <- diff(filter(d, date >= start_second_date)$R_mean)
  Rmax <- start_second_ix + which(Rdiffs < 0)[1] - 1 # going down again
  
  d %>% 
    mutate(
      first_low_ix = Rmin,
      second_max_ix = Rmax,
      first_low_date = d$date[Rmin],
      second_max_date = d$date[Rmax]
    )
}


#' Plot the R(t) crossing points
plot_crossings <- function(df, title = NULL) {
  
  them <- theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(hjust = 0.50)
  )
  
  p1 <- ggplot(df, aes(x = date, y = confirm)) +
    geom_bar(stat = 'identity') +
    ylab('') +
    ggtitle(title) +
    scale_x_date(date_breaks = '1 month') +
    geom_vline(aes(xintercept = first_low_date), col = 'skyblue') + 
    geom_vline(aes(xintercept = second_max_date), col = 'skyblue') + 
    theme_bw() + them
  
  p2 <- ggplot(df, aes(x = date, y = R_mean)) +
    geom_hline(yintercept = 1, linetype = 'dashed') +
    geom_line() +
    geom_line(aes(y = R_hi99), color = 'gray76') +
    geom_line(aes(y = R_lo01), color = 'gray76') +
    ylab('') +
    geom_vline(aes(xintercept = first_low_date), col = 'skyblue') + 
    geom_vline(aes(xintercept = second_max_date), col = 'skyblue') + 
    scale_x_date(date_breaks = '1 month') +
    theme_bw() + them
  
  gridExtra::grid.arrange(p1, p2, nrow = 2)
}


#' Plot cases with R(t) below it
plot_dat <- function(
  d, ylims, Rlimits, ylab1 = 'Confirmed cases',
  date_breaks = '1 months', ylab2 = expression(R[t]), xlab = 'Date', title = '') {
  
  colv <- cols[2]
  
  p1 <- ggplot(d, aes(x = date, y = infections_mean)) +
    geom_bar(
      aes(x = date, y = confirm), stat = 'identity', fill = 'gray76'
    ) +
    geom_line(size = 1) +
    geom_vline(aes(xintercept = first_low_date), col = colv) + 
    geom_vline(aes(xintercept = second_max_date), col = colv) + 
    xlab('') +
    ylab(ylab1) +
    scale_y_continuous(
      expand = c(0, 0)
    ) +
    scale_x_datetime(
      limits = c(as.POSIXct('2020-03-01'), as.POSIXct('2020-10-01')),
      date_breaks = date_breaks, date_labels = "%b", expand = c(0, 0)
    ) +
    ggtitle(title) +
    ptheme
  
  p2 <- ggplot(d, aes(x = date, y = R_mean)) +
    geom_hline(yintercept = 1, linetype = 'dotted', col = 'gray76') +
    geom_line(size = 1) +
    geom_line(aes(y = R_lo05), col = 'gray76', linetype = 'dashed') +
    geom_line(aes(y = R_hi95), col = 'gray76', linetype = 'dashed') +
    geom_vline(aes(xintercept = first_low_date), col = colv) + 
    geom_vline(aes(xintercept = second_max_date), col = colv) + 
    xlab(xlab) +
    ylab(ylab2) +
    scale_y_continuous(
      limits = Rlimits, expand = c(0, 0)
    ) +
    scale_x_datetime(
      limits = c(as.POSIXct('2020-03-01'), as.POSIXct('2020-10-01')),
      date_breaks = date_breaks, date_labels = "%b", expand = c(0, 0)
    ) +
    ptheme
  
  list('p1' = p1, 'p2' = p2)
}


#' Calculate the early warning indicators from data
#' cut_first_window = TRUE uses the Dakos et al. (2012) approach,
#' that is, does not use the first few data points smaller than the
#' rolling window to estimate the indicator (instead returns NA)
#' cut_first_window = FALSE uses the Drake et al. approach which
#' estimates the indicators also within the first rolling window
get_early_warning <- function(
  dat, windowsizes, bandwidths,
  nr_surrogates = 250, backward_only = TRUE,
  cut_first_window = FALSE
) {
  countries <- unique(dat$country)
  len <- length(countries)
  
  grid <- expand.grid(
    country = countries, windowsize = windowsizes, bandwidth = bandwidths
  )
  
  ews_names <- c(
    'variance', 'variance_first_diff', 'autocovariance', 'autocorrelation',
    'decay_time', 'mean', 'index_of_dispersion', 'coefficient_of_variation',
    'skewness', 'kurtosis'
  )
  
  # From: https://stackoverflow.com/questions/3903157/how-can-i-check-whether-a-function-call-results-in-a-warning
  # To collect convergence problems in ARIMA
  withWarnings <- function(expr) {
    myWarnings <- NULL
    wHandler <- function(w) {
      myWarnings <<- c(myWarnings, list(w))
      invokeRestart("muffleWarning")
    }
    val <- withCallingHandlers(expr, warning = wHandler)
    list(value = val, warnings = myWarnings)
  } 
  
  
  res <- foreach (i = seq(nrow(grid)), .combine = 'rbind') %dopar% {
  # for (i in seq(nrow(grid))) {
    
    config <- grid[i, ]
    country <- as.character(config$country)
    ws <- config$windowsize
    bw <- config$bandwidth
    
    d <- dat[dat$country == country, ]
    sel <- seq(d$first_low_ix[1], d$second_max_ix[1])
    
    y <- d$confirm[sel]
    len <- length(y)
    tindex <- seq(len)
    
    stats <- get_stats(
      y,
      center_trend = 'local_constant',
      center_kernel = 'uniform',
      # center_kernel = 'gaussian',
      center_bandwidth = bw,
      stat_trend = 'local_constant',
      stat_kernel = 'uniform',
      stat_bandwidth = ws,
      backward_only = backward_only
    )
    
    ews <- stats$stats
    ys <- y - stats$centered$center
    
    # get best fitting ARMA(p, q) model
    arima_res <- withWarnings(select_arima(ys))
    arima_coefs <- arima_res$value
    arima_warning <- as.numeric(length(arima_res$warnings) > 0)
    
    nr_arcoefs <- length(arima_coefs$ar1)
    nr_macoefs <- length(arima_coefs$ma1)
    
    surrogates <- tryCatch(
      replicate(nr_surrogates, arima.sim(model = arima_coefs, n = len)),
      error = function(e) {
        NULL
      }
    )
    
    if (!is.null(surrogates)) {
      stats_ews <- make_inference(surrogates, ews, ws, bw, cut_first_window)
    } else {
      stats_ews <- cbind(rep(NA, 10), rep(NA, 10), rep(NA, 10))
    }
    
    c(
      country, ws, bw, nr_arcoefs, nr_macoefs, arima_warning,
      stats_ews[, 1], stats_ews[, 2], stats_ews[, 3]
    )
  }
  
  # Only a single row
  if (is.null(dim(res))) {
    res <- data.frame(t(res))
  }
  
  colnames(res) <- c(
    'country_name',
    'windowsize',
    'bandwidth',
    'nr_arcoefs',
    'nr_macoefs',
    'arima_convergence_warning',
    paste0(ews_names, '_tau'),
    paste0(ews_names, '_pvalue'),
    paste0(ews_names, '_pvalue_abs')
  )
  
  rownames(res) <- NULL
  df <- data.frame(res)
  df[, -1] <- apply(df[, -1], 2, as.numeric)
  df$backward_only <- backward_only
  df$cut_first_window <- cut_first_window
  df
}
