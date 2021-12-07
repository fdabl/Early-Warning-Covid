source('Code/helpers.R')
cols <- brewer.pal(5, 'Set1')

######################################
#### Figure 4 & 5 (Illustration) #####
######################################

N <- 1e6
I0 <- N * 0.01
S0 <- N - I0
gamma <- 1 / 10
sigma <- 1 / 5.2
eta <- 2e-5
eta_start <- 50


#' Creates Figure 4 and 5 in the paper
create_illustration <- function(
  nr_constant, nr_increase, nr_addition = 150, nsims = 50, lowest_point = 0.50,
  ews1_seq = seq(0, 400, 100), ews2_seq = seq(-1, 1, 0.50), seq_step = 50, ews_seq_step = 50, seed = 2,
  bw = 50, ws = 50, ews = c('variance', 'autocorrelation'), xlab = 'Days', transition = TRUE,
  title = 'Second outbreak', ylab = 'Reported cases', ylab2 = expression(R[t]), ylab_ews = TRUE,
  colgray = 'gray76', logarithm = TRUE, total_ews = FALSE, ews_logarithm = FALSE,
  backward_only = TRUE, cut_first_window = FALSE
  ) {
  
  set.seed(seed)
  
  nr_decrease <- 25
  times <- seq(nr_decrease + nr_constant + nr_increase + nr_addition)
  
  Rt <- get_forcing(
    3, lowest_point, nr_decrease, nr_constant,
    ifelse(transition, nr_increase, 0), length(times)
  )
  beta_par_t <- Rt * gamma
  eta_par_t <- c(rep(0, eta_start), rep(eta, length(times) - eta_start))
  
  sim <- create_SEIR(
    times, beta_par_t = beta_par_t, eta_par_t = eta_par_t,
    params = c(sigma = sigma, gamma = gamma, eta = 0, rho = 1),
    init = c(S0 = S0, E0 = 0, I0 = I0, R0 = 0, N = N)
  )
  
  simres <- simulate(sim, nsim = nsims, format = 'data.frame')
  nnames <- colnames(simres)
  nnames[2] <- 'iteration'
  colnames(simres) <- nnames
  simres$iteration <- as.factor(simres$iteration)
  
  end_first <- nr_decrease + nr_constant
  start_second <- end_first + nr_increase
      
  simres$Rt <- rep(Rt, each = nsims)
  simres$Rt_low <- ifelse(total_ews, 1, end_first)
  simres$Rt_critical <- ifelse(total_ews, max(times) - 1, start_second)
  sel <- seq(end_first, start_second)
  
  # For plotting
  simres[simres$time == 1, 'reports'] <- simres[simres$time == 2, 'reports']
  
  simres_Rt <- simres[simres$iteration == 1, ]
  n_simres <- nrow(simres_Rt) + 1
  
  simres_mean <- simres %>% 
    group_by(time) %>% 
    summarize(mean_reports = mean(reports))
  
  if (bw == 'simulation') {
    df_ews <- simres %>% 
      group_by(iteration) %>%
      mutate(reports = reports - simres_mean$mean_reports)
  } else {
    df_ews <- simres %>% group_by(iteration)
  }
  
  
  shorten <- function(x) {
    to_remove <- ws - 1
    len <- length(x)
    
    # If there are enough data points left after
    # removing the first rolling window, do it
    # otherwise return a vector of NAs
    if ((len - to_remove) >= 2) {
      x[seq(to_remove)] <- NA
      return(x)
      
    } else {
      return(rep(NA, len))
    }
  }
  
  get_ews_short <- function(y, bw, ws, type, detrending = TRUE, backward_only = TRUE) {
    x <- get_ews(y, bw, ws, type, detrending, backward_only)
    
    if (cut_first_window) {
      return(lapply(x, shorten))
    }
    x
  }
  
  df_ews <- df_ews %>%
    summarize(
      index = seq(Rt_low[1], Rt_critical[1]),
      variance = get_ews_short(
        reports[index], bw = bw, ws = ws,
        type = 'uniform', backward_only = backward_only
      )[[ews[1]]],
      autocorrelation = get_ews_short(
        reports[index], bw = bw, ws = ws,
        type = 'uniform', backward_only = backward_only
      )[[ews[2]]]
    ) %>% 
    melt(c('index', 'iteration')) %>% 
    mutate(
      variable = factor(variable, labels = ews),
      Rt_low = end_first,
      Rt_critical = start_second
    )
  
  df_ews_mean <- df_ews %>% 
    group_by(variable, index) %>% 
    summarize(mean_ews = mean(value))
  
  simres$Rt_low <- end_first
  simres$Rt_critical <- start_second
  
  simres_Rt$Rt_low <- end_first
  simres_Rt$Rt_critical <- start_second
  
  p1 <- ggplot(simres, aes(x = time, y = reports)) +
    geom_line(color = colgray, aes(group = iteration)) +
    geom_line(data = simres_mean, aes(x = time, y = mean_reports), color = 'black', size = 1) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 3000)) +
    scale_x_continuous(expand = c(0, 0), limits = c(1, n_simres), breaks = c(NA, seq(25, n_simres, seq_step))) +
    geom_vline(aes(xintercept = Rt_low), col = cols[2], size = 1.25) + 
    geom_vline(aes(xintercept = Rt_critical), col = cols[2], size = 1.25) + 
    ggtitle(title) +
    xlab('') +
    ylab(ylab) +
    ptheme + #margin +
    theme(legend.position = 'none')
  
  if (logarithm) {
    p1 <- p1 + 
      scale_y_log10(expand = c(0, 0), limits = c(1, 5000), breaks = c(10, 100, 1000, 5000))
  }
  
  p2 <- ggplot(simres_Rt, aes(x = time, y = Rt)) +
    geom_hline(aes(yintercept = 1), linetype = 'dashed') +
    geom_line(color = 'gray66', size = 1.5) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 3)) +
    scale_x_continuous(expand = c(0, 0), limits = c(1, n_simres), breaks = c(NA, seq(25, n_simres, seq_step))) +
    geom_vline(aes(xintercept = Rt_low), col = cols[2], size = 1.25) + 
    geom_vline(aes(xintercept = Rt_critical), col = cols[2], size = 1.25) + 
    xlab('') +
    ylab(ylab2) +
    ptheme
  
  
  plot_ews <- function(df, ews, ylimits, ybreaks, ylab, ews_logarithm, xlab = 'Days', legend = TRUE) {
    df <- df[df$variable == ews, ]
    df_mean <- df_ews_mean[df_ews_mean$variable == ews, ]
    
    p <- ggplot(df, aes(x = index, y = value)) +
      geom_line(color = colgray, aes(group = iteration)) +
      geom_line(data = df_mean, aes(x = index, y = mean_ews), color = 'black', size = 1) +
      scale_y_continuous(
        expand = c(0, 0), limits = ylimits,
        breaks = ybreaks
      ) +
      scale_x_continuous(
        expand = c(0, 0), limits = c(min(df$index), max(df$index)),
        breaks = c(1, seq(25, max(df$index), ews_seq_step))
      ) +
      xlab(xlab) +
      ylab(ylab) +
      ptheme + 
      theme(
        legend.position = c(0.1, 0.8),
        legend.title = element_text(size = 10), 
        legend.key.size = unit(.80, 'lines'),
        legend.text = element_text(size = 8)
      ) +
      guides(color = guide_legend(override.aes = list(size = 0.50)))
    
    
    if (!legend) {
      p <- p + theme(legend.position = 'none')
    }
    
    if (ews_logarithm) {
      p <- p +
        scale_y_log10(
          expand = c(0, 0), limits = c(1, 1000000),
          breaks = c(10, 100, 1000, 10000, 100000, 1000000)
        ) +
        geom_vline(aes(xintercept = Rt_low), col = cols[2], size = 1.25) + 
        geom_vline(aes(xintercept = Rt_critical), col = cols[2], size = 1.25)
    }
    
    p
  }
    
  p_var <- plot_ews(
    df_ews, ews[1], c(min(ews1_seq), max(ews1_seq)), ews_logarithm = ews_logarithm,
    ews1_seq, ylab = ifelse(ylab_ews, str_to_title(ews[1]), ''), legend = FALSE, xlab = ''
  )
  p_aut <- plot_ews(
    df_ews, ews[2], c(min(ews2_seq), max(ews2_seq)), ews2_seq, ews_logarithm = FALSE,
    ylab = ifelse(ylab_ews, str_to_title(ews[2]), ''), legend = FALSE, xlab = xlab
  )
  
  list('A' = p1, 'B' = p2, 'C' = p_var, 'D' = p_aut)
}


# Combines figures to be like Figure 4 & 5 in the paper
combine_figures <- function(
  nr_constant, nr_increase, bw = 50, ws = 50, lowest_point = 0.50,
  ews1_seq = seq(0, 400, 100), ews2_seq = seq(-1, 1, 0.50), ews_seq_step = 50,
  ews = c('variance', 'autocorrelation'), seed = 1, colgray = 'gray86',
  nr_addition = 150, nsims = 50, seq_step = 50, total_ews = FALSE,
  ews_logarithm = FALSE, backward_only = TRUE, cut_first_window = FALSE
  ) {
  
  margin <- theme(plot.margin = unit(c(0.25, 0.5, 0.2, 0.2), 'cm'))
  margin_bottom <- theme(plot.margin = unit(c(0, 0.5, 0.2, 0.2), 'cm'))
  
  margin2 <- theme(plot.margin = unit(c(0.25, 1, 0.2, -0.5), 'cm'))
  margin2_bottom <- theme(plot.margin = unit(c(0, 1, 0.2, -0.5), 'cm'))
  
  p <- create_illustration(
    nr_constant = nr_constant, nr_increase = nr_increase, bw = bw, ws = ws,
    seed = seed, ews1_seq = ews1_seq, ews2_seq = ews2_seq, xlab = 'Days', ews = ews,
    colgray = colgray, nr_addition = nr_addition, nsims = nsims, seq_step = seq_step,
    lowest_point = lowest_point, total_ews = total_ews, ews_seq_step = ews_seq_step,
    ews_logarithm = ews_logarithm, backward_only = backward_only,
    cut_first_window = cut_first_window
  )
  
  pn <- create_illustration(
    nr_constant = nr_constant, nr_increase = nr_increase, bw = bw, ws = ws, ews = ews,
    seed = seed, ews1_seq = ews1_seq, ews2_seq = ews2_seq, transition = FALSE,
    title = 'No second outbreak', ylab = '', ylab2 = '', xlab = 'Days', ylab_ews  = FALSE,
    colgray = colgray,  nr_addition = nr_addition, nsims = nsims, seq_step = seq_step,
    lowest_point = lowest_point, total_ews = total_ews, ews_seq_step = ews_seq_step,
    ews_logarithm = ews_logarithm, backward_only = backward_only,
    cut_first_window = cut_first_window
  )
  
  gA <- ggplotGrob(p$A + margin)
  gB <- ggplotGrob(p$B + margin_bottom)
  gC <- ggplotGrob(p$C + margin_bottom)
  gD <- ggplotGrob(p$D + margin_bottom)
  
  gnA <- ggplotGrob(pn$A + margin2)
  gnB <- ggplotGrob(pn$B + margin2_bottom)
  gnC <- ggplotGrob(pn$C + margin2_bottom)
  gnD <- ggplotGrob(pn$D + margin2_bottom)
  
  ggarrange(
    rbind(gA, gB, gC, gD),
    rbind(gnA, gnB, gnC, gnD)
  )
}


nsims <- 50
width <- 8 - 1/4
height <- 10 - 3/4

########################
#### 200 Days constant
########################
pdf('Figures/Figure-4a-raw.pdf', width = width, height = height)
combine_figures(
  nsims = nsims, nr_constant = 200, nr_increase = 200, ws = 50, bw = 20,
  seq_step = 100, ews1_seq = seq(0, 200, 50), backward_only = TRUE,
  ews_seq_step = 25, nr_addition = 350 - 150, total_ews = FALSE, ews_logarithm = FALSE,
  cut_first_window = TRUE
)
dev.off()

########################
#### 50 Days constant
########################
pdf('Figures/Figure-4b-raw.pdf', width = width, height = height)
combine_figures(
  nsims = nsims, nr_constant = 50, nr_increase = 200, bw = 5,
  ws = 50, ews_seq_step = 25, ews1_seq = seq(0, 750, 250),
  nr_addition = 150, total_ews = FALSE, ews_logarithm = FALSE,
  cut_first_window = TRUE
)
dev.off()
