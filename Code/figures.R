source('Code/helpers.R')
cols <- brewer.pal(3, 'Set1')


# Load data with start and end dates of time-series between
# the first and second wave
df <- read.csv('Data/cases-and-rt-data-criterion.csv') %>% 
  mutate(
    date = as.POSIXct(date),
    first_low_date = as.POSIXct(first_low_date),
    second_max_date = as.POSIXct(second_max_date),
  ) %>% 
  filter(
    date >= as.Date('2020-03-01') & date <= as.Date('2020-10-01')
  )


#################################################
#### Figure 2 (Cases in European countries) #####
#################################################
dfa <- df %>% 
  group_by(country) %>% 
  mutate(
    normalized_infections_mean = infections_mean / max(infections_mean),
    normalized_confirm = confirm / max(confirm)
  )

break_fct <-  function(x){
  x <- ceiling(x / 1000) * 1000
  round(seq(0, max(x), length = 4))
}

p <- ggplot(dfa, aes(x = date, y = infections_mean)) +
  geom_bar(
    aes(x = date, y = confirm, fill = 'Reported cases'), stat = 'identity'
  ) +
  geom_line(aes(x = date, y = infections_mean, color = 'Estimated cases')) +
  geom_vline(aes(xintercept = first_low_date), col = cols[2]) +
  geom_vline(aes(xintercept = second_max_date), col = cols[2]) +
  facet_wrap(~ country, scales = 'free', ncol = 4) +
  scale_x_datetime(
    breaks = c(as.POSIXct('2020-03-01'), as.POSIXct('2020-06-01'), as.POSIXct('2020-10-01')),
    limits = c(as.POSIXct('2020-03-01'), as.POSIXct('2020-10-01')), #date_breaks = '2 month',
    date_labels = "%b"
  ) +
  scale_fill_manual(
    name = '',
    values = c('Reported cases' = 'gray76')
  ) +
  ylab('Cases') +
  xlab('Date') +
  ggtitle('Reported and estimated COVID-19 cases in European Countries') +
  ptheme2


pdf('Figures/Figure-2.pdf', width = 8-1/4, height = 11-3/4)
p
dev.off()


############################################
#### Figure 3 (Summary of EWS results) #####
#############################################
bw <- 4
ws <- 25
d <- read.csv('Results/ews-results.csv') %>% 
  filter(
    windowsize == ws, bandwidth == bw, backward_only == TRUE
  ) %>% 
  select(-windowsize, -bandwidth)

dp <- d %>% 
  select(country_name, ends_with('pvalue')) %>%
  melt() %>% 
  rename(ews_pvalue = variable, pvalue = value)

dk <- d %>% 
  select(ends_with('tau')) %>% 
  melt() %>% 
  rename(ews_value = variable, tau = value)

alpha <- 0.05
dj <- cbind(dp, dk) %>% 
  mutate(
    is_significant = pvalue < alpha | (1 - pvalue) < alpha,
    color = ifelse(is.na(is_significant), 'black', ifelse(is_significant, cols[1], 'black')),
    country_label = ifelse(is_significant, country_name, '')
  )

dj$ews <- factor(
  dj$ews_value,
  levels = rev(paste0(
    c(
    'mean', 'variance', 'coefficient_of_variation', 
    'index_of_dispersion', 'autocovariance', 'autocorrelation', 'decay_time',
    'skewness', 'kurtosis', 'variance_first_diff'
    ), '_tau'
  )),
  labels = rev(c(
    'Mean', 'Variance', 'Coef. of variation', 'Dispersion index', 'Autocovariance',
    'Autocorrelation', 'Decay time', 'Skewness', 'Kurtosis', 'Variance (first diff.)'
  )))


pdf('Figures/Figure-3.pdf', width = 9-1/4, height = 9-3/4)
psummary <- ggplot(
  dj, aes(x = factor(ews), y = tau, label = country_label)
) +
  geom_jitter(position = position_jitter(seed = 1, width = 0.2), aes(color = color)) +
  geom_text_repel(
    position = position_jitter(seed = 1, width = 0.2),
    size = 2.6
  ) +
  scale_color_identity() +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25))+
  scale_x_discrete()+
  coord_flip() +
  xlab('') +
  ylab(expression('Kendall\'s ' ~ tau)) +
  geom_segment(aes(y = -1, yend = 1, x = -Inf, xend = -Inf))+
  geom_segment(aes(x = 1, xend = 10, y = -Inf, yend = -Inf)) +
  ggtitle(expression('Kendall\'s ' ~ tau ~ ' for early warning indicators across countries')) +
  ptheme2
psummary
dev.off()


####################################
# Create plots for Appendix (A9-A12)
####################################
countries <- sort(unique(df$country))

len <- length(countries)
ylims <- df %>% 
  filter(country %in% countries) %>% 
  group_by(country) %>% 
  summarize(max = max(infections_mean)) %>% 
  select(max)

ylims <- ceiling(ylims$max * 1.5)

grob_tab <- c()
rt_plots <- list()
cases_plots <- c()

#' Plot country and R(t)
get_plot_grob <- function(
  country, up_ylim, xlab = NULL,
  ylab1 = NULL, ylab2 = NULL, margin = NULL) {
  
  p <- plot_dat(
    df[df$country == country, ],
    ylims = c(0, up_ylim), Rlimits = c(0, 3),
    date_breaks = '1 months', title = country,
    ylab1 = ifelse(is.null(ylab1), '', 'Confirmed cases'),
    ylab2 = ifelse(is.null(ylab2), '', expression(R[t])),
    xlab = ifelse(is.null(xlab), '', 'Date')
  )
  
  if (!is.null(margin)) {
    gA <- ggplotGrob(p$p1 + margin + theme(axis.title.y = element_text(size = 12)))
    gB <- ggplotGrob(p$p2 + margin + theme(axis.title.y = element_text(size = 12)))
  } else {
    gA <- ggplotGrob(p$p1 + theme(axis.title.y = element_text(size = 12)))
    gB <- ggplotGrob(p$p2 + theme(axis.title.y = element_text(size = 12)))
  }
  
  rbind(rbind(gA, gB))
}

margin1 <- theme(plot.margin = unit(c(0.2, 0, 0.2, 0.2), 'cm'))
margin2 <- theme(plot.margin = unit(c(0.2, 0.75, 0.2, 0), 'cm'))

rt_plots1 <- rbind(
  get_plot_grob(countries[1], ylims[1], ylab1 = T, ylab2 = T, margin = margin1),
  get_plot_grob(countries[3], ylims[3], ylab1 = T, ylab2 = T, margin = margin1),
  get_plot_grob(countries[5], ylims[5], ylab1 = T, ylab2 = T, xlab = T, margin = margin1)
)

rt_plots2 <- rbind(
  get_plot_grob(countries[2], ylims[2], margin = margin2),
  get_plot_grob(countries[4], ylims[4], margin = margin2),
  get_plot_grob(countries[6], ylims[6], xlab = 'Date', margin = margin2)
)

rt_plots3 <- rbind(
  get_plot_grob(countries[7], ylims[7], ylab1 = T, ylab2 = T, margin = margin1),
  get_plot_grob(countries[9], ylims[9], ylab1 = T, ylab2 = T, margin = margin1),
  get_plot_grob(countries[11], ylims[11], ylab1 = T, ylab2 = T, xlab = T, margin = margin1)
)

rt_plots4 <- rbind(
  get_plot_grob(countries[8], ylims[8], margin = margin2),
  get_plot_grob(countries[10], ylims[10], margin = margin2),
  get_plot_grob(countries[12], ylims[12], xlab = 'Date', margin = margin2)
)

rt_plots5 <- rbind(
  get_plot_grob(countries[13], ylims[13], ylab1 = T, ylab2 = T, margin = margin1),
  get_plot_grob(countries[15], ylims[15], ylab1 = T, ylab2 = T, margin = margin1),
  get_plot_grob(countries[17], ylims[17], ylab1 = T, ylab2 = T, xlab = T, margin = margin1)
)

rt_plots6 <- rbind(
  get_plot_grob(countries[14], ylims[14], margin = margin2),
  get_plot_grob(countries[16], ylims[16], margin = margin2),
  get_plot_grob(countries[18], ylims[18], xlab = 'Date', margin = margin2)
)

rt_plots7 <- rbind(
  get_plot_grob(countries[19], ylims[19], ylab1 = T, ylab2 = T, margin = margin1),
  get_plot_grob(countries[21], ylims[21], ylab1 = T, ylab2 = T, margin = margin1),
  get_plot_grob(countries[23], ylims[23], ylab1 = T, ylab2 = T, xlab = T, margin = margin1)
)

rt_plots8 <- rbind(
  get_plot_grob(countries[20], ylims[20], margin = margin2),
  get_plot_grob(countries[22], ylims[22], margin = margin2),
  get_plot_grob(countries[24], ylims[24], xlab = 'Date', margin = margin2)
)

margin3 <- theme(plot.margin = unit(c(4.25, 0.5, 4.25, 0.2), 'cm'))
gA <- get_plot_grob(countries[25], ylims[25], ylab1 = T, ylab2 = T, margin = margin1)
gB <- get_plot_grob(countries[26], ylims[26], xlab = 'Date', ylab1 = T, ylab2 = T, margin = margin1)
gC <- get_plot_grob(countries[27], ylims[27], xlab = 'Date', margin = margin2)

rt_plots1 <- cbind(rt_plots1, rt_plots2)
rt_plots2 <- cbind(rt_plots3, rt_plots4)
rt_plots3 <- cbind(rt_plots5, rt_plots6)
rt_plots4 <- cbind(rt_plots7, rt_plots8)
rt_plots5 <- cbind(rt_plots7, rt_plots8)

pdf('Figures/Figure-A08.pdf', width = 8-1/4, height = 11-3/4)
grid::grid.newpage(); grid::grid.draw(rt_plots1)
dev.off()

pdf('Figures/Figure-A09.pdf', width = 8-1/4, height = 11-3/4)
grid::grid.newpage(); grid::grid.draw(rt_plots2)
dev.off()

pdf('Figures/Figure-A10.pdf', width = 8-1/4, height = 11-3/4)
grid::grid.newpage(); grid::grid.draw(rt_plots3)
dev.off()

pdf('Figures/Figure-A11.pdf', width = 8-1/4, height = 11-3/4)
grid::grid.newpage(); grid::grid.draw(rt_plots4)
dev.off()

pdf('Figures/Figure-A12.pdf', width = 8-1/4, height = 8-3/4)
ggarrange(rbind(gA, gB), gC)
dev.off()


####################################
# Create plots for Appendix (A13-A22)
####################################
add_sensitivity <- function(p) {
  p +
    geom_tile() +
    # uncomment when plotting the mean
    # geom_point() +
    # scale_color_gradientn(
    scale_fill_gradientn(
      name = 'p-value',
      colors = brewer.pal(8, 'RdYlBu'), limits = c(0, 0.25),
      trans = 'sqrt'
    ) +
    xlab('Estimation window size') +
    ylab('Detrending window size') +
    facet_wrap(~ country_name, ncol = 4) +
    guides(
      fill = guide_colorbar(
        title.position = 'top', title.hjust = 0.5,
        title.size = 10, barheight = 10, barwidth = 0.75
      ),
      color = guide_colorbar(
        title.position = 'top', title.hjust = 0.5,
        title.size = 10, barheight = 10, barwidth = 0.75
      )
    ) +
    ptheme2 + 
    theme(
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      legend.position = 'right'
    )
}

# Load EWS results
df_ews <- read.csv('Results/ews-results.csv')

# Assess which EWS showed a significant increase
ews <- df_ews %>%
  mutate(
    autocorrelation_pvalue_censored = ifelse(autocorrelation_pvalue > 0.25, 0.25, autocorrelation_pvalue),
    variance_pvalue_censored = ifelse(variance_pvalue > 0.25, 0.25, variance_pvalue),
    variance_first_diff_pvalue_censored = ifelse(variance_first_diff_pvalue > 0.25, 0.25, variance_first_diff_pvalue),
    coefficient_of_variation_pvalue_censored = ifelse(coefficient_of_variation_pvalue > 0.25, 0.25, coefficient_of_variation_pvalue),
    decay_time_pvalue_censored = ifelse(decay_time_pvalue > 0.25, 0.25, decay_time_pvalue),
    
    mean_pvalue_censored = ifelse(mean_pvalue > 0.25, 0.25, mean_pvalue_abs),
    kurtosis_pvalue_censored = ifelse(kurtosis_pvalue > 0.25, 0.25, kurtosis_pvalue),
    skewness_pvalue_censored = ifelse(skewness_pvalue > 0.25, 0.25, skewness_pvalue),
    autocovariance_pvalue_censored = ifelse(autocovariance_pvalue > 0.25, 0.25, autocovariance_pvalue),
    index_of_dispersion_pvalue_censored = ifelse(index_of_dispersion_pvalue > 0.25, 0.25, index_of_dispersion_pvalue)
  )

pdf('Figures/Figure-A13.pdf', width = 8-1/4, height = 11-3/4)
p13 <- ggplot(ews, aes(x = factor(bandwidth), y = mean_pvalue_censored, color = mean_pvalue_censored)) +
  ggtitle('Sensitivity analysis for the mean')
add_sensitivity(p13) + ylab(expression(p ~ 'value')) + xlab('Detrending window size') + ylim(c(0, 0.25))
dev.off()

pdf('Figures/Figure-A14.pdf', width = 8-1/4, height = 11-3/4)
p14 <- ggplot(ews, aes(x = factor(windowsize), y = factor(bandwidth), fill = variance_pvalue_censored)) +
  ggtitle('Sensitivity analysis for the variance')
add_sensitivity(p14)
dev.off()

pdf('Figures/Figure-A15.pdf', width = 8-1/4, height = 11-3/4)
p15 <- ggplot(ews, aes(x = factor(windowsize), y = factor(bandwidth), fill = coefficient_of_variation_pvalue_censored)) +
  ggtitle('Sensitivity analysis for the coefficient of variation')
add_sensitivity(p15)
dev.off()

pdf('Figures/Figure-A16.pdf', width = 8-1/4, height = 11-3/4)
p16 <- ggplot(ews, aes(x = factor(windowsize), y = factor(bandwidth), fill = index_of_dispersion_pvalue_censored)) +
  ggtitle('Sensitivity analysis for the index of dispersion')
add_sensitivity(p16)
dev.off()

pdf('Figures/Figure-A17.pdf', width = 8-1/4, height = 11-3/4)
p17 <- ggplot(ews, aes(x = factor(windowsize), y = factor(bandwidth), fill = autocovariance_pvalue_censored)) +
  ggtitle('Sensitivity analysis for the autocovariance')
add_sensitivity(p17)
dev.off()

pdf('Figures/Figure-A18.pdf', width = 8-1/4, height = 11-3/4)
p18 <- ggplot(
  ews,
  aes(x = factor(windowsize), y = factor(bandwidth), fill = autocorrelation_pvalue_censored)
) + ggtitle('Sensitivity analysis for the autocorrelation')
add_sensitivity(p18)
dev.off()

pdf('Figures/Figure-A19.pdf', width = 8-1/4, height = 11-3/4)
p19 <- ggplot(ews, aes(x = factor(windowsize), y = factor(bandwidth), fill = decay_time_pvalue_censored)) +
  ggtitle('Sensitivity analysis for the decay time')
add_sensitivity(p19)
dev.off()

pdf('Figures/Figure-A20.pdf', width = 8-1/4, height = 11-3/4)
p20 <- ggplot(ews, aes(x = factor(windowsize), y = factor(bandwidth), fill = skewness_pvalue_censored)) +
  ggtitle('Sensitivity analysis for the skewness')
add_sensitivity(p20)
dev.off()

pdf('Figures/Figure-A21.pdf', width = 8-1/4, height = 11-3/4)
p21 <- ggplot(ews, aes(x = factor(windowsize), y = factor(bandwidth), fill = kurtosis_pvalue_censored)) +
  ggtitle('Sensitivity analysis for the kurtosis')
add_sensitivity(p21)
dev.off()

pdf('Figures/Figure-A22.pdf', width = 8-1/4, height = 11-3/4)
p22 <- ggplot(ews, aes(x = factor(windowsize), y = factor(bandwidth), fill = variance_first_diff_pvalue_censored)) +
  ggtitle('Sensitivity analysis for the first difference in variance')
add_sensitivity(p22)
dev.off()
