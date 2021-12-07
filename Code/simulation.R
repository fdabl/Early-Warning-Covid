source('Code/helpers.R')
cols <- brewer.pal(3, 'Set1')


rho <- 1
eta <- 2e-05
sigma <- 1/5.2
gamma <- 1/10
eta_start <- 50

nr_decrease <- 25
nr_constant <- c(25, 50, 100, 200)
nr_increase <- seq(25, 200, 5)

if (!file.exists('Results/simulation-results-AUC.RDS')) {
  
  # Simulation study comparing second waves to no second waves scenarios (AUC)
  registerDoParallel(cores = 9)
  res_auc <- simulate_ews(
    eta = eta,
    bw = 'dynamic', ws = 25, type = 'uniform',
    lowest_point = 0.50,
    nr_decrease = nr_decrease,
    nr_constant = nr_constant,
    nr_increase = nr_increase,
    eta_start = 50,
    nsims = 500, detrending = TRUE,
    cut_first_window = TRUE
  )
  
  saveRDS(res_auc, 'Results/simulation-results-AUC.RDS')
  
} else {
  
  res_auc <- readRDS('Results/simulation-results-AUC.RDS')
}


if (!file.exists('Results/simulation-results-TPR.RDS')) {
  
  # Simulation study using only second wave scenarios (ARMA / true positive rate)
  registerDoParallel(cores = 9)
  res_tpr <- simulate_ews(
    eta = eta,
    bw = 'dynamic', ws = 25, type = 'uniform',
    lowest_point = 0.50,
    nr_decrease = nr_decrease,
    nr_constant = nr_constant,
    nr_increase = nr_increase,
    eta_start = 50,
    nsims = 250, detrending = TRUE,
    nr_surrogates = 500, fit_ARMA = TRUE,
    cut_first_window = TRUE
  )
  
  saveRDS(res_tpr, 'Results/simulation-results-TPR.RDS')
  
} else {
  
  res_tpr <- readRDS('Results/simulation-results-TPR.RDS')
}


#' Compute AUC / TPR from simulation results
compute_metrics <- function(res, ews_names, AUC = FALSE) {

  outcome <- compute_AUC(res, ews_names, AUC = AUC, alpha = 0.05)
  outcome <- melt(outcome, c('nr_decrease', 'nr_constant', 'nr_increase', 'ylen_mean', 'ylen_sd'))
  
  outcome$ews <- factor(
    outcome$variable,
    levels = rev(ews_names),
    labels = rev(c(
      'Mean', 'Variance', 'Coef. of variation', 'Dispersion index', 'Autocovariance',
      'Autocorrelation', 'Decay time', 'Skewness', 'Kurtosis', 'Variance (first diff.)'
    )))
  
  outcome
}

ews_names_auc <- c(
  'mean', 'variance', 'coefficient_of_variation', 'index_of_dispersion',
  'autocovariance', 'autocorrelation', 'decay_time', 'skewness',
  'kurtosis', 'variance_first_diff'
)
ews_names_tpr <- paste0(ews_names_auc, '_pval')
auc <- compute_metrics(res_auc, ews_names_auc, AUC = TRUE)
tpr <- compute_metrics(res_tpr, ews_names_tpr, AUC = FALSE)


# Get nice color palette
pals <- rev(brewer.pal(10, 'RdYlBu'))
pals2 <- c(pals[seq(5)], 'gray75', pals[seq(6, 10)])
cc <- brewer.pal(5, 'Reds')
gg <- unikn::pal_grau
colss <- c(rev(gg), 'white', cc)


# Create AUC figure
pauc <- create_AUC_plot(
  na.omit(auc), legend = TRUE, ytext_size = 12, ylab = '', title = 'Area under the Curve',
  colpallete = pals2
) + scale_x_discrete(breaks = c(30, seq(50, 200, 25)))


pdf('Figures/Figure-5a-raw.pdf', width = 8-1/4, height = 11-3/4)
pauc
dev.off()


# Create TPR figure
ptpr <- create_AUC_plot(
  tpr, legend = TRUE, ytext_size = 12, ylab = '', title = 'True positive rate',
  legend_name = 'TPR', colpallete = pals2
) + scale_x_discrete(breaks = c(30, seq(50, 200, 25)))

pdf('Figures/Figure-5b-raw.pdf', width = 8-1/4, height = 11-3/4)
ptpr
dev.off()
