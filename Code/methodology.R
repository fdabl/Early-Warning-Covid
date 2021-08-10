source('Code/helpers.R')
cols <- brewer.pal(3, 'Set1')


#################################
#### Figure 1 (Methodology) #####
#################################

# Custom R(t) function
f <- function(x) {

  y <- c(2, 0.70, 0.60, 0.61, 0.70, 0.80, 1, 1.3, 1.2, 1)
  points <- data.frame(
    y = y, x = seq(0, length(y) - 1)
  )

  m <- lm(
    y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5) + I(x^6) + I(x^7) + I(x^8) + I(x^9),
    data = points
  )
  b <- coef(m)

  mat <- cbind(x^0, x^1, x^2, x^3, x^4, x^5, x^6, x^7, x^8, x^9)
  as.numeric(mat %*% b)
}

times <- seq(1, 400, by = 1)
x <- seq(0, 8, length.out = length(times))
Rt <- f(x)


N <- 100000
S0 <- 0.95 * N
I0 <- 0.05 * N
eta <- 5e-5
eta_start <- 100
sigma <- 1 / 5.2
gamma <- 1 / 10
beta_par_t <- Rt * gamma
eta_par_t <- c(rep(0, eta_start), rep(eta, length(times) - eta_start))

set.seed(1)
sim <- create_SEIR(
  times, beta_par_t = beta_par_t, eta_par_t = eta_par_t,
  params = c(sigma = 1 / 5.2, gamma = gamma, eta = eta, rho = 1),
  init = c(S0 = S0, E0 = 0, I0 = I0, R0 = 0, N = N)
)

simres <- simulate(sim, nsim = 1, format = 'data.frame')[-1, ]
simres$reports[1] <- 400
ix <- seq(1, nrow(simres) - 1, 5)
simres <- simres[ix, ]
Rtix <- Rt[ix]

pdf('Figures/Figure-1.pdf', width = 10-1/4, height = 8-1/4)

cex.lab <- 1.25
cex.axis <- 1.25
cex.text <- 1.25
colcases <- 'gray76'
layout.matrix <- rbind(c(1, 1), c(2, 3))
layout(mat = layout.matrix, heights = c(2, 2), widths = c(2, 2))

par(mar = c(5.5, 12, 3, 9))

barplot(
  simres$reports, pch = 20, axes = FALSE,
  xlab = '', ylab = 'Reported cases',
  ylim = c(0, 800), xlim = c(0, 80),
  col = colcases, xaxs = 'i', yaxs = 'i',
  cex.lab = cex.lab, space = 0
)
axis(1, cex.axis = cex.axis, at = seq(0, 80, 10))
axis(2, las = 2, cex.axis = cex.axis)
ix_min <- which.min(Rtix)
ix_max <- which.max(Rtix[-seq(20)]) + 20

par(new = TRUE)
plot(
  Rtix, type = 'l', col = 'black', lwd = 2, axes = FALSE,
  xlab = '', ylab = '', ylim = c(0, 3),
  xaxs = 'i', yaxs = 'i', cex.lab = cex.lab
)
axis(4, labels = seq(0, 3), at = seq(0, 3), cex.axis = cex.axis)
mtext(TeX('$R_t$'), side = 4, line = 1.5, cex = 1)
mtext('Time', side = 1, line = 2.5, cex = 1)
mtext('Illustration of the Methodology', side = 3, line = 1.25, cex = 1.5)
lines(c(0, 80), c(1, 1), lty = 'dashed', col = cols[2])
lines(c(ix_min, ix_min), c(0, 800), lwd = 2, col = cols[2])
lines(c(ix_max, ix_max), c(0, 800), lwd = 2, col = cols[2])

size <- 15
height <- 2
arrows(ix_max - 2 - size, height, ix_max - 2, height, col = 'black', length = 0.15, lwd = 1.5)
arrows(ix_min + 2 + size, height, ix_min + 2, height, col = 'black', length = 0.15, lwd = 1.5)
text(50, height, labels = 'Time period', cex = cex.text)

legend(
  26, 3,
  legend = c(expression(R[t]), 'Reported cases'),
  pch = c(NA, 15), lty = c(1, NA), col = c('black', 'gray66'), cex = 1, bty = 'n', lwd = 2
)


text(11.7, 2.8, labels = '(a) Select time period', cex = cex.text)

par(mar = c(5.1, 4.5, 0, 1.5))
tt <- seq(ix_min, ix_max)
y <- simres$reports[seq(ix_min, ix_max)]

stats <- get_stats(
  y,
  center_trend = 'local_constant',
  center_kernel = 'uniform',
  center_bandwidth = 3,
  stat_trend = 'local_constant',
  stat_kernel = 'uniform',
  stat_bandwidth = 10,
  backward_only = TRUE
)

ymean <- stats$stats$mean
yvariance <- stats$stats$variance

barplot(
  y, axes = FALSE,
  xlab = '', ylab = 'Reported cases',
  ylim = c(0, 80),
  col = colcases, xaxs = 'i', yaxs = 'i',
  cex.lab = cex.lab, space = 0
)

lines(ymean, pch = 20, col = 'black', lwd = 2)
axis(1, cex.axis = cex.axis, labels = seq(ix_min, ix_max, 10), at = seq(ix_min, ix_max, 10) - 25)
axis(2, las = 2, cex.axis = cex.axis)
mtext('Time', side = 1, line = 2.5, cex = 1)

text(
  13.9, 75, cex = cex.text,
  labels = TeX('(b) Detrending: $\\bar{y}_t = \\frac{1}{\\delta_1} \\, \\sum_{t - \\delta_1}^{t} y_j$')
)

legend(
  'topright',
  legend = c('Mean', 'Reported cases'),
  pch = c(NA, 15), lty = c(1, NA), col = c('black', 'gray66'), cex = 1, bty = 'n', lwd = 2
)

par(mar = c(5.1, 3.9, 0, 2.1))
ydetrended <- y - ymean
barplot(
  ydetrended, pch = 20, axes = FALSE,
  xlab = '', ylab = 'Detrended cases',
  ylim = c(-20, 30),
  col = colcases, xaxs = 'i', yaxs = 'i',
  cex.lab = cex.lab, space = 0
)
v <- stats$stats$autocorrelation
v <- stats$stats$variance
aut <- (v - min(v, na.rm = TRUE)) / (max(v, na.rm = TRUE) - min(v, na.rm = TRUE))
lines(aut * 20, lwd = 2, col = 'black')
axis(1, cex.axis = cex.axis, labels = seq(ix_min, ix_max, 10), at = seq(ix_min, ix_max, 10) - 25)
axis(2, las = 2, cex.axis = cex.axis)
mtext('Time', side = 1, line = 2.5, cex = 1)
text(
  20.25, 26, cex = cex.text,
  labels = TeX('(c) Estimation: $s_t = f(y_t, \\ldots, y_{t - \\delta_2}, \\bar{y}_t, \\ldots, \\bar{y}_{t - \\delta_2})$')
)
legend(
  'bottomright',
  legend = c('Variance (scaled)', 'Reported cases (detrended)'),
  pch = c(NA, 15), lty = c(1, NA), col = c('black', 'gray66'), cex = 1, bty = 'n', lwd = 2
)

dev.off()