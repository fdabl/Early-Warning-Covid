source('Code/helpers.R')


# Load or create the relevant data
# (Not provided in the repository because it is too big)
if (!file.exists('Data/European-Rt-Data.RDS')) {
  source('Code/setup-data.R')
  
} else {
  
  df <- read.csv('Data/cases-and-rt-data-criterion.csv') %>%
    mutate(
      date = as.POSIXct(date),
      first_low_date = as.POSIXct(first_low_date),
      second_max_date = as.POSIXct(second_max_date)
    )
}

         
# Estimate early warning indicators across bandwidths (i.e., detrending rolling window size)
# and estimation rolling window sizes
registerDoParallel(cores = 10)
bw <- seq(2, 20, 2)
ws <- seq(5, 50, 5)
df_ews <- get_early_warning(df, ws, bw, nr_surrogates = 500, backward_only = TRUE)

# Save results
write.csv(df_ews, 'Results/ews-results.csv', row.names = FALSE)


####################
#### Create Table 1
####################
dc <- df %>%
  group_by(
    country
  ) %>% 
  mutate(
    Rlen = second_max_ix - first_low_ix,
    Rmin = round(R_mean[first_low_ix], 2),
    Rmax = round(R_mean[second_max_ix], 2),
    Rforcing = round(Rmax - Rmin, 2)
  ) %>% 
  select(country, Rlen, Rmin, Rmax, Rforcing) %>% 
  distinct()


df_ews <- read.csv('Results/ews-results.csv') %>%
  filter(
    windowsize == 25, bandwidth == 4, backward_only == TRUE
  ) %>% 
  select(-windowsize, -bandwidth) %>% 
  rename(country = country_name) %>% 
  left_join(dc)

dp <- df_ews %>% 
  select(country, ends_with('pvalue')) %>% 
  melt() %>% 
  rename(ews_pvalue = variable, pvalue = value)

dk <- df_ews %>% 
  select(ends_with('tau')) %>% 
  melt() %>% 
  rename(ews_value = variable, tau = value)

alpha <- 0.05
dj <- cbind(dp, dk) %>% 
  mutate(
    is_significant = pvalue < alpha | (1 - pvalue) < alpha
  ) %>% 
  left_join(dc)


# Make a ranking of countries
ranking <- dj %>% 
  group_by(country) %>% 
  summarize(
    sig_large = sum(pvalue < alpha, na.rm = TRUE),
    sig_small = sum((1 - pvalue) < alpha, na.rm = TRUE)
  ) %>% 
  arrange(desc(sig_large), desc(sig_small)) %>% 
  data.frame

dd <- left_join(ranking, dc)

tab <- dd[, -7]
colnames(tab) <- c(
  'Country', 'Significant increase', 'Significant decrease',
  'Duration', expression(R[min]), expression(R[max])
)

tab
