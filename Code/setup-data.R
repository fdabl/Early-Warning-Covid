source('Code/helpers.R')

# EU without Spain but with UK
EU <- c(
  'Austria','Belgium','Bulgaria','Croatia','Cyprus',
  'Czechia','Denmark','Estonia','Finland','France',
  'Germany','Greece','Hungary','Ireland','Italy','Latvia',
  'Lithuania','Luxembourg','Malta','Netherlands','Poland',
  'Portugal','Romania','Slovakia','Slovenia',
  'Sweden','United Kingdom' # 'Spain'
)


# Get data from the WHO
dat <- get_national_data(EU, source = 'who')
  mutate(
    date = as.Date(date)
  ) %>%
  filter(
    date < '2020-09-30'
  ) %>% as_tibble


if (!file.exists('Data/European-Rt-Data.RDS')) {
  
  # Estimate R(t) for European countries
  res <- list()
  for (i in seq(length(EU))) {
    tryCatch(
      res[EU[i]] <- estimate_rt(dat, EU[i]),
      error = function(e) {
        print(e)
      })
    
    saveRDS(res, 'Data/European-Rt-Data.RDS')
  }
} else {
  
  res <- readRDS('Data/European-Rt-Data.RDS')
}


# List of offsets to help find the minimum and maximum R(t)
# values in between the first and the second wave for European countries
offsets <- list(
  'Germany' = c('2020-03-01', NULL),
  'Italy' = c('2020-05-01', NULL),
  'Lithuania' = c('2020-06-01', NULL),
  'Malta' = c('2020-06-01', NULL),
  'Netherlands' = c('2020-06-01', NULL),
  'United Kingdom' = c('2020-06-01', NULL),
  'Belgium' = c('2020-05-01', NULL),
  'Poland' = c(NULL, '2020-08-01'),
  'Cyprus' = c(NULL, '2020-08-01')
)


# Create the start and end dates for the time-series between
# the first and second wave in European countries
dat_list <- lapply(EU, function(name) {
  
  start_date <- NULL
  end_date <- NULL
  
  if (name %in% names(offsets)) {
    if (name %in% c('Poland', 'Cyprus')) {
      end_date <- offsets[[name]][1]
    } else {
      start_date <- offsets[[name]][1]
    }
  }
  
  print(paste(name, '-', 'start:', start_date, '- end:', end_date))
  
  create_dat(
    res[[name]], country = name,
    cutoff = 0.99, start_date = start_date, end_date = end_date
  )
})


# Save these data, which combine reported cases, R(t), and the start and end times to
# select the time-series between the first and second waves in European countries
df <- do.call('rbind', dat_list)
write.csv(df, 'Data/cases-and-rt-data-criterion.csv', row.names = FALSE)
