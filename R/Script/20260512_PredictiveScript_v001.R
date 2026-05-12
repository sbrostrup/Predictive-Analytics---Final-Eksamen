# ============================================================
# 0. Setup
# ============================================================

rm(list = ls())

library(tidyverse)
library(lubridate)
library(knitr)
library(ggplot2)
library(patchwork)
library(fpp3)
library(tseries)
library(feasts)
library(zoo)

theme_set(theme_minimal())


# ============================================================
# 1. Load raw data from FRED
# ============================================================

read_fred <- function(fred_code, variable_name) {
  url <- paste0(
    "https://fred.stlouisfed.org/graph/fredgraph.csv?id=",
    fred_code
  )
  
  read_csv(url, na = c(".", ""), show_col_types = FALSE) %>%
    rename(
      date = observation_date,
      !!variable_name := all_of(fred_code)
    ) %>%
    mutate(date = as.Date(date))
}

dk_unemployment_raw <- read_fred("LRHUADTTDKM156N", "dk_unemployment")
ez_unemployment_raw <- read_fred("LRHUADTTEZM156N", "ez_unemployment")
dk_hicp_raw <- read_fred("CP0000DKM086NEST", "dk_hicp")
ez_hicp_raw <- read_fred("CP0000EZ19M086NEST", "ez_hicp")


# ============================================================
# 2. Inspect raw data
# ============================================================

# The unemployment variables are already rates.
# The HICP variables are price index levels.
# We inspect the raw data before transforming HICP into inflation.

p_raw1 <- dk_unemployment_raw %>%
  ggplot(aes(x = date, y = dk_unemployment)) +
  geom_line() +
  labs(
    title = "Raw data: Danish unemployment",
    x = "Date",
    y = "Unemployment rate"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

p_raw2 <- ez_unemployment_raw %>%
  ggplot(aes(x = date, y = ez_unemployment)) +
  geom_line() +
  labs(
    title = "Raw data: Euro-area unemployment",
    x = "Date",
    y = "Unemployment rate"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

p_raw3 <- dk_hicp_raw %>%
  ggplot(aes(x = date, y = dk_hicp)) +
  geom_line() +
  labs(
    title = "Raw data: Danish HICP index",
    x = "Date",
    y = "HICP index"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

p_raw4 <- ez_hicp_raw %>%
  ggplot(aes(x = date, y = ez_hicp)) +
  geom_line() +
  labs(
    title = "Raw data: Euro-area HICP index",
    x = "Date",
    y = "HICP index"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

(p_raw1 + p_raw2) / (p_raw3 + p_raw4) +
  plot_annotation(
    title = "Inspection of raw data"
  )


# ============================================================
# 3. Transform HICP index into year-on-year inflation
# ============================================================

# The HICP series are price index levels, not inflation rates.
# Since the Phillips curve concerns inflation, we transform HICP into
# year-on-year inflation:
#
# inflation_t = 100 * (HICP_t / HICP_{t-12} - 1)
#
# This compares prices with the same month in the previous year and therefore
# gives an annual inflation rate.
#
# Unemployment is already measured as a rate, so it is kept unchanged.

dk_inflation <- dk_hicp_raw %>%
  arrange(date) %>%
  mutate(
    dk_inflation = 100 * (dk_hicp / lag(dk_hicp, 12) - 1)
  ) %>%
  drop_na(dk_inflation)

ez_inflation <- ez_hicp_raw %>%
  arrange(date) %>%
  mutate(
    ez_inflation = 100 * (ez_hicp / lag(ez_hicp, 12) - 1)
  ) %>%
  drop_na(ez_inflation)

dk_unemployment <- dk_unemployment_raw
ez_unemployment <- ez_unemployment_raw


# ============================================================
# 4. Add calendar variables for plotting
# ============================================================

dk_unemployment <- dk_unemployment %>%
  mutate(
    year = year(date),
    month = month(date),
    month_name = month(date, label = TRUE, abbr = TRUE)
  )

ez_unemployment <- ez_unemployment %>%
  mutate(
    year = year(date),
    month = month(date),
    month_name = month(date, label = TRUE, abbr = TRUE)
  )

dk_inflation <- dk_inflation %>%
  mutate(
    year = year(date),
    month = month(date),
    month_name = month(date, label = TRUE, abbr = TRUE)
  )

ez_inflation <- ez_inflation %>%
  mutate(
    year = year(date),
    month = month(date),
    month_name = month(date, label = TRUE, abbr = TRUE)
  )


# ============================================================
# 5. Inspect transformed data: trend
# ============================================================

# We now inspect trend in the variables that will be used in the analysis:
# unemployment rates and year-on-year inflation rates.

p1_trend <- dk_unemployment %>%
  ggplot(aes(x = date, y = dk_unemployment)) +
  geom_line() +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "red") +
  labs(
    title = "Danish unemployment rate",
    x = "Date",
    y = "Unemployment rate"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

p2_trend <- ez_unemployment %>%
  ggplot(aes(x = date, y = ez_unemployment)) +
  geom_line() +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "red") +
  labs(
    title = "Euro-area unemployment rate",
    x = "Date",
    y = "Unemployment rate"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

p3_trend <- dk_inflation %>%
  ggplot(aes(x = date, y = dk_inflation)) +
  geom_line() +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "red") +
  labs(
    title = "Danish inflation",
    x = "Date",
    y = "Inflation, % y/y"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

p4_trend <- ez_inflation %>%
  ggplot(aes(x = date, y = ez_inflation)) +
  geom_line() +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "red") +
  labs(
    title = "Euro-area inflation",
    x = "Date",
    y = "Inflation, % y/y"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

p1_trend + p2_trend + p3_trend + p4_trend +
  plot_layout(nrow = 1) +
  plot_annotation(
    title = "Trends in unemployment and inflation"
  )


# ============================================================
# 6. Inspect seasonality: seasonal plots
# ============================================================

p1_season <- dk_unemployment %>%
  ggplot(aes(x = month, y = dk_unemployment, group = year)) +
  geom_line(alpha = 0.4) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  labs(
    title = "Seasonality: Danish unemployment",
    x = "Month",
    y = "Unemployment rate"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

p2_season <- ez_unemployment %>%
  ggplot(aes(x = month, y = ez_unemployment, group = year)) +
  geom_line(alpha = 0.4) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  labs(
    title = "Seasonality: Euro-area unemployment",
    x = "Month",
    y = "Unemployment rate"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

p3_season <- dk_inflation %>%
  ggplot(aes(x = month, y = dk_inflation, group = year)) +
  geom_line(alpha = 0.4) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  labs(
    title = "Seasonality: Danish inflation",
    x = "Month",
    y = "Inflation, % y/y"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

p4_season <- ez_inflation %>%
  ggplot(aes(x = month, y = ez_inflation, group = year)) +
  geom_line(alpha = 0.4) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  labs(
    title = "Seasonality: Euro-area inflation",
    x = "Month",
    y = "Inflation, % y/y"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

p1_season + p2_season + p3_season + p4_season +
  plot_layout(nrow = 1) +
  plot_annotation(
    title = "Seasonality in unemployment and inflation"
  )


# ============================================================
# 7. Inspect seasonality: seasonal subseries plots
# ============================================================

# Convert to tsibble format for gg_subseries().

dk_unemployment_ts <- dk_unemployment %>%
  mutate(date = yearmonth(date)) %>%
  as_tsibble(index = date)

ez_unemployment_ts <- ez_unemployment %>%
  mutate(date = yearmonth(date)) %>%
  as_tsibble(index = date)

dk_inflation_ts <- dk_inflation %>%
  mutate(date = yearmonth(date)) %>%
  as_tsibble(index = date)

ez_inflation_ts <- ez_inflation %>%
  mutate(date = yearmonth(date)) %>%
  as_tsibble(index = date)

p1_sub <- dk_unemployment_ts %>%
  gg_subseries(dk_unemployment) +
  labs(
    title = "Seasonal subseries: Danish unemployment",
    x = "Year",
    y = "Unemployment rate"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

p2_sub <- ez_unemployment_ts %>%
  gg_subseries(ez_unemployment) +
  labs(
    title = "Seasonal subseries: Euro-area unemployment",
    x = "Year",
    y = "Unemployment rate"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

p3_sub <- dk_inflation_ts %>%
  gg_subseries(dk_inflation) +
  labs(
    title = "Seasonal subseries: Danish inflation",
    x = "Year",
    y = "Inflation, % y/y"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

p4_sub <- ez_inflation_ts %>%
  gg_subseries(ez_inflation) +
  labs(
    title = "Seasonal subseries: Euro-area inflation",
    x = "Year",
    y = "Inflation, % y/y"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

(p1_sub + p2_sub) / (p3_sub + p4_sub) +
  plot_annotation(
    title = "Seasonal subseries plots: Denmark and Euro area"
  )


# ============================================================
# 8. Combine variables into one monthly tsibble
# ============================================================

macro_data <- dk_unemployment %>%
  dplyr::select(date, dk_unemployment) %>%
  left_join(
    ez_unemployment %>% dplyr::select(date, ez_unemployment),
    by = "date"
  ) %>%
  left_join(
    dk_inflation %>% dplyr::select(date, dk_inflation),
    by = "date"
  ) %>%
  left_join(
    ez_inflation %>% dplyr::select(date, ez_inflation),
    by = "date"
  ) %>%
  drop_na() %>%
  mutate(date = yearmonth(date)) %>%
  as_tsibble(index = date)


# ============================================================
# 9. Test stationarity in levels
# ============================================================

# ADF test:
# H0 = non-stationarity
#
# KPSS test:
# H0 = stationarity

# Danish unemployment
adf.test(macro_data$dk_unemployment)
kpss.test(macro_data$dk_unemployment)

# Euro-area unemployment
adf.test(macro_data$ez_unemployment)
kpss.test(macro_data$ez_unemployment)

# Danish inflation
adf.test(macro_data$dk_inflation)
kpss.test(macro_data$dk_inflation)

# Euro-area inflation
adf.test(macro_data$ez_inflation)
kpss.test(macro_data$ez_inflation)


# ============================================================
# 10. Visual confirmation for non-stationary data
# ============================================================

# Danish unemployment
macro_data %>%
  gg_tsdisplay(dk_unemployment, plot_type = "partial") +
  labs(title = "Time series display: Danish unemployment")

# Euro-area unemployment
macro_data %>%
  gg_tsdisplay(ez_unemployment, plot_type = "partial") +
  labs(title = "Time series display: Euro-area unemployment")

# Danish inflation
macro_data %>%
  gg_tsdisplay(dk_inflation, plot_type = "partial") +
  labs(title = "Time series display: Danish inflation")

# Euro-area inflation
macro_data %>%
  gg_tsdisplay(ez_inflation, plot_type = "partial") +
  labs(title = "Time series display: Euro-area inflation")


# ============================================================
# 11. First-order differencing for stationarity
# ============================================================

# First differencing is a statistical transformation used to obtain
# stationarity. This is separate from the economic transformation from
# HICP index to inflation rate.

macro_data <- macro_data %>%
  mutate(
    diff_dk_unemployment = difference(dk_unemployment),
    diff_ez_unemployment = difference(ez_unemployment),
    diff_dk_inflation = difference(dk_inflation),
    diff_ez_inflation = difference(ez_inflation)
  )

macro_data_diff <- macro_data %>%
  drop_na(
    diff_dk_unemployment,
    diff_ez_unemployment,
    diff_dk_inflation,
    diff_ez_inflation
  )


# ============================================================
# 12. Test stationarity after first differencing
# ============================================================

# Danish unemployment
adf.test(macro_data_diff$diff_dk_unemployment)
kpss.test(macro_data_diff$diff_dk_unemployment)

# Euro-area unemployment
adf.test(macro_data_diff$diff_ez_unemployment)
kpss.test(macro_data_diff$diff_ez_unemployment)

# Danish inflation
adf.test(macro_data_diff$diff_dk_inflation)
kpss.test(macro_data_diff$diff_dk_inflation)

# Euro-area inflation
adf.test(macro_data_diff$diff_ez_inflation)
kpss.test(macro_data_diff$diff_ez_inflation)


# ============================================================
# 13. Visual stationarity check after first differencing
# ============================================================

# Danish unemployment
macro_data_diff %>%
  gg_tsdisplay(diff_dk_unemployment, plot_type = "partial") +
  labs(title = "Time series display: Differenced Danish unemployment")

# Euro-area unemployment
macro_data_diff %>%
  gg_tsdisplay(diff_ez_unemployment, plot_type = "partial") +
  labs(title = "Time series display: Differenced Euro-area unemployment")

# Danish inflation
macro_data_diff %>%
  gg_tsdisplay(diff_dk_inflation, plot_type = "partial") +
  labs(title = "Time series display: Differenced Danish inflation")

# Euro-area inflation
macro_data_diff %>%
  gg_tsdisplay(diff_ez_inflation, plot_type = "partial") +
  labs(title = "Time series display: Differenced Euro-area inflation")


# ============================================================
# 14. Seasonal differencing diagnostic
# ============================================================

# Seasonal plots and ACF/PACF plots may suggest some annual dependence.
# We therefore inspect seasonal differencing as a diagnostic.
#
# Since seasonal differencing alone did not provide clearer stationarity
# evidence than first differencing in our tests, it is not used as the main
# transformation. Remaining seasonality can be handled later in the forecasting
# model through seasonal ARIMA terms, seasonal lags, or monthly dummies.

macro_data <- macro_data %>%
  mutate(
    seas_diff_dk_unemployment = difference(dk_unemployment, lag = 12),
    seas_diff_ez_unemployment = difference(ez_unemployment, lag = 12),
    seas_diff_dk_inflation = difference(dk_inflation, lag = 12),
    seas_diff_ez_inflation = difference(ez_inflation, lag = 12)
  )

macro_data_seas <- macro_data %>%
  drop_na(
    seas_diff_dk_unemployment,
    seas_diff_ez_unemployment,
    seas_diff_dk_inflation,
    seas_diff_ez_inflation
  )


# Stationarity tests after seasonal differencing

# Danish unemployment
adf.test(macro_data_seas$seas_diff_dk_unemployment)
kpss.test(macro_data_seas$seas_diff_dk_unemployment)

# Euro-area unemployment
adf.test(macro_data_seas$seas_diff_ez_unemployment)
kpss.test(macro_data_seas$seas_diff_ez_unemployment)

# Danish inflation
adf.test(macro_data_seas$seas_diff_dk_inflation)
kpss.test(macro_data_seas$seas_diff_dk_inflation)

# Euro-area inflation
adf.test(macro_data_seas$seas_diff_ez_inflation)
kpss.test(macro_data_seas$seas_diff_ez_inflation)


# ============================================================
# 15. Visual diagnostic after seasonal differencing
# ============================================================

# Danish unemployment
macro_data_seas %>%
  gg_tsdisplay(seas_diff_dk_unemployment, plot_type = "partial") +
  labs(title = "Time series display: Seasonally differenced Danish unemployment")

# Euro-area unemployment
macro_data_seas %>%
  gg_tsdisplay(seas_diff_ez_unemployment, plot_type = "partial") +
  labs(title = "Time series display: Seasonally differenced Euro-area unemployment")

# Danish inflation
macro_data_seas %>%
  gg_tsdisplay(seas_diff_dk_inflation, plot_type = "partial") +
  labs(title = "Time series display: Seasonally differenced Danish inflation")

# Euro-area inflation
macro_data_seas %>%
  gg_tsdisplay(seas_diff_ez_inflation, plot_type = "partial") +
  labs(title = "Time series display: Seasonally differenced Euro-area inflation")


# ============================================================
# 16. Optional: create zoo object for later GETS modelling
# ============================================================

# The gets package works with zoo objects. This object is created for later
# use in GETS / dynamic modelling and can be removed if not needed.

macro_zoo_data <- dk_unemployment_raw %>%
  dplyr::select(date, dk_unemployment) %>%
  left_join(
    ez_unemployment_raw %>% dplyr::select(date, ez_unemployment),
    by = "date"
  ) %>%
  left_join(
    dk_inflation %>% dplyr::select(date, dk_inflation),
    by = "date"
  ) %>%
  left_join(
    ez_inflation %>% dplyr::select(date, ez_inflation),
    by = "date"
  ) %>%
  drop_na()

macro_zoo <- zoo(
  macro_zoo_data %>%
    dplyr::select(
      dk_unemployment,
      ez_unemployment,
      dk_inflation,
      ez_inflation
    ),
  order.by = macro_zoo_data$date
)

# Inspect the zoo object
head(macro_zoo) %>%
  kable(digits = 3, align = "c")

# Plot the zoo object
plot(
  macro_zoo,
  main = "Denmark and Euro-area unemployment and inflation",
  col = "darkblue",
  lwd = 2
)

# ============================================================
# Extra diagnostic: seasonal + first differencing
# ============================================================

# Combined differencing:
# Delta Delta_12 y_t = (y_t - y_{t-12}) - (y_{t-1} - y_{t-13})
#
# This is a stronger transformation and is used only as a diagnostic.

macro_data <- macro_data %>%
  mutate(
    diff_seas_dk_unemployment = difference(seas_diff_dk_unemployment),
    diff_seas_ez_unemployment = difference(seas_diff_ez_unemployment),
    diff_seas_dk_inflation = difference(seas_diff_dk_inflation),
    diff_seas_ez_inflation = difference(seas_diff_ez_inflation)
  )

macro_data_diff_seas <- macro_data %>%
  drop_na(
    diff_seas_dk_unemployment,
    diff_seas_ez_unemployment,
    diff_seas_dk_inflation,
    diff_seas_ez_inflation
  )

# Stationarity tests after both seasonal and first differencing
adf.test(macro_data_diff_seas$diff_seas_dk_unemployment)
kpss.test(macro_data_diff_seas$diff_seas_dk_unemployment)

adf.test(macro_data_diff_seas$diff_seas_ez_unemployment)
kpss.test(macro_data_diff_seas$diff_seas_ez_unemployment)

adf.test(macro_data_diff_seas$diff_seas_dk_inflation)
kpss.test(macro_data_diff_seas$diff_seas_dk_inflation)

adf.test(macro_data_diff_seas$diff_seas_ez_inflation)
kpss.test(macro_data_diff_seas$diff_seas_ez_inflation)

# Visual diagnostics after both seasonal and first differencing
macro_data_diff_seas %>%
  gg_tsdisplay(diff_seas_dk_unemployment, plot_type = "partial") +
  labs(title = "Seasonal + first differenced Danish unemployment")

macro_data_diff_seas %>%
  gg_tsdisplay(diff_seas_dk_inflation, plot_type = "partial") +
  labs(title = "Seasonal + first differenced Danish inflation")




# ============================================================
# 18. Univariate ARIMA benchmark
# ============================================================

# The ARIMA benchmark is estimated on first-differenced Danish inflation,
# since this series is stationary according to the ADF and KPSS tests.

fit_arima_benchmark <- macro_data_diff %>%
  model(
    arima_benchmark = ARIMA(diff_dk_inflation)
  )

# Inspect selected ARIMA model
report(fit_arima_benchmark)

# Information criteria
glance(fit_arima_benchmark) %>%
  dplyr::select(.model, AIC, AICc, BIC) %>%
  kable(digits = 2, align = "c")


# ============================================================
# 19. Residual diagnostics for ARIMA benchmark
# ============================================================

# Plot residuals, ACF of residuals, and histogram
fit_arima_benchmark %>%
  gg_tsresiduals()


# Ljung-Box test for residual autocorrelation
# H0: residuals are independently distributed / no autocorrelation

augment(fit_arima_benchmark) %>%
  features(.innov, ljung_box, lag = 24, dof = 0) %>%
  kable(digits = 4, align = "c")



# ============================================================
# 20. Dynamic regression model: Phillips curve specification
# ============================================================

# The dynamic regression model tests whether Danish unemployment adds
# explanatory power for Danish inflation beyond inflation's own time-series
# dynamics.
#
# Dependent variable:
#   diff_dk_inflation
#
# Explanatory variable:
#   diff_dk_unemployment
#
# The ARIMA errors allow the model to capture remaining autocorrelation and
# possible seasonal dynamics in the residuals.

fit_dyn <- macro_data_diff %>%
  model(
    dyn_reg = ARIMA(diff_dk_inflation ~ diff_dk_unemployment)
  )

# Print model estimates
report(fit_dyn)





# ============================================================
# 21. Compare ARIMA benchmark and dynamic regression
# ============================================================

fit_compare <- macro_data_diff %>%
  model(
    arima_benchmark = ARIMA(diff_dk_inflation),
    dyn_reg = ARIMA(diff_dk_inflation ~ diff_dk_unemployment)
  )

# Inspect selected models
fit_compare

# Compare information criteria
glance(fit_compare) %>%
  dplyr::select(.model, AIC, AICc, BIC) %>%
  arrange(AICc) %>%
  kable(digits = 2, align = "c")




# ============================================================
# 22. Residual diagnostics for dynamic regression
# ============================================================

# Plot residuals, residual ACF, and histogram
fit_dyn %>%
  gg_tsresiduals()

# Ljung-Box test for residual autocorrelation
# H0: no residual autocorrelation

augment(fit_dyn) %>%
  features(.innov, ljung_box, lag = 24, dof = 0) %>%
  kable(digits = 4, align = "c")





