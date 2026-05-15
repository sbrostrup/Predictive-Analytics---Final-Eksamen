# ============================================================
# 0. Setup
# ============================================================

rm(list = ls())

if (!is.null(dev.list())) dev.off()

cat("\014")  # Clear console in RStudio

library(tidyverse)
library(lubridate)
library(knitr)
library(ggplot2)
library(patchwork)
library(fpp3)
library(tseries)
library(feasts)
library(zoo)
library(strucchange)
library(gets)
library(lmtest)
library(httr)
library(jsonlite)


# ============================================================
# 1. Load data from Eurostat API
# ============================================================

# Helper function: fetch and parse the Eurostat JSON-stat API.
# Returns a long-format tibble with columns: geo, dato, vaerdi.
hent_eurostat <- function(url) {
  svar <- httr::GET(url)
  
  if (httr::http_error(svar)) {
    stop("API-fejl: ", httr::status_code(svar))
  }
  
  data <- jsonlite::fromJSON(
    httr::content(svar, as = "text", encoding = "UTF-8")
  )
  
  tider     <- names(data$dimension$time$category$index)
  geo_koder <- names(data$dimension$geo$category$index)
  
  n_tider <- length(tider)
  n_geo   <- length(geo_koder)
  
  matrix_vals <- matrix(
    NA_real_,
    nrow = n_geo,
    ncol = n_tider,
    dimnames = list(geo_koder, tider)
  )
  
  for (navn in names(data$value)) {
    idx     <- as.integer(navn)
    geo_idx <- idx %/% n_tider
    tid_idx <- idx %%  n_tider
    
    matrix_vals[geo_idx + 1, tid_idx + 1] <- data$value[[navn]]
  }
  
  as.data.frame(matrix_vals) |>
    tibble::rownames_to_column("geo") |>
    tidyr::pivot_longer(
      cols = -geo,
      names_to = "periode",
      values_to = "vaerdi"
    ) |>
    dplyr::mutate(dato = lubridate::ym(periode)) |>
    dplyr::select(geo, dato, vaerdi) |>
    dplyr::arrange(geo, dato)
}


# --- 1.1 Inflation: HICP annual rate of change (YoY %) ---

inflation_gammel <- hent_eurostat(paste0(
  "https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data/",
  "prc_hicp_manr?geo=EU27_2020&geo=DK&coicop=CP00&unit=RCH_A&lang=en"
))

inflation_ny <- hent_eurostat(paste0(
  "https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data/",
  "prc_hicp_minr?geo=EU27_2020&geo=DK&coicop18=TOTAL&unit=RCH_A&lang=en",
  "&sinceTimePeriod=2026-01"
))

inflation <- dplyr::bind_rows(inflation_gammel, inflation_ny) |>
  dplyr::distinct(geo, dato, .keep_all = TRUE) |>
  dplyr::arrange(geo, dato)


# --- 1.2 Unemployment: Not Seasonally Adjusted monthly rate ---

arbejdsloes <- hent_eurostat(paste0(
  "https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data/",
  "une_rt_m?geo=EU27_2020&geo=DK&sex=T&age=TOTAL&s_adj=NSA&unit=PC_ACT&lang=en"
))


# --- 1.3 Pivot to wide format and join ---

inflation_wide <- inflation |>
  tidyr::pivot_wider(names_from = geo, values_from = vaerdi) |>
  dplyr::rename(
    dk_inflation = DK,
    eu_inflation = EU27_2020
  )

arbejdsloes_wide <- arbejdsloes |>
  tidyr::pivot_wider(names_from = geo, values_from = vaerdi) |>
  dplyr::rename(
    dk_unemployment = DK,
    eu_unemployment = EU27_2020
  )

df <- dplyr::full_join(inflation_wide, arbejdsloes_wide, by = "dato") |>
  dplyr::arrange(dato) |>
  tidyr::drop_na()

cat("Rows:", nrow(df), "\n")
cat("Period:", format(min(df$dato), "%b %Y"), "to",
    format(max(df$dato), "%b %Y"), "\n")


# ============================================================
# 2. Data Preparation
# ============================================================

combined <- df |>
  dplyr::mutate(yearmonth = tsibble::yearmonth(dato)) |>
  dplyr::select(
    yearmonth,
    dk_inflation,
    eu_inflation,
    dk_unemployment,
    eu_unemployment
  ) |>
  tsibble::as_tsibble(index = yearmonth) |>
  tidyr::drop_na()

combined_dk <- combined |>
  dplyr::select(yearmonth, dk_inflation, dk_unemployment)

head(combined) |> knitr::kable(align = "c", digits = 2)
tail(combined) |> knitr::kable(align = "c", digits = 2)

cat("Sample period:", as.character(min(combined$yearmonth)),
    "to", as.character(max(combined$yearmonth)), "\n")
cat("Number of monthly observations:", nrow(combined), "\n")


# ============================================================
# 3. Data Visualisation
# ============================================================

p_dk_infl <- combined |>
  feasts::autoplot(dk_inflation) +
  geom_smooth(method = "lm", se = FALSE, colour = "red", linetype = "dashed") +
  labs(
    title = "Danish Inflation (YoY, %)",
    x = "Month",
    y = "Inflation (%)"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

p_dk_unemp <- combined |>
  feasts::autoplot(dk_unemployment) +
  geom_smooth(method = "lm", se = FALSE, colour = "red", linetype = "dashed") +
  labs(
    title = "Danish Unemployment Rate (%)",
    x = "Month",
    y = "Unemployment (%)"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

p_eu_infl <- combined |>
  feasts::autoplot(eu_inflation) +
  geom_smooth(method = "lm", se = FALSE, colour = "red", linetype = "dashed") +
  labs(
    title = "EU27 Inflation (YoY, %)",
    x = "Month",
    y = "Inflation (%)"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

p_eu_unemp <- combined |>
  feasts::autoplot(eu_unemployment) +
  geom_smooth(method = "lm", se = FALSE, colour = "red", linetype = "dashed") +
  labs(
    title = "EU27 Unemployment Rate (%)",
    x = "Month",
    y = "Unemployment (%)"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

(p_dk_infl | p_dk_unemp) / (p_eu_infl | p_eu_unemp)


# --- Seasonal plots ---

combined |>
  feasts::gg_season(dk_inflation) +
  labs(
    title = "Seasonal Plot: Danish Inflation (YoY, %)",
    x = "Month",
    y = "Inflation (%)"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

combined |>
  feasts::gg_subseries(dk_inflation) +
  labs(
    title = "Seasonal Subseries Plot: Danish Inflation (YoY, %)",
    x = "Month",
    y = "Inflation (%)"
  ) +
  theme(plot.title = element_text(hjust = 0.5))


# --- STL decomposition ---

# The STL decomposition indicates only weak seasonality in Danish YoY inflation.
# This is expected because the inflation variable is already measured as a
# year-on-year percentage change, which removes much of the regular monthly
# seasonal pattern.

combined |>
  fabletools::model(
    STL(
      dk_inflation ~ trend(window = 25) +
        season(window = "periodic"),
      robust = TRUE
    )
  ) |>
  fabletools::components() |>
  autoplot() +
  labs(title = "STL Decomposition: Danish Inflation")


# --- ACF / PACF displays ---

# Since the ADF and KPSS tests below indicate that Danish inflation is
# stationary in levels, we identify AR and MA terms using the ACF/PACF of
# dk_inflation in levels, not the first-differenced series.

combined |>
  feasts::gg_tsdisplay(dk_inflation, plot_type = "partial") +
  labs(title = "Danish Inflation: Time Series, ACF and PACF")

combined |>
  feasts::gg_tsdisplay(dk_unemployment, plot_type = "partial") +
  labs(title = "Danish Unemployment: Time Series, ACF and PACF")

combined |>
  feasts::gg_tsdisplay(eu_inflation, plot_type = "partial") +
  labs(title = "EU27 Inflation: Time Series, ACF and PACF")

combined |>
  feasts::gg_tsdisplay(eu_unemployment, plot_type = "partial") +
  labs(title = "EU27 Unemployment: Time Series, ACF and PACF")


# ============================================================
# 4. Stationarity Tests
# ============================================================

# Important:
# ADF H0: the series has a unit root / is non-stationary.
# KPSS H0: the series is stationary.
#
# Our results indicate:
#   - Danish inflation: stationary in levels.
#   - Danish unemployment: non-stationary in levels.
#   - EU inflation: mixed/borderline evidence.
#   - EU unemployment: non-stationary in levels.
#
# Therefore, Danish inflation does not need first differencing for the
# baseline ARIMA model. We model it in levels with d = 0.
# The original inflation series does not need to be white noise. Rather,
# the residuals from the fitted ARIMA/ARMA model should behave like white noise.

cat("\n=== ADF Tests on levels (H0: Non-stationary) ===\n")
tseries::adf.test(na.omit(combined$dk_inflation))
tseries::adf.test(na.omit(combined$dk_unemployment))
tseries::adf.test(na.omit(combined$eu_inflation))
tseries::adf.test(na.omit(combined$eu_unemployment))

cat("\n=== KPSS Tests on levels (H0: Stationary) ===\n")
tseries::kpss.test(na.omit(combined$dk_inflation))
tseries::kpss.test(na.omit(combined$dk_unemployment))
tseries::kpss.test(na.omit(combined$eu_inflation))
tseries::kpss.test(na.omit(combined$eu_unemployment))


# --- First differences, kept for optional robustness checks ---

# We still create first differences because they may be useful later,
# especially for unemployment and for robustness checks.
# However, the baseline ARIMA model below uses dk_inflation in levels.

combined_dk <- combined_dk |>
  dplyr::mutate(
    diff_dk_inflation    = difference(dk_inflation),
    diff_dk_unemployment = difference(dk_unemployment)
  )

combined <- combined |>
  dplyr::mutate(
    diff_dk_inflation    = difference(dk_inflation),
    diff_eu_inflation    = difference(eu_inflation),
    diff_dk_unemployment = difference(dk_unemployment),
    diff_eu_unemployment = difference(eu_unemployment)
  )


# --- Optional: quickly inspect first-differenced inflation again ---

# Uncomment this block if you want to inspect differenced inflation.
# This is not used for the baseline ARIMA model because dk_inflation is
# stationary in levels according to ADF and KPSS.

# combined |>
#   feasts::gg_tsdisplay(diff_dk_inflation, plot_type = "partial") +
#   labs(title = "First-Differenced Danish Inflation: ACF and PACF")

# cat("\n=== Optional ADF/KPSS on first-differenced Danish inflation ===\n")
# tseries::adf.test(na.omit(combined$diff_dk_inflation))
# tseries::kpss.test(na.omit(combined$diff_dk_inflation))


# --- Optional: inspect all first differences ---

# combined |>
#   feasts::gg_tsdisplay(diff_eu_inflation, plot_type = "partial") +
#   labs(title = "First-Differenced EU27 Inflation: ACF and PACF")

# combined |>
#   feasts::gg_tsdisplay(diff_dk_unemployment, plot_type = "partial") +
#   labs(title = "First-Differenced Danish Unemployment: ACF and PACF")

# combined |>
#   feasts::gg_tsdisplay(diff_eu_unemployment, plot_type = "partial") +
#   labs(title = "First-Differenced EU27 Unemployment: ACF and PACF")


# ============================================================
# 5. Structural Break Test: QLR
# ============================================================

# Even though Danish inflation appears stationary in levels, the series
# contains large temporary shocks, especially around 2021-2023.
# A structural break test is useful as a robustness check.

dk_infl_ts <- ts(
  combined$dk_inflation,
  start = c(
    lubridate::year(min(combined$yearmonth)),
    lubridate::month(min(combined$yearmonth))
  ),
  frequency = 12
)

qlr <- strucchange::Fstats(dk_infl_ts ~ 1, from = 0.1, to = 0.9)
strucchange::sctest(qlr, type = "supF")

Box.test(residuals(lm(dk_infl_ts ~ 1)), lag = 12, type = "Ljung-Box")


# Find det estimerede brudpunkt
breakpoints(dk_infl_ts ~ 1)

# Eller se hvilken dato sup-F optræder
qlr$breakpoint

#Dansk inflation har to strukturelle brud:
#  → December 2012  (overgang til lavinflationsperiode)
#  → August 2021    (inflationschok)


bp <- strucchange::breakpoints(dk_infl_ts ~ 1)
summary(bp)

plot(
  qlr,
  main = "QLR Test: Structural Breaks in Danish Inflation",
  xlab = "Time",
  ylab = "F-statistic"
)

lines(bp, col = "red")

bp_date <- combined$yearmonth[bp$breakpoints]
cat("Estimated structural break date(s):", as.character(bp_date), "\n")


# ============================================================
# 6. ARIMA Baseline Model
# ============================================================

# Since Danish inflation is stationary in levels, the baseline model is
# estimated as an ARMA/SARMA model, i.e. ARIMA with d = 0.
#
# The ACF/PACF suggests strong persistence:
#   - ACF decays gradually
#   - PACF has a large spike at lag 1
#
# This points toward an AR-type model, such as AR(1), AR(2), or ARMA(1,1).
# We compare several candidate models and then check whether the residuals
# behave like white noise.

test_end   <- max(combined_dk$yearmonth)
test_start <- test_end - 23

train <- combined_dk |>
  dplyr::filter(yearmonth < test_start)

test <- combined_dk |>
  dplyr::filter(yearmonth >= test_start)

cat("Training period:", as.character(min(train$yearmonth)),
    "to", as.character(max(train$yearmonth)),
    "(", nrow(train), "obs )\n")

cat("Test period    :", as.character(min(test$yearmonth)),
    "to", as.character(max(test$yearmonth)),
    "(", nrow(test), "obs )\n")


# --- 6.1 Automatic model selection with no differencing ---

# This model allows fable to choose AR and MA terms, but fixes d = 0.
# We also allow simple seasonal ARMA terms, but no seasonal differencing.

fit_auto_d0 <- train |>
  fabletools::model(
    auto_d0 = fable::ARIMA(
      dk_inflation ~ pdq(p = 0:6, d = 0, q = 0:6) +
        PDQ(P = 0:2, D = 0, Q = 0:2),
      stepwise = FALSE,
      approximation = FALSE
    )
  )

report(fit_auto_d0)


# --- 6.2 Manual candidate ARMA/SARMA models ---

# These candidate models are chosen based on:
#   - stationarity tests: d = 0
#   - ACF/PACF: strong persistence and likely AR structure
#   - possible weak annual dynamics around lag 12

fit_candidates <- train |>
  fabletools::model(
    ar1      = fable::ARIMA(dk_inflation ~ pdq(1,0,0)),
    ar2      = fable::ARIMA(dk_inflation ~ pdq(2,0,0)),
    ar3      = fable::ARIMA(dk_inflation ~ pdq(3,0,0)),
    ma1      = fable::ARIMA(dk_inflation ~ pdq(0,0,1)),
    ma2      = fable::ARIMA(dk_inflation ~ pdq(0,0,2)),
    arma11   = fable::ARIMA(dk_inflation ~ pdq(1,0,1)),
    arma12   = fable::ARIMA(dk_inflation ~ pdq(1,0,2)),
    arma21   = fable::ARIMA(dk_inflation ~ pdq(2,0,1)),
    arma22   = fable::ARIMA(dk_inflation ~ pdq(2,0,2)),
    sar_ar1  = fable::ARIMA(dk_inflation ~ pdq(1,0,0) + PDQ(1,0,0)),
    sar_ma1  = fable::ARIMA(dk_inflation ~ pdq(1,0,0) + PDQ(0,0,1)),
    sarma11  = fable::ARIMA(dk_inflation ~ pdq(1,0,1) + PDQ(1,0,0)),
    sarma12  = fable::ARIMA(dk_inflation ~ pdq(1,0,1) + PDQ(0,0,1))
  )

fit_candidates |>
  fabletools::glance() |>
  dplyr::select(.model, AIC, AICc, BIC, sigma2) |>
  dplyr::arrange(AICc) |>
  knitr::kable(digits = 3, align = "c")


# --- 6.3 Compare automatic d = 0 model with manual candidates ---

dplyr::bind_rows(
  fit_auto_d0 |>
    fabletools::glance() |>
    dplyr::mutate(model_group = "Automatic d = 0"),
  
  fit_candidates |>
    fabletools::glance() |>
    dplyr::mutate(model_group = "Manual candidates")
) |>
  dplyr::select(model_group, .model, AIC, AICc, BIC, sigma2) |>
  dplyr::arrange(AICc) |>
  knitr::kable(digits = 3, align = "c")


# --- 6.4 Residual diagnostics for automatic d = 0 model ---

# A good ARIMA/ARMA model does not require the original series to be white noise.
# Instead, the residuals should behave approximately like white noise.

fit_auto_d0 |>
  feasts::gg_tsresiduals() +
  labs(title = "ARMA/SARMA Baseline: Residual Diagnostics")

dof_auto_d0 <- fit_auto_d0 |>
  broom::tidy() |>
  dplyr::filter(stringr::str_detect(term, "^ar|^ma|^sar|^sma")) |>
  nrow()

fit_auto_d0 |>
  fabletools::augment() |>
  feasts::features(.innov, ljung_box, dof = dof_auto_d0, lag = 24) |>
  knitr::kable(digits = 3, align = "c")


# --- 6.5 Forecast from automatic d = 0 model ---

fc_auto_d0 <- fit_auto_d0 |>
  fabletools::forecast(h = nrow(test))

fc_auto_d0 |>
  autoplot(
    combined_dk |> dplyr::filter(yearmonth >= yearmonth("2020 Jan")),
    level = c(80, 95)
  ) +
  labs(
    title = "ARMA/SARMA Baseline: Forecast vs. Actual Danish Inflation",
    x = "Month",
    y = "Inflation (YoY, %)"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

fc_auto_d0 |>
  fabletools::accuracy(test) |>
  dplyr::select(.model, RMSE, MAE, MAPE, MASE) |>
  knitr::kable(digits = 3, align = "c")


# ============================================================
# Stop here for now
# ============================================================

# The next step would be to select the preferred ARMA/SARMA specification
# based on AICc/BIC, residual diagnostics, and forecast accuracy.
# After that, the dynamic regression / Phillips curve model can be updated.