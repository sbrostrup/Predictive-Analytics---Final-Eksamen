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
library(strucchange) # Structural break tests (QLR)
library(gets)        # GETS model selection + indicator saturation
library(lmtest)      # Granger causality tests

# theme_set(theme_minimal())  # Commented out to keep fpp3 default plot style


# ============================================================
# 1. Load data from Eurostat CSV exports
# ============================================================

# Data files (Eurostat "linear" CSV exports) live alongside the script.
# Each file has columns: DATAFLOW, LAST UPDATE, freq, ..., geo, TIME_PERIOD,
# OBS_VALUE, OBS_FLAG, CONF_STATUS. The geo column holds the full English
# label rather than the short code, so we map it back to DK / EU27_2020 to
# keep the downstream pipeline unchanged.

script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile, mustWork = FALSE)),
  error = function(e) getwd()
)
if (!nzchar(script_dir) || is.na(script_dir)) script_dir <- getwd()

geo_label_to_code <- c(
  "Denmark"                                      = "DK",
  "European Union - 27 countries (from 2020)"    = "EU27_2020"
)

# Helper: read a Eurostat linear CSV and return a long tibble with the same
# (geo, dato, vaerdi) contract the rest of the script expects.
laes_eurostat_csv <- function(sti) {
  read_csv(sti, show_col_types = FALSE) |>
    filter(geo %in% names(geo_label_to_code)) |>
    transmute(
      geo    = unname(geo_label_to_code[geo]),
      dato   = ym(TIME_PERIOD),
      vaerdi = as.numeric(OBS_VALUE)
    ) |>
    filter(!is.na(dato)) |>
    arrange(geo, dato)
}

# --- 1.1 Inflation: HICP annual rate of change (YoY %) ---
# The minr CSV already covers the full sample (1997-01 onwards) for both DK
# and EU27, so no manr/minr stitching is needed.

inflation <- laes_eurostat_csv(
  file.path(script_dir, "prc_hicp_minr__custom_21465689_linear.csv")
)

# --- 1.2 Unemployment: Not Seasonally Adjusted (NSA) monthly rate ---
# NSA data preserves the original unadjusted series, consistent with the
# raw HICP data. The regular seasonal pattern in NSA unemployment is visible
# in the ACF after first differencing but does not affect stationarity.
# When used as a regressor in the dynamic regression model, the ARIMA error
# structure absorbs any residual seasonal dynamics.

arbejdsloes <- laes_eurostat_csv(
  file.path(script_dir, "une_rt_m__custom_21465700_linear.csv")
)

# --- 1.3 Pivot to wide format and join ---

inflation_wide <- inflation |>
  pivot_wider(names_from = geo, values_from = vaerdi) |>
  rename(dk_inflation = DK, eu_inflation = EU27_2020)

arbejdsloes_wide <- arbejdsloes |>
  pivot_wider(names_from = geo, values_from = vaerdi) |>
  rename(dk_unemployment = DK, eu_unemployment = EU27_2020)

df <- full_join(inflation_wide, arbejdsloes_wide, by = "dato") |>
  arrange(dato) |>
  drop_na()

cat("Rows:", nrow(df), "\n")
cat("Period:", format(min(df$dato), "%b %Y"), "to",
    format(max(df$dato), "%b %Y"), "\n")


# ============================================================
# 2. Data Preparation
# ============================================================

# Key improvement over the FRED version:
#   Inflation is delivered directly as HICP YoY % change by Eurostat —
#   no manual lag-12 computation needed. Both DK and EU27 series are
#   therefore directly comparable and on the same scale.
#
# Unemployment is sourced as Not Seasonally Adjusted (NSA), preserving
# the unadjusted series. The regular seasonal calendar pattern remains
# visible in the ACF of first-differenced NSA unemployment (spike at lag 12),
# but this does not affect stationarity — ADF and KPSS confirm I(0) after
# first differencing. When used as a regressor in the dynamic regression,
# the ARIMA error structure absorbs any residual seasonal dynamics.
#
# Transformation choice (log / Box-Cox):
#   Both inflation and unemployment are already expressed as percentage
#   values — scale-free relative measures. Further log-transformation
#   is non-standard for these series and would complicate interpretation.
#   We work directly with the percentage-point values throughout.

combined <- df %>%
  mutate(yearmonth = yearmonth(dato)) %>%
  select(yearmonth, dk_inflation, eu_inflation,
         dk_unemployment, eu_unemployment) %>%
  as_tsibble(index = yearmonth) %>%
  drop_na()

# Danish-only tsibble: used for ARIMA, dynamic regression, and TSCV.
# Since all four Eurostat series share the same sample window, combined_dk
# is simply a column-subset of combined with no additional truncation.
combined_dk <- combined %>%
  select(yearmonth, dk_inflation, dk_unemployment)

head(combined) %>% kable(align = "c", digits = 2)
tail(combined) %>% kable(align = "c", digits = 2)

cat("Sample period:", as.character(min(combined$yearmonth)),
    "to", as.character(max(combined$yearmonth)), "\n")
cat("Number of monthly observations:", nrow(combined), "\n")


# ============================================================
# 3. Data Visualisation
# ============================================================

# --- 3.1 Time series plots (levels) ---

# A linear trend line (red, dashed) makes any upward or downward drift visible.

p_dk_infl <- combined %>%
  autoplot(dk_inflation) +
  geom_smooth(method = "lm", se = FALSE, colour = "red", linetype = "dashed") +
  labs(title = "Danish Inflation (YoY, %)",
       x = "Month", y = "Inflation (%)") +
  theme(plot.title = element_text(hjust = 0.5))

p_dk_unemp <- combined %>%
  autoplot(dk_unemployment) +
  geom_smooth(method = "lm", se = FALSE, colour = "red", linetype = "dashed") +
  labs(title = "Danish Unemployment Rate (%)",
       x = "Month", y = "Unemployment (%)") +
  theme(plot.title = element_text(hjust = 0.5))

p_eu_infl <- combined %>%
  autoplot(eu_inflation) +
  geom_smooth(method = "lm", se = FALSE, colour = "red", linetype = "dashed") +
  labs(title = "EU27 Inflation (YoY, %)",
       x = "Month", y = "Inflation (%)") +
  theme(plot.title = element_text(hjust = 0.5))

p_eu_unemp <- combined %>%
  autoplot(eu_unemployment) +
  geom_smooth(method = "lm", se = FALSE, colour = "red", linetype = "dashed") +
  labs(title = "EU27 Unemployment Rate (%)",
       x = "Month", y = "Unemployment (%)") +
  theme(plot.title = element_text(hjust = 0.5))

(p_dk_infl | p_dk_unemp) / (p_eu_infl | p_eu_unemp)


# --- 3.2 Seasonal plots ---

# gg_season() overlays each calendar year as a separate line.
# Persistent crossing of lines in specific months signals seasonality.

combined %>%
  gg_season(dk_inflation) +
  labs(title = "Seasonal Plot: Danish Inflation (YoY, %)",
       x = "Month", y = "Inflation (%)") +
  theme(plot.title = element_text(hjust = 0.5))

combined %>%
  gg_season(dk_unemployment) +
  labs(title = "Seasonal Plot: Danish Unemployment (%)",
       x = "Month", y = "Unemployment (%)") +
  theme(plot.title = element_text(hjust = 0.5))

combined %>%
  gg_season(eu_inflation) +
  labs(title = "Seasonal Plot: EU27 Inflation (YoY, %)",
       x = "Month", y = "Inflation (%)") +
  theme(plot.title = element_text(hjust = 0.5))

combined %>%
  gg_season(eu_unemployment) +
  labs(title = "Seasonal Plot: EU27 Unemployment (%)",
       x = "Month", y = "Unemployment (%)") +
  theme(plot.title = element_text(hjust = 0.5))


# --- 3.3 Seasonal subseries plots ---

combined %>%
  gg_subseries(dk_inflation) +
  labs(title = "Seasonal Subseries Plot: Danish Inflation (YoY, %)",
       x = "Month", y = "Inflation (%)") +
  theme(plot.title = element_text(hjust = 0.5))

combined %>%
  gg_subseries(dk_unemployment) +
  labs(title = "Seasonal Subseries Plot: Danish Unemployment (%)",
       x = "Month", y = "Unemployment (%)") +
  theme(plot.title = element_text(hjust = 0.5))

combined %>%
  gg_subseries(eu_inflation) +
  labs(title = "Seasonal Subseries Plot: EU27 Inflation (YoY, %)",
       x = "Month", y = "Inflation (%)") +
  theme(plot.title = element_text(hjust = 0.5))

combined %>%
  gg_subseries(eu_unemployment) +
  labs(title = "Seasonal Subseries Plot: EU27 Unemployment (%)",
       x = "Month", y = "Unemployment (%)") +
  theme(plot.title = element_text(hjust = 0.5))


# --- 3.4 ACF / PACF displays ---

combined %>%
  gg_tsdisplay(dk_inflation, plot_type = "partial") +
  labs(title = "Danish Inflation: Time Series, ACF and PACF")

combined %>%
  gg_tsdisplay(dk_unemployment, plot_type = "partial") +
  labs(title = "Danish Unemployment: Time Series, ACF and PACF")

combined %>%
  gg_tsdisplay(eu_inflation, plot_type = "partial") +
  labs(title = "EU27 Inflation: Time Series, ACF and PACF")

combined %>%
  gg_tsdisplay(eu_unemployment, plot_type = "partial") +
  labs(title = "EU27 Unemployment: Time Series, ACF and PACF")


# ============================================================
# 4. Stationarize the Data
# ============================================================

# We test stationarity with two complementary tests:
#
#   ADF (Augmented Dickey-Fuller)  H0: Series has a unit root (non-stationary)
#   KPSS                           H0: Series is stationary
#
# A series is confirmed non-stationary when:
#   - ADF fails to reject H0 (p > 0.05), AND
#   - KPSS rejects H0 (p < 0.05)


# --- 4.1 Stationarity tests on levels ---

cat("\n=== ADF Tests on levels (H0: Non-stationary) ===\n")
adf.test(na.omit(combined$dk_inflation))
adf.test(na.omit(combined$dk_unemployment))
adf.test(na.omit(combined$eu_inflation))
adf.test(na.omit(combined$eu_unemployment))

cat("\n=== KPSS Tests on levels (H0: Stationary) ===\n")
kpss.test(na.omit(combined$dk_inflation))
kpss.test(na.omit(combined$dk_unemployment))
kpss.test(na.omit(combined$eu_inflation))
kpss.test(na.omit(combined$eu_unemployment))


# --- 4.2 Transformations: first differences ---

# All four series are I(1): non-stationary in levels, stationary after one
# first difference. We therefore apply a standard first difference to each.
#
# The unemployment series are sourced as Not Seasonally Adjusted (NSA), so
# a residual seasonal pattern remains visible in the ACF of the differenced
# series (spike at lag 12). This is not a stationarity problem — ADF and KPSS
# confirm both differenced unemployment series are I(0) — but rather a nuisance
# that is absorbed by the seasonal component of the ARIMA error structure when
# the series are used as regressors in the dynamic regression model.

combined_dk <- combined_dk %>%
  mutate(
    diff_dk_inflation    = difference(dk_inflation),
    diff_dk_unemployment = difference(dk_unemployment)
  )

combined <- combined %>%
  mutate(
    diff_dk_inflation    = difference(dk_inflation),
    diff_eu_inflation    = difference(eu_inflation),
    diff_dk_unemployment = difference(dk_unemployment),
    diff_eu_unemployment = difference(eu_unemployment)
  )


# --- 4.3 Re-test after first differencing ---

cat("\n=== ADF Tests on first differences (H0: Non-stationary) ===\n")
adf.test(na.omit(combined$diff_dk_inflation))
adf.test(na.omit(combined$diff_eu_inflation))
adf.test(na.omit(combined$diff_dk_unemployment))
adf.test(na.omit(combined$diff_eu_unemployment))

cat("\n=== KPSS Tests on first differences (H0: Stationary) ===\n")
kpss.test(na.omit(combined$diff_dk_inflation))
kpss.test(na.omit(combined$diff_eu_inflation))
kpss.test(na.omit(combined$diff_dk_unemployment))
kpss.test(na.omit(combined$diff_eu_unemployment))


# --- 4.4 Visual confirmation after transformation ---

combined %>%
  gg_tsdisplay(diff_dk_inflation, plot_type = "partial") +
  labs(title = "First-Differenced Danish Inflation: ACF and PACF")

combined %>%
  gg_tsdisplay(diff_eu_inflation, plot_type = "partial") +
  labs(title = "First-Differenced EU27 Inflation: ACF and PACF")

combined %>%
  gg_tsdisplay(diff_dk_unemployment, plot_type = "partial") +
  labs(title = "First-Differenced Danish Unemployment: ACF and PACF")

combined %>%
  gg_tsdisplay(diff_eu_unemployment, plot_type = "partial") +
  labs(title = "First-Differenced EU27 Unemployment: ACF and PACF")


# --- 4.5 Structural break test: Quandt Likelihood Ratio (QLR) ---

# The QLR test scans all candidate break dates and returns the one with the
# highest F-statistic. H0: No structural break.
# We trim 15% from each end as is conventional.

dk_infl_ts <- ts(
  combined$dk_inflation,
  start     = c(year(min(combined$yearmonth)),
                month(min(combined$yearmonth))),
  frequency = 12
)

qlr <- Fstats(dk_infl_ts ~ 1, from = 0.15, to = 0.85)
sctest(qlr, type = "supF")

bp <- breakpoints(dk_infl_ts ~ 1)
summary(bp)

plot(qlr,
     main = "QLR Test: Structural Breaks in Danish Inflation",
     xlab = "Time",
     ylab = "F-statistic")
lines(bp, col = "red")

bp_date <- combined$yearmonth[bp$breakpoints]
cat("Estimated structural break date(s):", as.character(bp_date), "\n")


# ============================================================
# 5. ARIMA Baseline Model
# ============================================================

# The ARIMA model serves as the baseline: it forecasts Danish inflation
# using only its own past values. A stronger model (dynamic regression)
# should outperform this benchmark.
#
# Strategy:
#   1. Split the sample into a training set and a test set.
#   2. Auto-select the best ARIMA on the training set.
#   3. Estimate a grid of manual candidates and compare AICc / BIC.
#   4. Check residual diagnostics.
#   5. Forecast over the test period and evaluate accuracy.


# --- 5.1 Train / test split ---

# Test set: the most recent 24 months of available Danish data.
# Training set: everything before the test period.
# A 24-month (2-year) hold-out is consistent with the evaluation window
# used throughout the course problem sets.

test_end   <- max(combined_dk$yearmonth)
test_start <- test_end - 23   # 24-month window

train <- combined_dk %>% filter(yearmonth <  test_start)
test  <- combined_dk %>% filter(yearmonth >= test_start)

cat("Training period:", as.character(min(train$yearmonth)),
    "to", as.character(max(train$yearmonth)),
    "(", nrow(train), "obs )\n")
cat("Test period    :", as.character(min(test$yearmonth)),
    "to", as.character(max(test$yearmonth)),
    "(", nrow(test),  "obs )\n")


# --- 5.2 Automatic model selection ---

fit_auto <- train %>%
  model(auto = ARIMA(dk_inflation, stepwise = FALSE, approximation = FALSE))

report(fit_auto)


# --- 5.3 Manual candidate models ---

# Informed by the ACF / PACF of diff_dk_inflation:
#   - ACF spikes → MA terms likely
#   - PACF spikes → AR terms likely
#   - Seasonal spike at lag 12 → possible seasonal MA term
#
# We estimate ARIMA(p,1,q) and ARIMA(p,1,q)(P,0,Q)[12] candidates.

fit_candidates <- train %>%
  model(
    arima010   = ARIMA(dk_inflation ~ pdq(0,1,0)),
    arima011   = ARIMA(dk_inflation ~ pdq(0,1,1)),
    arima012   = ARIMA(dk_inflation ~ pdq(0,1,2)),
    arima110   = ARIMA(dk_inflation ~ pdq(1,1,0)),
    arima111   = ARIMA(dk_inflation ~ pdq(1,1,1)),
    arima112   = ARIMA(dk_inflation ~ pdq(1,1,2)),
    arima210   = ARIMA(dk_inflation ~ pdq(2,1,0)),
    arima211   = ARIMA(dk_inflation ~ pdq(2,1,1)),
    arima212   = ARIMA(dk_inflation ~ pdq(2,1,2)),
    arima011_s = ARIMA(dk_inflation ~ pdq(0,1,1) + PDQ(0,0,1)),
    arima111_s = ARIMA(dk_inflation ~ pdq(1,1,1) + PDQ(0,0,1)),
    arima211_s = ARIMA(dk_inflation ~ pdq(2,1,1) + PDQ(0,0,1))
  )

fit_candidates %>%
  glance() %>%
  select(.model, AIC, AICc, BIC, sigma2) %>%
  arrange(AICc) %>%
  kable(digits = 3, align = "c")


# --- 5.4 Select best model and check residuals ---

# The auto model is used as the final baseline (best AICc).
# The manual grid serves as a robustness check.

report(fit_auto)

fit_auto %>%
  gg_tsresiduals() +
  labs(title = "ARIMA Baseline: Residual Diagnostics")

# Ljung-Box: dof = number of estimated AR + MA parameters.
# Update dof based on the model reported above.
dof_arima <- fit_auto %>%
  tidy() %>%
  filter(str_detect(term, "^ar|^ma|^sar|^sma")) %>%
  nrow()

fit_auto %>%
  augment() %>%
  features(.innov, ljung_box, dof = dof_arima, lag = 24) %>%
  kable(digits = 3, align = "c")


# --- 5.5 Forecast and accuracy evaluation ---

fc_arima <- fit_auto %>%
  forecast(h = nrow(test))

fc_arima %>%
  autoplot(
    combined_dk %>% filter(yearmonth >= yearmonth("2020 Jan")),
    level = c(80, 95)
  ) +
  labs(
    title = "ARIMA Baseline: Forecast vs. Actual Danish Inflation",
    x     = "Month",
    y     = "Inflation (YoY, %)"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

fc_arima %>%
  accuracy(test) %>%
  select(.model, RMSE, MAE, MAPE, MASE) %>%
  kable(digits = 3, align = "c")


# ============================================================
# 6. Model 2: GETS Dynamic Regression (Phillips Curve)
# ============================================================

# The Phillips curve posits a negative relationship between unemployment and
# inflation: when unemployment falls the labour market tightens, wage growth
# accelerates, and eventually consumer prices rise. We model this as a dynamic
# regression of Danish inflation on Danish unemployment.
#
# Two-step strategy:
#   Step 1 — GETS model selection: start from a General Unrestricted Model
#            (GUM) with many lags and use indicator saturation (IIS, SIS) to
#            detect outliers and structural breaks while selecting the most
#            informative lag structure.
#   Step 2 — fpp3 dynamic regression: estimate the selected specification as
#            an ARIMA model with external regressors (ARIMAX) within the fpp3
#            framework, producing forecasts directly comparable to the ARIMA
#            baseline.


# --- 6.1 GETS: Data preparation ---

gets_train <- train %>%
  as_tibble() %>%
  arrange(yearmonth) %>%
  filter(!is.na(diff_dk_inflation), !is.na(diff_dk_unemployment))

y_tr <- gets_train$diff_dk_inflation
x_tr <- gets_train$diff_dk_unemployment
n_tr <- length(y_tr)


# --- 6.2 GETS: Build GUM regressor matrix (lags 0–6 of Δunemployment) ---

max_lag_gets <- 6
x_mat <- sapply(0:max_lag_gets, function(k) {
  c(rep(NA_real_, k), x_tr[seq_len(n_tr - k)])
})
colnames(x_mat) <- c("d_unemp_l0", paste0("d_unemp_l", seq_len(max_lag_gets)))

keep  <- complete.cases(x_mat)
y_gum <- y_tr[keep]
x_gum <- x_mat[keep, , drop = FALSE]


# --- 6.3 GETS: GUM with indicator saturation and model selection ---

# isat() fits the GUM and applies GETS reduction with indicator saturation
# in one call.
#
#   ar = 1:6    — autoregressive lags of Δinflation
#   mxreg       — contemporaneous and lags 0–6 of Δunemployment
#   iis = TRUE  — Impulse Indicator Saturation: detects individual outliers
#   sis = TRUE  — Step Indicator Saturation: detects permanent level shifts
#   t.pval      — retention threshold (|t| > 1.96)
#
# Retained sis## indicators should correspond to the QLR break dates.

gets_fit <- isat(
  y      = y_gum,
  ar     = 1:6,
  mxreg  = x_gum,
  iis    = TRUE,
  sis    = TRUE,
  t.pval = 0.05
)

print(gets_fit)


# --- 6.4 fpp3: Add lagged unemployment columns to the tsibble ---

combined_dk <- combined_dk %>%
  mutate(
    d_unemp_l0 = diff_dk_unemployment,
    d_unemp_l1 = lag(diff_dk_unemployment, 1),
    d_unemp_l2 = lag(diff_dk_unemployment, 2),
    d_unemp_l3 = lag(diff_dk_unemployment, 3)
  )

# Rebuild train / test from the updated combined_dk.
train <- combined_dk %>% filter(yearmonth <  test_start)
test  <- combined_dk %>% filter(yearmonth >= test_start)


# --- 6.5 fpp3: Fit dynamic regression with ARIMA errors ---

# ARIMA(dk_inflation ~ xreg) selects d = 1 automatically (dk_inflation is I(1)).
# The model is equivalent to regressing Δdk_inflation on the unemployment
# regressors with ARIMA(p,0,q) errors — both sides are I(0).

fit_dynreg <- train %>%
  model(
    dynreg = ARIMA(
      dk_inflation ~ d_unemp_l0 + d_unemp_l1 + d_unemp_l2 + d_unemp_l3,
      stepwise      = FALSE,
      approximation = FALSE
    )
  )

report(fit_dynreg)


# --- 6.6 Residual diagnostics ---

fit_dynreg %>%
  gg_tsresiduals() +
  labs(title = "Dynamic Regression: Residual Diagnostics")

# dof = estimated ARIMA parameters (AR + MA only; xreg excluded).
dof_dynreg <- fit_dynreg %>%
  tidy() %>%
  filter(str_detect(term, "^ar|^ma|^sar|^sma")) %>%
  nrow()

fit_dynreg %>%
  augment() %>%
  features(.innov, ljung_box, dof = dof_dynreg, lag = 24) %>%
  kable(digits = 3, align = "c")


# --- 6.7 Forecast and accuracy comparison ---

# Conditional forecast: conditions on actual test-period unemployment,
# isolating the Phillips curve mechanism from any unemployment forecast error.

fc_dynreg <- fit_dynreg %>%
  forecast(new_data = test)

fc_dynreg %>%
  autoplot(
    combined_dk %>% filter(yearmonth >= yearmonth("2020 Jan")),
    level = c(80, 95)
  ) +
  labs(
    title = "Dynamic Regression: Forecast vs. Actual Danish Inflation",
    x     = "Month",
    y     = "Inflation (YoY, %)"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

# Head-to-head accuracy: ARIMA baseline vs. dynamic regression.
bind_rows(
  fc_arima  %>% accuracy(test),
  fc_dynreg %>% accuracy(test)
) %>%
  select(.model, RMSE, MAE, MAPE, MASE) %>%
  arrange(RMSE) %>%
  kable(digits = 3, align = "c")

view(combined)

# ============================================================
# 7. Granger Causality: EU27 as Analytical Element
# ============================================================

# Denmark maintains a fixed exchange rate peg to the Euro under ERM II.
# EU27 monetary policy therefore directly shapes Danish financial conditions,
# and EU27 inflation may lead Danish inflation through import prices, energy
# markets, and inflation expectations. We test this using Granger causality.
#
# H0: X does NOT Granger-cause Y.
# A small p-value → reject H0 → X contains predictive information for Y.
# All series are first-differenced (I(0)) with lag order = 6.


# --- 7.1 Visual motivation: co-movement of DK and EU27 inflation ---

combined %>%
  select(yearmonth, dk_inflation, eu_inflation) %>%
  pivot_longer(
    cols      = c(dk_inflation, eu_inflation),
    names_to  = "series",
    values_to = "inflation"
  ) %>%
  autoplot(inflation, aes(colour = series)) +
  scale_colour_manual(
    values = c(dk_inflation = "steelblue", eu_inflation = "firebrick"),
    labels = c(dk_inflation = "Denmark", eu_inflation = "EU27")
  ) +
  labs(
    title  = "Danish vs. EU27 Inflation (YoY, %)",
    x      = "Month",
    y      = "Inflation (%)",
    colour = NULL
  ) +
  theme(plot.title = element_text(hjust = 0.5))

# Cross-correlation function: positive lags indicate EU27 leads DK.
combined %>%
  CCF(diff_dk_inflation, diff_eu_inflation, lag_max = 24) %>%
  autoplot() +
  labs(
    title = "Cross-Correlation: Δ DK Inflation and Δ EU27 Inflation",
    x     = "Lag (months; positive = EU27 leads)",
    y     = "Correlation"
  ) +
  theme(plot.title = element_text(hjust = 0.5))


# --- 7.2 Granger causality tests ---

granger_data <- combined %>%
  as_tibble() %>%
  filter(
    !is.na(diff_dk_inflation),
    !is.na(diff_eu_inflation),
    !is.na(diff_dk_unemployment),
    !is.na(diff_eu_unemployment)
  )

lag_order <- 6

cat("\n=== Granger Causality Tests — lag order:", lag_order, "===\n")
cat("(H0: X does NOT Granger-cause Y)\n\n")

# Test 1: EU27 inflation → DK inflation (main ERM II spillover hypothesis)
cat("--- (1) EU27 inflation → DK inflation ---\n")
grangertest(
  diff_dk_inflation ~ diff_eu_inflation,
  order = lag_order,
  data  = granger_data
)

# Test 2: EU27 unemployment → DK inflation
cat("--- (2) EU27 unemployment → DK inflation ---\n")
grangertest(
  diff_dk_inflation ~ diff_eu_unemployment,
  order = lag_order,
  data  = granger_data
)

# Test 3: DK unemployment → DK inflation (domestic Phillips curve check)
cat("--- (3) DK unemployment → DK inflation ---\n")
grangertest(
  diff_dk_inflation ~ diff_dk_unemployment,
  order = lag_order,
  data  = granger_data
)

# Test 4: EU27 inflation → DK unemployment
cat("--- (4) EU27 inflation → DK unemployment ---\n")
grangertest(
  diff_dk_unemployment ~ diff_eu_inflation,
  order = lag_order,
  data  = granger_data
)


# ============================================================
# 8. Time Series Cross-Validation (TSCV)
# ============================================================

# TSCV evaluates forecast accuracy across many expanding training windows
# rather than a single train/test split. At each step, one observation is
# added to the training set and a one-step-ahead forecast is produced.
# Averaging across all steps gives a robust estimate of forecast performance
# (rolling forecast origin evaluation, as described in the course framework).
#
# We apply TSCV to the ARIMA baseline. The dynamic regression model
# conditions on observed future unemployment at each step; TSCV for that
# model would require an unemployment forecast at every rolling origin,
# which is outside the scope of this analysis.


# --- 8.1 Build expanding-window dataset ---

# .init = 60: minimum training window of 60 months (5 years).
# .step = 1:  expand by one observation per fold.

combined_cv <- combined_dk %>%
  stretch_tsibble(.init = 60, .step = 1)

cat("Number of CV folds:", max(combined_cv$.id), "\n")
cat("First fold ends   :", as.character(combined_cv %>%
      filter(.id == 1) %>% pull(yearmonth) %>% max()), "\n")
cat("Last fold ends    :", as.character(combined_cv %>%
      filter(.id == max(.id)) %>% pull(yearmonth) %>% max()), "\n")


# --- 8.2 Fit ARIMA on every fold and forecast one step ahead ---

# stepwise = TRUE and approximation = TRUE for computational feasibility.
fit_cv <- combined_cv %>%
  model(
    arima_cv = ARIMA(dk_inflation, stepwise = TRUE, approximation = TRUE)
  )

fc_cv <- fit_cv %>%
  forecast(h = 1)


# --- 8.3 Accuracy across all folds ---

cv_accuracy <- fc_cv %>%
  accuracy(combined_dk) %>%
  select(.model, RMSE, MAE, MAPE, MASE)

cv_accuracy %>%
  kable(digits = 3, align = "c")


# --- 8.4 Compare TSCV vs. fixed test set accuracy ---

# TSCV averages 1-step-ahead errors over the full history.
# The fixed test set evaluates multi-step forecasting over the specific
# 24-month hold-out window. Together they give a complete picture.

bind_rows(
  cv_accuracy %>%
    mutate(.model = "ARIMA — TSCV (avg. 1-step, full history)"),
  fc_arima %>%
    accuracy(test) %>%
    select(.model, RMSE, MAE, MAPE, MASE) %>%
    mutate(.model = "ARIMA — fixed test set (24-step)")
) %>%
  kable(digits = 3, align = "c")

# Note on interpretation:
#   TSCV RMSE = average 1-step forecast error across the full sample.
#   Test set RMSE = average multi-step error during the specific hold-out period.
#   The two numbers are not directly comparable (different horizons and periods)
#   but together give a complete picture of model performance.
