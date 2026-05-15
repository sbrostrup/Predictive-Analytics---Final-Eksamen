# ============================================================
# 0. Setup
# ============================================================

rm(list = ls())

library(tidyverse)
library(lubridate)
library(knitr)
library(fpp3)
library(tseries)
library(vars)
library(dplyr)
library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")


# ============================================================
# 1. Indlæs data
# ============================================================

find_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0)
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)))
  frame_files <- vapply(sys.frames(), function(x) {
    ofile <- x$ofile; if (is.null(ofile)) "" else ofile
  }, character(1))
  frame_files <- frame_files[nzchar(frame_files)]
  if (length(frame_files) > 0)
    return(dirname(normalizePath(frame_files[length(frame_files)], mustWork = FALSE)))
  getwd()
}

find_data_dir <- function(start_dir) {
  candidate_names <- c("Data - csv.", "Data - csv")
  current_dir <- normalizePath(start_dir, mustWork = FALSE)
  repeat {
    existing <- file.path(current_dir, candidate_names)[dir.exists(file.path(current_dir, candidate_names))]
    if (length(existing) == 1) return(existing)
    parent_dir <- dirname(current_dir)
    if (identical(parent_dir, current_dir)) break
    current_dir <- parent_dir
  }
  stop("Kunne ikke finde mappen 'Data - csv'.")
}

geo_label_to_code <- c(
  "Denmark"                                   = "DK",
  "European Union - 27 countries (from 2020)" = "EU27_2020"
)

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

data_dir <- find_data_dir(find_script_dir())

inflation_wide <- laes_eurostat_csv(
  file.path(data_dir, "prc_hicp_minr__custom_21465689_linear.csv")
) |>
  pivot_wider(names_from = geo, values_from = vaerdi) |>
  rename(dk_inflation = DK, eu_inflation = EU27_2020)

arbejdsloes_wide <- laes_eurostat_csv(
  file.path(data_dir, "une_rt_m__custom_21465700_linear.csv")
) |>
  pivot_wider(names_from = geo, values_from = vaerdi) |>
  rename(dk_unemployment = DK, eu_unemployment = EU27_2020)

df <- full_join(inflation_wide, arbejdsloes_wide, by = "dato") |>
  arrange(dato) |>
  drop_na()


# ============================================================
# 2. Dataforberedelse
# ============================================================

combined_dk <- df %>%
  mutate(yearmonth = yearmonth(dato)) %>%
  select(yearmonth, dk_inflation, dk_unemployment) %>%
  as_tsibble(index = yearmonth) %>%
  drop_na() %>%
  mutate(
    diff_dk_inflation    = difference(dk_inflation),
    diff_dk_unemployment = difference(dk_unemployment)
  )

cat("Sample:", as.character(min(combined_dk$yearmonth)),
    "to", as.character(max(combined_dk$yearmonth)),
    "(", nrow(combined_dk), "obs )\n")


# ============================================================
# 3. Stationaritetstests
# ============================================================

# ADF  H0: Ikke-stationær  → p < 0.05 = stationær
# KPSS H0: Stationær       → p < 0.05 = ikke-stationær

cat("\n--- Niveauer ---\n")
adf.test(na.omit(combined_dk$dk_inflation))
kpss.test(na.omit(combined_dk$dk_inflation))


cat("\n--- Første differencer ---\n")
adf.test(na.omit(combined_dk$diff_dk_inflation))
kpss.test(na.omit(combined_dk$diff_dk_inflation))


# ============================================================
# 4. Train / test split  (seneste 12 måneder som test)
# ============================================================

test_set  <- combined_dk %>% slice_tail(n = 12)
train_set <- combined_dk %>% filter(yearmonth < min(test_set$yearmonth))

cat("\nTrain:", as.character(min(train_set$yearmonth)),
    "to", as.character(max(train_set$yearmonth)),
    "(", nrow(train_set), "obs )\n")
cat("Test :", as.character(min(test_set$yearmonth)),
    "to", as.character(max(test_set$yearmonth)),
    "(", nrow(test_set), "obs )\n")


# ============================================================
# 5. Strukturelt brud dummy
#
# QLR og SIS er enige: bruddet starter oktober 2021.
# SIS finder at niveauet falder igen september 2023.
# Break-dummyen er derfor aktiv i intervallet okt 2021 – aug 2023.
# ============================================================

break_start <- yearmonth("2021 Oct")
break_end   <- yearmonth("2023 Jan")

train_models <- train_set %>%
  mutate(break_dummy = as.integer(yearmonth >= break_start & yearmonth <= break_end))

test_models <- test_set %>%
  mutate(break_dummy = as.integer(yearmonth >= break_start & yearmonth <= break_end))


# ============================================================
# 6. Fit modeller
# ============================================================

fit <- train_models %>%
  model(
    # ARIMA uden xreg
    arima_auto          = ARIMA(dk_inflation,
                                stepwise = FALSE, approximation = FALSE),
    arima_auto_seasonal = ARIMA(dk_inflation ~ PDQ(),
                                stepwise = FALSE, approximation = FALSE),
    
    # ARIMA med break-dummy
    arima_break         = ARIMA(dk_inflation ~ break_dummy,
                                stepwise = FALSE, approximation = FALSE),
    arima_break_seasonal = ARIMA(dk_inflation ~ break_dummy + PDQ(),
                                 stepwise = FALSE, approximation = FALSE),
    
    # Dynamisk regression: unemployment som xreg
    dyn_unemp           = ARIMA(dk_inflation ~ dk_unemployment,
                                stepwise = FALSE, approximation = FALSE),
    dyn_unemp_seasonal  = ARIMA(dk_inflation ~ dk_unemployment + PDQ(),
                                stepwise = FALSE, approximation = FALSE),
    
    # Dynamisk regression: unemployment + break
    dyn_unemp_break          = ARIMA(dk_inflation ~ dk_unemployment + break_dummy,
                                     stepwise = FALSE, approximation = FALSE),
    dyn_unemp_break_seasonal = ARIMA(dk_inflation ~ dk_unemployment + break_dummy + PDQ(),
                                     stepwise = FALSE, approximation = FALSE),
    
    # TSLM med break
    tslm_break        = TSLM(dk_inflation ~ trend() + break_dummy),
    tslm_break_season = TSLM(dk_inflation ~ trend() + break_dummy + season())
  )


# ============================================================
# 7. In-sample: informationskriterier
# ============================================================

fit %>%
  glance() %>%
  select(.model, AIC, AICc, BIC, sigma2) %>%
  arrange(AICc) %>%
  kable(digits = 3, align = "c", caption = "Informationskriterier")


# ============================================================
# 8. Koefficienter
# ============================================================

fit %>%
  tidy() %>%
  kable(digits = 4, align = "c", caption = "Estimerede koefficienter")



# ============================================================
# 9. Ljung-Box test — alle modeller (på træningssæt)
# ============================================================

# Antal ARMA-parametre per model bruges til dof-korrektion
# (vi sætter dof = 0 og lader fable håndtere det internt,
#  da tidy() ikke altid returnerer AR/MA-led for TSLM)

model_names <- c(
  "arima_auto", "arima_auto_seasonal",
  "arima_break", "arima_break_seasonal",
  "dyn_unemp", "dyn_unemp_seasonal",
  "dyn_unemp_break", "dyn_unemp_break_seasonal",
  "tslm_break", "tslm_break_season"
)

map_dfr(model_names, function(m) {
  fit %>%
    select(all_of(m)) %>%
    augment() %>%
    features(.innov, ljung_box, lag = 12, dof = 0) %>%
    mutate(.model = m)
}) %>%
  select(.model, lb_stat, lb_pvalue) %>%
  mutate(
    ok = ifelse(lb_pvalue > 0.05, "✓ Ingen autokor.", "✗ Autokorrelation")
  ) %>%
  arrange(lb_pvalue) %>%
  kable(digits = 4, align = "c", caption = "Ljung-Box test (lag = 24, træningssæt)")


# ============================================================
# 9. Residualdiagnostik
# ============================================================

walk(c("arima_auto", "arima_break", "dyn_unemp", "dyn_unemp_break",
       "tslm_break", "tslm_break_season"), function(m) {
         print(
           fit %>%
             select(all_of(m)) %>%
             gg_tsresiduals() +
             labs(title = paste("Residualer:", m))
         )
       })


# ============================================================
# 10. Forecast og accuracy (24-måneders test)
# ============================================================

# Modeller med xreg kræver new_data; rene ARIMA kan bruge h
fc <- bind_rows(
  fit %>% select(arima_auto)           %>% forecast(h = nrow(test_models)),
  fit %>% select(arima_auto_seasonal)  %>% forecast(h = nrow(test_models)),
  fit %>% select(arima_break)          %>% forecast(new_data = test_models),
  fit %>% select(arima_break_seasonal) %>% forecast(new_data = test_models),
  fit %>% select(dyn_unemp)            %>% forecast(new_data = test_models),
  fit %>% select(dyn_unemp_seasonal)   %>% forecast(new_data = test_models),
  fit %>% select(dyn_unemp_break)          %>% forecast(new_data = test_models),
  fit %>% select(dyn_unemp_break_seasonal) %>% forecast(new_data = test_models),
  fit %>% select(tslm_break)           %>% forecast(new_data = test_models),
  fit %>% select(tslm_break_season)    %>% forecast(new_data = test_models)
)

fc %>%
  accuracy(test_set) %>%
  select(.model, RMSE, MAE, MAPE) %>%
  arrange(RMSE) %>%
  kable(digits = 3, align = "c", caption = "Forecast accuracy — 12-måneders test")


# ============================================================
# 11. Forecast plot
# ============================================================

fc %>%
  autoplot(
    combined_dk %>% filter(yearmonth >= yearmonth("2022 Jan")),
    level = c(80, 95)
  ) +
  facet_wrap(~ .model, ncol = 2) +
  labs(
    title = "Forecasts: alle modeller",
    x = "Måned", y = "Dansk inflation (YoY, %)"
  ) +
  theme(plot.title = element_text(hjust = 0.5))










