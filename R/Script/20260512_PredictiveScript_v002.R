# ============================================================
# 0. Setup
# ============================================================
# rm(list = ls())
# dev.off()

library(tidyverse)
library(lubridate)
library(knitr)
library(patchwork)
library(fpp3)
library(tseries)
library(strucchange)
library(gets)
library(vars)
library(dplyr)
library(conflicted)

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("lag", "dplyr")


# ============================================================
# 1. Indlæs data
# ============================================================

find_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)))
  }
  
  frame_files <- vapply(sys.frames(), function(x) {
    ofile <- x$ofile
    if (is.null(ofile)) "" else ofile
  }, character(1))
  
  frame_files <- frame_files[nzchar(frame_files)]
  
  if (length(frame_files) > 0) {
    return(dirname(normalizePath(frame_files[length(frame_files)], mustWork = FALSE)))
  }
  
  getwd()
}

find_data_dir <- function(start_dir) {
  candidate_names <- c("Data - csv.", "Data - csv")
  current_dir <- normalizePath(start_dir, mustWork = FALSE)
  
  repeat {
    existing <- file.path(current_dir, candidate_names)[
      dir.exists(file.path(current_dir, candidate_names))
    ]
    
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
  rename(
    dk_inflation = DK,
    eu_inflation = EU27_2020
  )

arbejdsloes_wide <- laes_eurostat_csv(
  file.path(data_dir, "une_rt_m__custom_21465700_linear.csv")
) |>
  pivot_wider(names_from = geo, values_from = vaerdi) |>
  rename(
    dk_unemployment = DK,
    eu_unemployment = EU27_2020
  )

df <- full_join(inflation_wide, arbejdsloes_wide, by = "dato") |>
  arrange(dato) |>
  drop_na()

cat("Rows:", nrow(df), "\n")
cat("Period:", format(min(df$dato), "%b %Y"),
    "to", format(max(df$dato), "%b %Y"), "\n")


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

cat("\n--- Inflation i niveauer ---\n")
adf.test(na.omit(combined_dk$dk_inflation))
kpss.test(na.omit(combined_dk$dk_inflation))

cat("\n--- Inflation i første differencer ---\n")
adf.test(na.omit(combined_dk$diff_dk_inflation))
kpss.test(na.omit(combined_dk$diff_dk_inflation))


# ============================================================
# 4. Strukturelle brudtest
# ============================================================

# Vi tester først for et mean shift i dansk inflation.
# Derefter tester vi for brud i relationen:
#
#   dk_inflation ~ dk_unemployment
#
# Til ARIMA-modellen nedenfor bruger vi én pre/post-break dummy:
#
#   post_break = 0 før breakdatoen
#   post_break = 1 fra breakdatoen og frem
#
# Det er en klassisk step dummy / level shift dummy.

qlr_data <- combined_dk |>
  as_tibble() |>
  transmute(
    yearmonth,
    dato = as.Date(yearmonth),
    dk_inflation,
    dk_unemployment
  ) |>
  drop_na()

dk_infl_ts <- ts(
  qlr_data$dk_inflation,
  start     = c(year(min(qlr_data$yearmonth)),
                month(min(qlr_data$yearmonth))),
  frequency = 12
)


# ------------------------------------------------------------
# 4.1 QLR / supF: mean shift i dansk inflation
# ------------------------------------------------------------

cat("\n=== QLR supF test: Mean shift i dansk inflation ===\n")

qlr_mean <- Fstats(
  dk_infl_ts ~ 1,
  from = 0.15,
  to   = 0.85
)

qlr_mean_test <- sctest(qlr_mean, type = "supF")
print(qlr_mean_test)

bp_mean_1 <- breakpoints(dk_infl_ts ~ 1, breaks = 1)

summary(bp_mean_1)

bp_mean_1_date <- qlr_data$yearmonth[bp_mean_1$breakpoints]

cat("Estimeret brudsdato, mean shift:",
    as.character(bp_mean_1_date), "\n")

plot(
  qlr_mean,
  main = "QLR Test: Mean Shift i dansk inflation",
  xlab = "Tid",
  ylab = "F-statistik"
)
lines(bp_mean_1, col = "red")


# ------------------------------------------------------------
# 4.2 QLR / supF: brud i inflation-arbejdsløshed relation
# ------------------------------------------------------------

cat("\n=== QLR supF test: DK inflation ~ DK unemployment ===\n")

qlr_formula <- dk_inflation ~ dk_unemployment

qlr_reg <- Fstats(
  qlr_formula,
  data = qlr_data,
  from = 0.15,
  to   = 0.85
)

qlr_reg_test <- sctest(qlr_reg, type = "supF")
print(qlr_reg_test)

bp_reg_1 <- breakpoints(
  qlr_formula,
  data = qlr_data,
  breaks = 1
)

summary(bp_reg_1)

bp_reg_1_date <- qlr_data$yearmonth[bp_reg_1$breakpoints]

cat("Estimeret brudsdato, regression:",
    as.character(bp_reg_1_date), "\n")

plot(
  qlr_reg,
  main = "QLR Test: Brud i inflation-arbejdsløshed relation",
  xlab = "Tid",
  ylab = "F-statistik"
)
lines(bp_reg_1, col = "red")


# ------------------------------------------------------------
# 4.3 Pre/post-break regressioner
# ------------------------------------------------------------

cat("\n--- Phillips-kurve FØR regression-break ---\n")
summary(
  lm(
    dk_inflation ~ dk_unemployment,
    data = qlr_data,
    subset = yearmonth < bp_reg_1_date
  )
)

cat("\n--- Phillips-kurve EFTER regression-break ---\n")
summary(
  lm(
    dk_inflation ~ dk_unemployment,
    data = qlr_data,
    subset = yearmonth >= bp_reg_1_date
  )
)


# ------------------------------------------------------------
# 4.4 Bai-Perron: multiple brud
# ------------------------------------------------------------

cat("\n=== Bai-Perron: multiple brud i regressionen ===\n")

bp_reg_multi <- breakpoints(
  qlr_formula,
  data = qlr_data,
  breaks = 4
)

summary(bp_reg_multi)
confint(bp_reg_multi)

bp_reg_multi_dates <- qlr_data$yearmonth[bp_reg_multi$breakpoints]

cat("Estimerede brudsdatoer:",
    paste(as.character(bp_reg_multi_dates), collapse = ", "), "\n")

qlr_regimes <- qlr_data |>
  mutate(
    regime = cut(
      row_number(),
      breaks = c(0, bp_reg_multi$breakpoints, n()),
      labels = paste("Regime", seq_along(c(bp_reg_multi$breakpoints, n()))),
      include.lowest = TRUE
    )
  )

qlr_regimes |>
  group_by(regime) |>
  summarise(
    start              = as.character(min(yearmonth)),
    end                = as.character(max(yearmonth)),
    n                  = n(),
    mean_inflation     = mean(dk_inflation, na.rm = TRUE),
    mean_unemployment  = mean(dk_unemployment, na.rm = TRUE),
    unemployment_slope = coef(lm(dk_inflation ~ dk_unemployment))[2],
    .groups = "drop"
  ) |>
  kable(
    digits = 3,
    align = "c",
    caption = "Regime-oversigt: Bai-Perron"
  )

ggplot(qlr_regimes, aes(x = dato, y = dk_inflation, colour = regime)) +
  geom_line(linewidth = 0.8) +
  geom_vline(
    xintercept = as.Date(bp_reg_multi_dates),
    linetype = "dashed",
    colour = "red"
  ) +
  labs(
    title = "Dansk inflation: strukturelle regimer",
    x = "Måned",
    y = "Inflation (YoY, %)",
    colour = NULL
  ) +
  theme(plot.title = element_text(hjust = 0.5))


# ------------------------------------------------------------
# 4.5 Valg af breakdato til ARIMA-modeller
# ------------------------------------------------------------

# Til ARIMA-modellen bruger vi mean-shift breaket, fordi ARIMA-modellen
# nedenfor handler om inflationens niveau, ikke om regressionsrelationen
# til arbejdsløshed.
#
# Hvis du hellere vil bruge regression-breaket, kan du ændre linjen til:
#   break_date <- bp_reg_1_date

break_date <- bp_mean_1_date

cat("\nBreakdato anvendt i ARIMA-modeller:",
    as.character(break_date), "\n")


# ============================================================
# 5. Train / test split
# ============================================================

# Vi bruger seneste 12 måneder som test.
# Skift slice_tail(n = 12) til n = 24, hvis du ønsker 2 års test.

test_set <- combined_dk %>%
  slice_tail(n = 12)

train_set <- combined_dk %>%
  filter(yearmonth < min(test_set$yearmonth))

cat("\nTrain:", as.character(min(train_set$yearmonth)),
    "to", as.character(max(train_set$yearmonth)),
    "(", nrow(train_set), "obs )\n")

cat("Test :", as.character(min(test_set$yearmonth)),
    "to", as.character(max(test_set$yearmonth)),
    "(", nrow(test_set), "obs )\n")

# ============================================================
# 6. Pre/post-break dummies: to structural breaks
# ============================================================

# Vi bruger de to breakdatoer fra Bai-Perron-testen.
# Hvis bp_reg_multi_dates indeholder flere end to breaks, tager vi de to første.
# Alternativt kan du hardcode dem:
# break_date_1 <- yearmonth("2013 Jan")
# break_date_2 <- yearmonth("2021 Dec")

break_dates <- sort(bp_reg_multi_dates)
break_date_1 <- break_dates[1]
break_date_2 <- break_dates[2]


cat("\nBreakdato 1 anvendt i ARIMA-modeller:",
    as.character(break_date_1), "\n")

cat("Breakdato 2 anvendt i ARIMA-modeller:",
    as.character(break_date_2), "\n")


# ------------------------------------------------------------
# Lav to post-break dummies
# ------------------------------------------------------------

train_models <- train_set %>%
  mutate(
    post_break_1 = as.integer(yearmonth >= break_date_1),
    post_break_2 = as.integer(yearmonth >= break_date_2)
  )

test_models <- test_set %>%
  mutate(
    post_break_1 = as.integer(yearmonth >= break_date_1),
    post_break_2 = as.integer(yearmonth >= break_date_2)
  )


# ------------------------------------------------------------
# Tjek dummyerne omkring begge breakdatoer
# ------------------------------------------------------------

combined_dk %>%
  mutate(
    post_break_1 = as.integer(yearmonth >= break_date_1),
    post_break_2 = as.integer(yearmonth >= break_date_2)
  ) %>%
  filter(
    (yearmonth >= break_date_1 - 3 & yearmonth <= break_date_1 + 3) |
      (yearmonth >= break_date_2 - 3 & yearmonth <= break_date_2 + 3)
  ) %>%
  select(yearmonth, dk_inflation, post_break_1, post_break_2) %>%
  kable(
    digits = 3,
    align = "c",
    caption = "Tjek af to post-break dummies"
  )


# ============================================================
# 7. ARIMA-modeltilpasning: ignorerer vs. håndterer to breaks
# ============================================================

# Vi sammenligner:
#
#   1. arima_ignore_break:
#      Auto ARIMA uden structural break-dummies.
#
#   2. arima_handle_break:
#      Auto ARIMA med to post-break dummies.
#
# Fortolkning:
#   post_break_1 måler niveauskift fra break 1 og frem.
#   post_break_2 måler yderligere niveauskift fra break 2 og frem.
#
# Dermed får modellen tre implicitte regimer:
#   Regime 1: før break 1
#   Regime 2: mellem break 1 og break 2
#   Regime 3: efter break 2

fit_arima <- train_models %>%
  model(
    arima_ignore_break = ARIMA(
      dk_inflation,
      stepwise = FALSE,
      approximation = FALSE
    ),
    
    arima_handle_break = ARIMA(
      dk_inflation ~ post_break_1 + post_break_2,
      stepwise = FALSE,
      approximation = FALSE
    )
  )


# ------------------------------------------------------------
# 7.1 Se hvilke ARIMA-modeller der er valgt
# ------------------------------------------------------------

fit_arima

fit_arima %>%
  select(arima_ignore_break) %>%
  report()

fit_arima %>%
  select(arima_handle_break) %>%
  report()


# ------------------------------------------------------------
# 7.2 Information criteria
# ------------------------------------------------------------

fit_arima %>%
  glance() %>%
  select(.model, AIC, AICc, BIC, sigma2) %>%
  arrange(AICc) %>%
  kable(
    digits = 3,
    align = "c",
    caption = "ARIMA: ignorerer vs. håndterer to structural breaks"
  )


# ------------------------------------------------------------
# 7.3 Koefficienter
# ------------------------------------------------------------

fit_arima %>%
  tidy() %>%
  kable(
    digits = 4,
    align = "c",
    caption = "Koefficienter: ARIMA med og uden to breaks"
  )


# ------------------------------------------------------------
# 7.4 Ljung-Box residualtest
# ------------------------------------------------------------

ljung_box_arima <- function(model_table, model_name, lag_value = 12) {
  
  dof_m <- model_table %>%
    select(all_of(model_name)) %>%
    tidy() %>%
    filter(str_detect(term, "^ar|^ma|^sar|^sma")) %>%
    nrow()
  
  model_table %>%
    select(all_of(model_name)) %>%
    augment() %>%
    features(
      .innov,
      ljung_box,
      lag = lag_value,
      dof = dof_m
    ) %>%
    mutate(
      .model = model_name,
      dof = dof_m,
      passes_ljung_box = lb_pvalue > 0.05
    )
}

bind_rows(
  ljung_box_arima(fit_arima, "arima_ignore_break"),
  ljung_box_arima(fit_arima, "arima_handle_break")
) %>%
  select(.model, lb_stat, lb_pvalue, dof, passes_ljung_box) %>%
  arrange(lb_pvalue) %>%
  kable(
    digits = 4,
    align = "c",
    caption = "Ljung-Box residualtest: ARIMA med og uden to breaks"
  )


# ------------------------------------------------------------
# 7.5 Residualplots
# ------------------------------------------------------------

fit_arima %>%
  select(arima_ignore_break) %>%
  gg_tsresiduals() +
  labs(title = "Residualdiagnostik: ARIMA ignorerer breaks")

fit_arima %>%
  select(arima_handle_break) %>%
  gg_tsresiduals() +
  labs(title = "Residualdiagnostik: ARIMA med to post-break dummies")



# ============================================================
# 8. ARIMA på første differencer:
# ignorerer vs. håndterer to structural breaks
# ============================================================

train_diff <- train_models %>%
  mutate(
    diff_dk_inflation = difference(dk_inflation)
  ) %>%
  filter(!is.na(diff_dk_inflation))

test_diff <- test_models %>%
  mutate(
    diff_dk_inflation = difference(dk_inflation)
  ) %>%
  filter(!is.na(diff_dk_inflation))


fit_arima_diff <- train_diff %>%
  model(
    arima_diff_ignore_break = ARIMA(
      diff_dk_inflation,
      stepwise = FALSE,
      approximation = FALSE
    ),
    
    arima_diff_handle_break = ARIMA(
      diff_dk_inflation ~ post_break_1 + post_break_2,
      stepwise = FALSE,
      approximation = FALSE
    )
  )


# ------------------------------------------------------------
# 8.1 Se modeller
# ------------------------------------------------------------

fit_arima_diff

fit_arima_diff %>%
  select(arima_diff_ignore_break) %>%
  report()

fit_arima_diff %>%
  select(arima_diff_handle_break) %>%
  report()


# ------------------------------------------------------------
# 8.2 Information criteria
# ------------------------------------------------------------

fit_arima_diff %>%
  glance() %>%
  select(.model, AIC, AICc, BIC, sigma2) %>%
  arrange(AICc) %>%
  kable(
    digits = 3,
    align = "c",
    caption = "ARIMA på ∆inflation: med og uden to breaks"
  )


# ------------------------------------------------------------
# 8.3 Koefficienter
# ------------------------------------------------------------

fit_arima_diff %>%
  tidy() %>%
  kable(
    digits = 4,
    align = "c",
    caption = "Koefficienter: ARIMA på ∆inflation med og uden to breaks"
  )


# ------------------------------------------------------------
# 8.4 Ljung-Box residualtest
# ------------------------------------------------------------

bind_rows(
  ljung_box_arima(fit_arima_diff, "arima_diff_ignore_break"),
  ljung_box_arima(fit_arima_diff, "arima_diff_handle_break")
) %>%
  select(.model, lb_stat, lb_pvalue, dof, passes_ljung_box) %>%
  arrange(lb_pvalue) %>%
  kable(
    digits = 4,
    align = "c",
    caption = "Ljung-Box: ARIMA på ∆inflation med og uden to breaks"
  )


# ------------------------------------------------------------
# 8.5 Residualplots
# ------------------------------------------------------------

fit_arima_diff %>%
  select(arima_diff_ignore_break) %>%
  gg_tsresiduals() +
  labs(title = "Residualdiagnostik: ARIMA på ∆inflation uden breaks")

fit_arima_diff %>%
  select(arima_diff_handle_break) %>%
  gg_tsresiduals() +
  labs(title = "Residualdiagnostik: ARIMA på ∆inflation med to post-break dummies")