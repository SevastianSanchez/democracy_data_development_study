# =====================================================================
# TESTING ONLY: Option B quick check (not for the paper)
# ---------------------------------------------------------------------
# Question: does statistical capacity move after a country enters a
# V-Dem autocratization episode, compared with countries that never do?
#
# What this script does:
#   1. Builds autocratization onsets from V-Dem ERT (data/input_data/ert.csv)
#   2. PART A, SCI 2004-2020 (long panel, low/middle-income countries only):
#        - raw means by event time (treated vs clean controls)
#        - stacked event study (Cengiz et al. 2019 style), country-clustered SEs
#   3. PART B, SPI 2016-2023 (all incomes, short panel): same, shorter window
#   4. Optional: Callaway & Sant'Anna (2021) estimator via the `did` package
#      if it is installed (install.packages("did"))
#
# Outputs go to testing_functions_n_mods/optionB_testing/outputs/
# Only base R + dplyr/tidyr/ggplot2/sandwich/lmtest are required.
# =====================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(sandwich)
library(lmtest)

out_dir <- "testing_functions_n_mods/optionB_testing/outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------
# Settings (change here)
# ---------------------------------------------------------------------
PRE  <- 4   # years before onset shown (SCI)
POST <- 3   # years after onset shown (SCI)
REF  <- -1  # reference year (the year before onset)

# ---------------------------------------------------------------------
# 1. Autocratization onsets from V-Dem ERT
# ---------------------------------------------------------------------
ert <- read.csv("data/input_data/ert.csv") %>%
  select(country_code = country_text_id, year, aut_ep) %>%
  arrange(country_code, year) %>%
  group_by(country_code) %>%
  mutate(aut_ep = ifelse(is.na(aut_ep), 0, aut_ep),
         onset  = as.integer(aut_ep == 1 & lag(aut_ep, default = 0) == 0)) %>%
  ungroup() %>%
  filter(year >= 1990, year <= 2023)   # onsets computed on the full series first

# every onset (country, onset year g); a country can have more than one
onsets <- ert %>% filter(onset == 1) %>% transmute(country_code, g = year)

# a country-year "in an episode" flag, used to define clean controls
in_ep <- ert %>% select(country_code, year, aut_ep)

# ---------------------------------------------------------------------
# Helper: build a stacked dataset and estimate the event study
#   panel   : data frame with country_code, year, y
#   cohorts : onset years to include
#   pre/post: event window
# Treated  = countries with an onset in year g and no episode-year
#            during the pre-window (a "clean" onset)
# Controls = countries with NO episode-year anywhere in [g-pre, g+post]
# ---------------------------------------------------------------------
build_stack <- function(panel, cohorts, pre, post) {
  bind_rows(lapply(cohorts, function(g) {
    win <- (g - pre):(g + post)
    ep_win <- in_ep %>% filter(year %in% win) %>%
      group_by(country_code) %>% summarise(any_ep = any(aut_ep == 1), .groups = "drop")
    ep_pre <- in_ep %>% filter(year %in% (g - pre):(g - 1)) %>%
      group_by(country_code) %>% summarise(pre_ep = any(aut_ep == 1), .groups = "drop")

    treated  <- onsets %>% filter(g == !!g) %>%
      left_join(ep_pre, by = "country_code") %>%
      filter(is.na(pre_ep) | !pre_ep) %>% pull(country_code)
    controls <- ep_win %>% filter(!any_ep) %>% pull(country_code)

    panel %>%
      filter(year %in% win, country_code %in% c(treated, controls)) %>%
      mutate(stack = g,
             treated = as.integer(country_code %in% treated),
             event_time = year - g)
  }))
}

run_event_study <- function(stacked, pre, post, label) {
  # keep treated countries observed at the reference year
  d <- stacked %>% filter(!is.na(y))
  ets <- setdiff(-pre:post, REF)
  for (k in ets) {
    nm <- paste0("et_", ifelse(k < 0, paste0("m", abs(k)), k))
    d[[nm]] <- as.integer(d$treated == 1 & d$event_time == k)
  }
  et_vars <- grep("^et_", names(d), value = TRUE)
  d$unit   <- interaction(d$stack, d$country_code, drop = TRUE)
  d$period <- interaction(d$stack, d$year, drop = TRUE)
  f <- as.formula(paste("y ~", paste(et_vars, collapse = " + "),
                        "+ factor(unit) + factor(period)"))
  m  <- lm(f, data = d)
  vc <- vcovCL(m, cluster = ~country_code)
  ct <- coeftest(m, vcov. = vc)[et_vars, , drop = FALSE]
  res <- data.frame(term = et_vars, est = ct[, 1], se = ct[, 2]) %>%
    mutate(event_time = as.integer(sub("m", "-", sub("et_", "", term)))) %>%
    bind_rows(data.frame(term = "ref", est = 0, se = 0, event_time = REF)) %>%
    arrange(event_time) %>%
    mutate(lo = est - 1.96 * se, hi = est + 1.96 * se, sample = label)
  attr(res, "n_treated")  <- n_distinct(d$country_code[d$treated == 1])
  attr(res, "n_controls") <- n_distinct(d$country_code[d$treated == 0])
  res
}

raw_means <- function(stacked, label) {
  stacked %>% filter(!is.na(y)) %>%
    group_by(event_time, group = ifelse(treated == 1, "Autocratizing", "Never in episode (same years)")) %>%
    summarise(mean_y = mean(y), n = n(), .groups = "drop") %>%
    mutate(sample = label)
}

plot_pair <- function(raw, es, ylab, title, file) {
  p1 <- ggplot(raw, aes(event_time, mean_y, colour = group)) +
    geom_vline(xintercept = -0.5, linetype = "dashed", colour = "grey50") +
    geom_line() + geom_point() +
    labs(x = "Years since autocratization onset", y = ylab, colour = NULL,
         title = "Raw averages") +
    theme_minimal() + theme(legend.position = "bottom")
  p2 <- ggplot(es, aes(event_time, est)) +
    geom_hline(yintercept = 0, colour = "grey40") +
    geom_vline(xintercept = -0.5, linetype = "dashed", colour = "grey50") +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) +
    geom_line() + geom_point() +
    labs(x = "Years since autocratization onset",
         y = paste("Difference vs controls,", ylab),
         title = sprintf("Stacked event study (%d treated, %d controls; 95%% CI)",
                         attr(es, "n_treated"), attr(es, "n_controls"))) +
    theme_minimal()
  ggsave(file.path(out_dir, paste0(file, "_raw.png")), p1, width = 7, height = 4.5, dpi = 150)
  ggsave(file.path(out_dir, paste0(file, "_es.png")),  p2, width = 7, height = 4.5, dpi = 150)
  message("Saved: ", file, "_raw.png and ", file, "_es.png")
}

# ---------------------------------------------------------------------
# PART A: SCI 2004-2020 (low/middle-income countries only)
# ---------------------------------------------------------------------
sci <- read.csv("data/input_data/SCI_All_Dim_TS.csv", check.names = FALSE) %>%
  transmute(country_code, year = Year,
            y = suppressWarnings(as.numeric(`IQ.SCI.OVRL`))) %>%
  filter(!is.na(y))

sci_cohorts <- (min(sci$year) + PRE):(max(sci$year) - POST)  # 2008-2017
sci_stack   <- build_stack(sci, sci_cohorts, PRE, POST)
sci_es      <- run_event_study(sci_stack, PRE, POST, "SCI 2004-2020")
sci_raw     <- raw_means(sci_stack, "SCI 2004-2020")
plot_pair(sci_raw, sci_es, "SCI (0-100)",
          "SCI", "optionB_SCI")
write.csv(sci_es, file.path(out_dir, "optionB_SCI_event_study.csv"), row.names = FALSE)

# ---------------------------------------------------------------------
# PART B: SPI 2016-2023 (all income groups; short window)
# ---------------------------------------------------------------------
SPI_PRE <- 2; SPI_POST <- 3
spi <- read.csv("data/output_datasets/panel_data.csv") %>%
  select(country_code, year, y = spi_comp) %>%
  filter(!is.na(y))

spi_cohorts <- (min(spi$year) + SPI_PRE):(max(spi$year) - SPI_POST)  # 2018-2020
spi_stack   <- build_stack(spi, spi_cohorts, SPI_PRE, SPI_POST)
spi_es      <- run_event_study(spi_stack, SPI_PRE, SPI_POST, "SPI 2016-2023")
spi_raw     <- raw_means(spi_stack, "SPI 2016-2023")
plot_pair(spi_raw, spi_es, "SPI (0-100)",
          "SPI", "optionB_SPI")
write.csv(spi_es, file.path(out_dir, "optionB_SPI_event_study.csv"), row.names = FALSE)

# ---------------------------------------------------------------------
# Print a compact summary
# ---------------------------------------------------------------------
summ <- function(es) {
  cat(sprintf("  %s: %d treated, %d controls\n", unique(es$sample),
              attr(es, "n_treated"), attr(es, "n_controls")))
  print(es %>% select(event_time, est, se, lo, hi) %>% mutate(across(-event_time, ~round(.x, 2))),
        row.names = FALSE)
}
cat("\n=== Option B quick check ===\n"); summ(sci_es); summ(spi_es)

# ---------------------------------------------------------------------
# OPTIONAL: Callaway & Sant'Anna (2021), if `did` is installed
# ---------------------------------------------------------------------
if (requireNamespace("did", quietly = TRUE)) {
  # first onset inside the SCI window; countries already in an episode at
  # the start of the window are dropped (neither clean treated nor control)
  start_ep <- in_ep %>% filter(year == min(sci$year), aut_ep == 1) %>% pull(country_code)
  g_first  <- onsets %>% filter(g > min(sci$year), g <= max(sci$year)) %>%
    group_by(country_code) %>% summarise(g = min(g), .groups = "drop")
  cs_dat <- sci %>%
    filter(!country_code %in% start_ep) %>%
    left_join(g_first, by = "country_code") %>%
    mutate(g = ifelse(is.na(g), 0, g), id = as.integer(factor(country_code)))
  cs <- did::att_gt(yname = "y", tname = "year", idname = "id", gname = "g",
                    data = cs_dat, control_group = "nevertreated",
                    base_period = "universal")
  cs_es <- did::aggte(cs, type = "dynamic", min_e = -PRE, max_e = POST, na.rm = TRUE)
  print(summary(cs_es))
  ggsave(file.path(out_dir, "optionB_SCI_callaway_santanna.png"),
         did::ggdid(cs_es), width = 7, height = 4.5, dpi = 150)
} else {
  message("`did` not installed: skipped Callaway & Sant'Anna check.")
}

# =====================================================================
# PART C: SPI Pillar 2 (data services) vs Pillar 5 (data infrastructure)
# ---------------------------------------------------------------------
# Mechanism test: if erosion works through dissemination, Pillar 2 should
# fall after onset while Pillar 5 should not.
# Data: World Bank SPI_index.csv (github.com/worldbank/SPI, 03_output_data),
# saved here as SPI_index_worldbank_github.csv. Pillar scores 2016-2024.
# Also: Pillar 2.2 online access (ODIN openness, 2015-2024).
# =====================================================================
spi_raw <- read.csv("testing_functions_n_mods/optionB_testing/SPI_index_worldbank_github.csv",
                    fileEncoding = "UTF-8-BOM") %>%
  rename(country_code = iso3c, year = date)

run_outcome <- function(var, label, pre, post, file) {
  pan <- spi_raw %>% transmute(country_code, year, y = .data[[var]]) %>% filter(!is.na(y))
  coh <- (min(pan$year) + pre):(max(pan$year) - post)
  st  <- build_stack(pan, coh, pre, post)
  es  <- run_event_study(st, pre, post, label)
  plot_pair(raw_means(st, label), es, paste(label, "(0-100)"), label, file)
  es
}

p2_es  <- run_outcome("SPI.INDEX.PIL2",   "Pillar 2 data services",       2, 3, "optionB_SPI_P2")
p5_es  <- run_outcome("SPI.INDEX.PIL5",   "Pillar 5 data infrastructure", 2, 3, "optionB_SPI_P5")
p22_es <- run_outcome("SPI.DIM2.2.INDEX", "Pillar 2.2 online access (ODIN)", 3, 3, "optionB_SPI_P2_2")

# side-by-side plot of P2 vs P5
both <- bind_rows(p2_es, p5_es)
p_both <- ggplot(both, aes(event_time, est, colour = sample, fill = sample)) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_vline(xintercept = -0.5, linetype = "dashed", colour = "grey50") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line() + geom_point() +
  labs(x = "Years since autocratization onset", y = "Difference vs controls (0-100)",
       colour = NULL, fill = NULL,
       title = sprintf("SPI Pillar 2 vs Pillar 5 (%d treated, %d controls; 95%% CI)",
                       attr(p2_es, "n_treated"), attr(p2_es, "n_controls"))) +
  theme_minimal() + theme(legend.position = "bottom")
ggsave(file.path(out_dir, "optionB_SPI_P2_vs_P5.png"), p_both, width = 7, height = 4.5, dpi = 150)

cat("\n=== Part C: Pillar 2 vs Pillar 5 ===\n"); summ(p2_es); summ(p5_es); summ(p22_es)

# =====================================================================
# PART D: all five SPI pillars, same window and sample (2016-2024)
# ---------------------------------------------------------------------
# Same stacked design as Part C for every pillar, so the panels are
# directly comparable. Pillars 1 and 3 exist back to 2004/2005 in the
# World Bank file, but are restricted to 2016+ here on purpose.
# =====================================================================
pillars <- c(SPI.INDEX.PIL1 = "Pillar 1: Data use",
             SPI.INDEX.PIL2 = "Pillar 2: Data services",
             SPI.INDEX.PIL3 = "Pillar 3: Data products",
             SPI.INDEX.PIL4 = "Pillar 4: Data sources",
             SPI.INDEX.PIL5 = "Pillar 5: Data infrastructure")

pillar_es <- bind_rows(lapply(names(pillars), function(v) {
  pan <- spi_raw %>% transmute(country_code, year, y = .data[[v]]) %>%
    filter(!is.na(y), year >= 2016)
  st  <- build_stack(pan, 2018:2021, 2, 3)
  es  <- run_event_study(st, 2, 3, pillars[[v]])
  es$n_treated <- attr(es, "n_treated"); es$n_controls <- attr(es, "n_controls")
  es
})) %>%
  mutate(facet = sprintf("%s\n(%d treated, %d controls)", sample, n_treated, n_controls))

p_all <- ggplot(pillar_es, aes(event_time, est)) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_vline(xintercept = -0.5, linetype = "dashed", colour = "grey50") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) +
  geom_line() + geom_point() +
  facet_wrap(~facet, ncol = 3) +
  labs(x = "Years since autocratization onset",
       y = "Difference vs never-autocratizing countries (0-100)",
       title = "SPI pillars around V-Dem autocratization onsets (stacked event study, 95% CI)",
       subtitle = "SPI 2016-2024; onset cohorts 2018-2021; reference year = -1; SEs clustered by country") +
  theme_minimal(base_size = 10)
ggsave(file.path(out_dir, "optionB_SPI_all_pillars.png"), p_all, width = 10, height = 6.5, dpi = 150)
write.csv(pillar_es, file.path(out_dir, "optionB_SPI_all_pillars.csv"), row.names = FALSE)

cat("\n=== Part D: all pillars ===\n")
print(pillar_es %>% select(sample, event_time, est, se, lo, hi) %>%
        mutate(across(c(est, se, lo, hi), ~round(.x, 2))), row.names = FALSE)

# =====================================================================
# PART E: parallel-trends checks for the five pillars
# ---------------------------------------------------------------------
# 1. Raw averages, treated vs controls, by event time (visual check)
# 2. Pre-trend test: are the pre-onset coefficients jointly zero?
#    (Wald test with country-clustered SEs)
# 3. Longer pre-period: 3 years before onset (cohorts 2019-2021 only)
# =====================================================================
pretest <- function(pan, cohorts, pre, post, label) {
  st <- build_stack(pan, cohorts, pre, post) %>% filter(!is.na(y))
  ets <- setdiff(-pre:post, REF)
  for (k in ets) st[[paste0("et_", ifelse(k < 0, paste0("m", abs(k)), k))]] <-
    as.integer(st$treated == 1 & st$event_time == k)
  etv <- grep("^et_", names(st), value = TRUE)
  st$unit <- interaction(st$stack, st$country_code, drop = TRUE)
  st$period <- interaction(st$stack, st$year, drop = TRUE)
  m  <- lm(as.formula(paste("y ~", paste(etv, collapse = "+"), "+ factor(unit) + factor(period)")), data = st)
  vc <- vcovCL(m, cluster = ~country_code)
  leads <- grep("^et_m", etv, value = TRUE)
  b <- coef(m)[leads]; V <- vc[leads, leads, drop = FALSE]
  W <- as.numeric(t(b) %*% solve(V) %*% b)
  data.frame(outcome = label, pre_years = pre,
             n_treated = n_distinct(st$country_code[st$treated == 1]),
             leads_tested = paste(sub("et_m", "-", leads), collapse = ","),
             wald_chi2 = round(W, 2), df = length(leads),
             p_value = round(pchisq(W, length(leads), lower.tail = FALSE), 3))
}

pt_tab <- list(); raw_all <- list(); es_long <- list()
for (v in names(pillars)) {
  pan <- spi_raw %>% transmute(country_code, year, y = .data[[v]]) %>% filter(!is.na(y), year >= 2016)
  pt_tab[[v]] <- bind_rows(pretest(pan, 2018:2021, 2, 3, pillars[[v]]),
                           pretest(pan, 2019:2021, 3, 3, pillars[[v]]))
  raw_all[[v]] <- raw_means(build_stack(pan, 2018:2021, 2, 3), pillars[[v]])
  es_long[[v]] <- run_event_study(build_stack(pan, 2019:2021, 3, 3), 3, 3, pillars[[v]])
}
pt_tab <- bind_rows(pt_tab); raw_all <- bind_rows(raw_all); es_long <- bind_rows(es_long)

p_raw <- ggplot(raw_all, aes(event_time, mean_y, colour = group)) +
  geom_vline(xintercept = -0.5, linetype = "dashed", colour = "grey50") +
  geom_line() + geom_point(size = 1.2) +
  facet_wrap(~sample, ncol = 3, scales = "free_y") +
  labs(x = "Years since autocratization onset", y = "Average pillar score (0-100)", colour = NULL,
       title = "Raw averages: autocratizing vs never-autocratizing countries",
       subtitle = "Parallel trends: the two lines should move together before the dashed line") +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(out_dir, "optionB_SPI_all_pillars_raw.png"), p_raw, width = 10, height = 6.5, dpi = 150)

p_long <- ggplot(es_long, aes(event_time, est)) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_vline(xintercept = -0.5, linetype = "dashed", colour = "grey50") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) + geom_line() + geom_point() +
  facet_wrap(~sample, ncol = 3) +
  labs(x = "Years since autocratization onset", y = "Difference vs controls (0-100)",
       title = "Longer pre-period: 3 years before onset (cohorts 2019-2021)") +
  theme_minimal(base_size = 10)
ggsave(file.path(out_dir, "optionB_SPI_all_pillars_pre3.png"), p_long, width = 10, height = 6.5, dpi = 150)
write.csv(pt_tab, file.path(out_dir, "optionB_pretrend_tests.csv"), row.names = FALSE)

cat("\n=== Part E: pre-trend tests ===\n"); print(pt_tab, row.names = FALSE)
cat("\nLonger pre-period estimates:\n")
print(es_long %>% select(sample, event_time, est, lo, hi) %>% mutate(across(c(est, lo, hi), ~round(.x, 2))), row.names = FALSE)

# =====================================================================
# PART F: SCI sub-dimensions, 2004-2020 (low/middle-income countries)
# ---------------------------------------------------------------------
# Longer-run version of the dissemination vs. other-capacity contrast.
# SCI has three dimensions:
#   IQ.SCI.PRDC = periodicity and timeliness (closest to dissemination)
#   IQ.SCI.MTHD = statistical methodology (standards, classifications)
#   IQ.SCI.SRCE = source data (censuses, surveys, vital registration)
# =====================================================================
sci_dims <- c(IQ.SCI.PRDC = "SCI: Periodicity & timeliness",
              IQ.SCI.MTHD = "SCI: Methodology",
              IQ.SCI.SRCE = "SCI: Source data",
              IQ.SCI.OVRL = "SCI: Overall")
sci_all <- read.csv("data/input_data/SCI_All_Dim_TS.csv", check.names = FALSE)

sci_dim_es <- list(); sci_dim_raw <- list(); sci_dim_pt <- list()
for (v in names(sci_dims)) {
  pan <- sci_all %>% transmute(country_code, year = Year,
                               y = suppressWarnings(as.numeric(.data[[v]]))) %>% filter(!is.na(y))
  coh <- (min(pan$year) + PRE):(max(pan$year) - POST)
  st  <- build_stack(pan, coh, PRE, POST)
  sci_dim_es[[v]]  <- run_event_study(st, PRE, POST, sci_dims[[v]]) %>%
    mutate(n_treated = n_distinct(st$country_code[st$treated == 1]),
           n_controls = n_distinct(st$country_code[st$treated == 0]))
  sci_dim_raw[[v]] <- raw_means(st, sci_dims[[v]])
  sci_dim_pt[[v]]  <- pretest(pan, coh, PRE, POST, sci_dims[[v]])
}
sci_dim_es <- bind_rows(sci_dim_es) %>%
  mutate(facet = sprintf("%s\n(%d treated, %d controls)", sample, n_treated, n_controls))
sci_dim_raw <- bind_rows(sci_dim_raw); sci_dim_pt <- bind_rows(sci_dim_pt)

p_sci_es <- ggplot(sci_dim_es, aes(event_time, est)) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_vline(xintercept = -0.5, linetype = "dashed", colour = "grey50") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) + geom_line() + geom_point() +
  facet_wrap(~facet, ncol = 2) +
  labs(x = "Years since autocratization onset", y = "Difference vs controls (0-100)",
       title = "SCI dimensions around V-Dem autocratization onsets (2004-2020, 95% CI)",
       subtitle = "Stacked event study; onset cohorts 2008-2017; low/middle-income countries") +
  theme_minimal(base_size = 10)
p_sci_raw <- ggplot(sci_dim_raw, aes(event_time, mean_y, colour = group)) +
  geom_vline(xintercept = -0.5, linetype = "dashed", colour = "grey50") +
  geom_line() + geom_point(size = 1.2) + facet_wrap(~sample, ncol = 2, scales = "free_y") +
  labs(x = "Years since autocratization onset", y = "Average score (0-100)", colour = NULL,
       title = "SCI dimensions: raw averages (parallel-trends check)") +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(out_dir, "optionB_SCI_dimensions_es.png"),  p_sci_es,  width = 9, height = 6.5, dpi = 150)
ggsave(file.path(out_dir, "optionB_SCI_dimensions_raw.png"), p_sci_raw, width = 9, height = 6.5, dpi = 150)

cat("\n=== Part F: SCI dimensions ===\n"); print(sci_dim_pt, row.names = FALSE)
print(sci_dim_es %>% select(sample, event_time, est, lo, hi) %>% mutate(across(c(est, lo, hi), ~round(.x, 2))), row.names = FALSE)

# =====================================================================
# PART G: HRV transparency index (1993-2015, all income groups)
# ---------------------------------------------------------------------
# Dissemination proxy for the pre-2016 period (Hollyer, Rosendorff &
# Vreeland; HRV2 update 2024). Download from QoG:
#   https://datafinder.qog.gu.se/dataset/hrv  (variable hrv_index)
# Save as testing_functions_n_mods/optionB_testing/hrv.csv with columns
# ccodealp (ISO3), year, hrv_index. Skipped if the file is missing.
# Current hrv.csv = HRV v1 (1980-2010, 125 countries) extracted from the
# QoG Standard Time-Series dataset, Jan 2026 edition.
# =====================================================================
hrv_path <- "testing_functions_n_mods/optionB_testing/hrv.csv"
if (file.exists(hrv_path)) {
  hrv <- read.csv(hrv_path) %>% transmute(country_code = ccodealp, year, y = hrv_index) %>%
    filter(!is.na(y), year >= 1990)
  coh <- (min(hrv$year) + PRE):(max(hrv$year) - POST)
  st  <- build_stack(hrv, coh, PRE, POST)
  hrv_es <- run_event_study(st, PRE, POST, "HRV transparency index")
  plot_pair(raw_means(st, "HRV transparency index"), hrv_es, "HRV index", "HRV", "optionB_HRV")
  print(pretest(hrv, coh, PRE, POST, "HRV transparency index"), row.names = FALSE)
  summ(hrv_es)
} else message("hrv.csv not found: skipped Part G (see header for download).")
