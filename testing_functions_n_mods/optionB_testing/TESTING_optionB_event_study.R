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
