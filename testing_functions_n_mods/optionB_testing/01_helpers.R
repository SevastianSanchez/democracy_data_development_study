# =====================================================================
# 01_helpers.R  |  Option B preliminary tests: settings and functions
# ---------------------------------------------------------------------
# PURPOSE
#   Shared setup for 02_run_tests.R. Nothing here produces output on its
#   own. It (1) loads packages and sets folders, (2) builds the list of
#   V-Dem autocratization onsets (the "treatment"), and (3) defines the
#   functions every test uses: build the stacked sample, estimate the
#   event study, test pre-trends, and draw labeled plots.
#
# RUN FROM: the repository root (sourced by 02_run_tests.R)
#
# DESIGN IN ONE PARAGRAPH (stacked event study, Cengiz et al. 2019)
#   For every onset year g, we build a small panel covering the years
#   g-pre ... g+post. "Treated" countries start a V-Dem autocratization
#   episode in year g and were not in an episode during the pre-window.
#   "Clean controls" are countries that are in no episode at any point
#   in that window. The stacks are appended and we regress the outcome
#   on treated x event-time dummies, with stack-by-country and
#   stack-by-year fixed effects. The year before onset (-1) is the
#   reference, so every coefficient reads as "gap between treated and
#   controls in year k, compared with the gap in the year before onset".
#   Standard errors are clustered by country.
# =====================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(sandwich)    # clustered standard errors
  library(lmtest)      # coeftest()
  library(patchwork)   # put the raw-trends and event-study panels side by side
})

# ---------------------------------------------------------------------
# 1. Folders and global settings
# ---------------------------------------------------------------------
STUDY_DIR <- "testing_functions_n_mods/optionB_testing"
DATA_DIR  <- file.path(STUDY_DIR, "data")
FIG_DIR   <- file.path(STUDY_DIR, "outputs", "figures")
TAB_DIR   <- file.path(STUDY_DIR, "outputs", "tables")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)

REF <- -1   # reference year: the year before onset (fixed at 0 in every plot)

# Plot colours (validated colour-blind-safe pair) and shared caption text
COL_TREATED <- "#2a78d6"   # blue   = autocratizing countries / first series
COL_CONTROL <- "#eb6834"   # orange = clean controls / second series
GROUP_COLS  <- c("Autocratizing (treated)" = COL_TREATED,
                 "Clean controls (no episode in window)" = COL_CONTROL)

CAPTION_ES <- paste(
  "How to read the event study: each point is the gap between autocratizing countries and clean controls",
  "in that year, minus the gap in the year before onset (-1, fixed at 0).",
  "Shaded band = 95% confidence interval, standard errors clustered by country.",
  "Dashed line = onset of a V-Dem autocratization episode.",
  "Points near 0 before the dashed line support parallel trends; points below 0 after it mean the outcome fell",
  "relative to controls.", sep = "\n")

CAPTION_RAW <- paste(
  "How to read the raw trends: plain averages of the outcome by year relative to onset, before any model.",
  "Parallel-trends check: the two lines should move together before the dashed line.", sep = "\n")

# ---------------------------------------------------------------------
# 2. Treatment: V-Dem autocratization onsets (ERT, aut_ep)
# ---------------------------------------------------------------------
# onset = first year of an autocratization episode (aut_ep switches 0 -> 1).
# Onsets are computed on the full 1900-2024 series BEFORE keeping 1990+,
# so a country already mid-episode in 1990 is not counted as a new onset.
vdem <- read.csv(file.path(DATA_DIR, "vdem_ert_episodes.csv"))

ert <- vdem %>%
  select(country_code, year, aut_ep, v2x_polyarchy) %>%
  arrange(country_code, year) %>%
  group_by(country_code) %>%
  mutate(aut_ep = ifelse(is.na(aut_ep), 0, aut_ep),
         onset  = as.integer(aut_ep == 1 & lag(aut_ep, default = 0) == 0)) %>%
  ungroup() %>%
  filter(year >= 1990)

# every onset (country, onset year g); a country can have more than one
onsets <- ert %>% filter(onset == 1) %>% transmute(country_code, g = year)

# country-year "in an episode" flag, used to define clean treated and controls
in_ep <- ert %>% select(country_code, year, aut_ep)

# ---------------------------------------------------------------------
# 3. Functions
# ---------------------------------------------------------------------

# make_panel(): pick one outcome column from a data file and return the
# three columns every test uses (country_code, year, y).
make_panel <- function(df, var, from = -Inf, to = Inf) {
  df %>%
    transmute(country_code, year, y = suppressWarnings(as.numeric(.data[[var]]))) %>%
    filter(!is.na(y), year >= from, year <= to)
}

# build_stack(): one stacked dataset with a block for each onset year g.
#   treated  = onset in g, and no episode-year in g-pre ... g-1
#   controls = no episode-year anywhere in g-pre ... g+post
build_stack <- function(panel, cohorts, pre, post) {
  bind_rows(lapply(cohorts, function(g) {
    win <- (g - pre):(g + post)
    ep_win <- in_ep %>% filter(year %in% win) %>%
      group_by(country_code) %>% summarise(any_ep = any(aut_ep == 1), .groups = "drop")
    ep_pre <- in_ep %>% filter(year %in% (g - pre):(g - 1)) %>%
      group_by(country_code) %>% summarise(pre_ep = any(aut_ep == 1), .groups = "drop")
    treated <- onsets %>% filter(g == !!g) %>%
      left_join(ep_pre, by = "country_code") %>%
      filter(is.na(pre_ep) | !pre_ep) %>% pull(country_code)
    controls <- ep_win %>% filter(!any_ep) %>% pull(country_code)
    panel %>%
      filter(year %in% win, country_code %in% c(treated, controls)) %>%
      mutate(stack = g, treated = as.integer(country_code %in% treated),
             event_time = year - g)
  }))
}

# fit_event_study(): build the stack, estimate the event study, and run the
# pre-trend test from the same model. Returns a list with
#   $est      coefficients by event time (est, se, 95% CI)
#   $pretest  Wald test that all pre-onset coefficients are jointly 0
#   $raw      raw averages by event time, treated vs controls
#   $n_treated, $n_controls, $label, $cohorts, $years
fit_event_study <- function(panel, cohorts, pre, post, label) {
  d   <- build_stack(panel, cohorts, pre, post) %>% filter(!is.na(y))
  ets <- setdiff(-pre:post, REF)
  for (k in ets) {
    d[[paste0("et_", ifelse(k < 0, paste0("m", abs(k)), k))]] <-
      as.integer(d$treated == 1 & d$event_time == k)
  }
  et_vars  <- grep("^et_", names(d), value = TRUE)
  d$unit   <- interaction(d$stack, d$country_code, drop = TRUE)
  d$period <- interaction(d$stack, d$year, drop = TRUE)
  m  <- lm(as.formula(paste("y ~", paste(et_vars, collapse = " + "),
                            "+ factor(unit) + factor(period)")), data = d)
  vc <- vcovCL(m, cluster = ~country_code)
  ct <- coeftest(m, vcov. = vc)[et_vars, , drop = FALSE]

  est <- data.frame(term = et_vars, est = ct[, 1], se = ct[, 2], row.names = NULL) %>%
    mutate(event_time = as.integer(sub("m", "-", sub("et_", "", term)))) %>%
    bind_rows(data.frame(term = "ref", est = 0, se = 0, event_time = REF)) %>%
    arrange(event_time) %>%
    mutate(lo = est - 1.96 * se, hi = est + 1.96 * se, outcome = label)

  leads <- grep("^et_m", et_vars, value = TRUE)
  b <- coef(m)[leads]; V <- vc[leads, leads, drop = FALSE]
  W <- as.numeric(t(b) %*% solve(V) %*% b)

  n_tr <- n_distinct(d$country_code[d$treated == 1])
  n_co <- n_distinct(d$country_code[d$treated == 0])
  pretest <- data.frame(outcome = label, pre_years = pre, n_treated = n_tr,
                        leads_tested = paste(sub("et_m", "-", leads), collapse = ","),
                        wald_chi2 = round(W, 2), df = length(leads),
                        p_value = round(pchisq(W, length(leads), lower.tail = FALSE), 3))

  raw <- d %>%
    group_by(event_time,
             group = ifelse(treated == 1, names(GROUP_COLS)[1], names(GROUP_COLS)[2])) %>%
    summarise(mean_y = mean(y), n = n(), .groups = "drop") %>%
    mutate(outcome = label)

  list(est = est, pretest = pretest, raw = raw, n_treated = n_tr, n_controls = n_co,
       label = label, cohorts = range(cohorts), years = range(panel$year),
       pre = pre, post = post)
}

# label_n(): "Outcome (22 treated, 121 controls; pre-trend p = 0.13)"
label_n <- function(fit) {
  sprintf("%s\n(%d treated, %d controls; pre-trend p = %.2f)", fit$label,
          fit$n_treated, fit$n_controls, fit$pretest$p_value)
}

# design_note(): one-line description of the sample, used in subtitles
design_note <- function(fit) {
  sprintf("Data %d-%d | onset cohorts %d-%d | window %d to +%d years | reference year -1",
          fit$years[1], fit$years[2], fit$cohorts[1], fit$cohorts[2], -fit$pre, fit$post)
}

# Building blocks for the two panel types ----------------------------
panel_es <- function(est, ylab, title = NULL) {
  ggplot(est, aes(event_time, est)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
    geom_vline(xintercept = -0.5, linetype = "dashed", colour = "grey50") +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey60", alpha = 0.25) +
    geom_line(colour = "grey15", linewidth = 0.7) +
    geom_point(colour = "grey15", size = 2) +
    scale_x_continuous(breaks = seq(-10, 10, 1)) +
    labs(x = "Years since autocratization onset", y = ylab, title = title) +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
}

panel_raw <- function(raw, ylab, title = NULL) {
  ggplot(raw, aes(event_time, mean_y, colour = group, shape = group)) +
    geom_vline(xintercept = -0.5, linetype = "dashed", colour = "grey50") +
    geom_line(linewidth = 0.7) + geom_point(size = 2) +
    scale_colour_manual(values = GROUP_COLS) +
    scale_shape_manual(values = c(16, 17)) +
    scale_x_continuous(breaks = seq(-10, 10, 1)) +
    labs(x = "Years since autocratization onset", y = ylab, colour = NULL, shape = NULL,
         title = title) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"))
}

# plot_pair(): raw trends (left) + event study (right) for ONE outcome
plot_pair <- function(fit, ylab, title, file, question) {
  p_raw <- panel_raw(fit$raw, paste("Average", ylab), "A. Raw trends (parallel-trends check)")
  p_es  <- panel_es(fit$est, paste("Gap vs controls,", ylab),
                    sprintf("B. Event study (%d treated, %d controls)",
                            fit$n_treated, fit$n_controls)) +
    labs(subtitle = sprintf("Pre-trend test p = %.2f (%s)", fit$pretest$p_value,
                            ifelse(fit$pretest$p_value > 0.1, "no sign of diverging pre-trends",
                                   "pre-trends differ: read with care")))
  p <- (p_raw | p_es) +
    plot_annotation(title = title,
                    subtitle = paste0("Question: ", question, "\n", design_note(fit)),
                    caption = paste(CAPTION_RAW, CAPTION_ES, sep = "\n"),
                    theme = theme(plot.title = element_text(face = "bold", size = 13),
                                  plot.caption = element_text(hjust = 0, size = 8, colour = "grey30")))
  save_fig(p, file, 12, 6)
  invisible(p)
}

# plot_facets(): event studies for SEVERAL outcomes, one small panel each
plot_facets <- function(fits, ylab, title, file, question, ncol = 3, w = 11, h = 7.5) {
  est <- bind_rows(lapply(fits, function(f) mutate(f$est, facet = label_n(f))))
  est$facet <- factor(est$facet, levels = unique(est$facet))
  p <- panel_es(est, ylab) + facet_wrap(~facet, ncol = ncol) +
    labs(title = title,
         subtitle = paste0("Question: ", question, "\n", design_note(fits[[1]])),
         caption = CAPTION_ES) +
    theme(plot.title = element_text(face = "bold", size = 13),
          plot.caption = element_text(hjust = 0, size = 8, colour = "grey30"),
          strip.text = element_text(face = "bold"))
  save_fig(p, file, w, h)
  invisible(p)
}

# plot_raw_facets(): raw trends for SEVERAL outcomes (parallel-trends check)
plot_raw_facets <- function(fits, ylab, title, file, question, ncol = 3, w = 11, h = 7.5) {
  raw <- bind_rows(lapply(fits, function(f) mutate(f$raw, facet = f$label)))
  raw$facet <- factor(raw$facet, levels = unique(raw$facet))
  p <- panel_raw(raw, ylab) + facet_wrap(~facet, ncol = ncol, scales = "free_y") +
    labs(title = title,
         subtitle = paste0("Question: ", question, "\n", design_note(fits[[1]])),
         caption = CAPTION_RAW) +
    theme(plot.title = element_text(face = "bold", size = 13),
          plot.caption = element_text(hjust = 0, size = 8, colour = "grey30"),
          strip.text = element_text(face = "bold"))
  save_fig(p, file, w, h)
  invisible(p)
}

# plot_overlay(): two or three event studies on the same axes
plot_overlay <- function(fits, ylab, title, file, question, w = 9, h = 6) {
  est <- bind_rows(lapply(fits, function(f) mutate(f$est, series = label_n(f))))
  est$series <- factor(est$series, levels = unique(est$series))
  cols <- c("#2a78d6", "#eb6834", "#1baf7a")[seq_along(fits)]
  p <- ggplot(est, aes(event_time, est, colour = series, fill = series, shape = series)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
    geom_vline(xintercept = -0.5, linetype = "dashed", colour = "grey50") +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.12, colour = NA) +
    geom_line(linewidth = 0.8) + geom_point(size = 2.2) +
    scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
    scale_x_continuous(breaks = seq(-10, 10, 1)) +
    labs(x = "Years since autocratization onset", y = ylab,
         colour = NULL, fill = NULL, shape = NULL, title = title,
         subtitle = paste0("Question: ", question, "\n", design_note(fits[[1]])),
         caption = CAPTION_ES) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom", legend.direction = "vertical",
          panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 13),
          plot.caption = element_text(hjust = 0, size = 8, colour = "grey30"))
  save_fig(p, file, w, h)
  invisible(p)
}

save_fig <- function(p, file, w, h) {
  ggsave(file.path(FIG_DIR, file), p, width = w, height = h, dpi = 150, bg = "white")
  message("  saved figure: outputs/figures/", file)
}

# summary_row(): one line per test for the overview table
summary_row <- function(fit, test_id) {
  e <- fit$est
  pick <- function(k) {
    r <- e[e$event_time == k, ]
    if (nrow(r) == 0) return(NA_character_)
    sprintf("%.2f [%.2f, %.2f]%s", r$est, r$lo, r$hi, ifelse(r$lo > 0 | r$hi < 0, " *", ""))
  }
  data.frame(test = test_id, outcome = fit$label,
             data_years = paste(fit$years, collapse = "-"),
             onset_cohorts = paste(fit$cohorts, collapse = "-"),
             n_treated = fit$n_treated, n_controls = fit$n_controls,
             pretrend_p = fit$pretest$p_value,
             t0 = pick(0), t1 = pick(1), t2 = pick(2), t3 = pick(3))
}

# print_fit(): compact console summary
print_fit <- function(fit) {
  cat(sprintf("  %s: %d treated, %d controls, pre-trend p = %.3f\n", fit$label,
              fit$n_treated, fit$n_controls, fit$pretest$p_value))
  print(fit$est %>% select(event_time, est, se, lo, hi) %>%
          mutate(across(-event_time, ~round(.x, 2))), row.names = FALSE)
}
