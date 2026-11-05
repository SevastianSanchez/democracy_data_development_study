# =====================================================================
# 02_run_tests.R  |  Option B preliminary tests: run every test
# ---------------------------------------------------------------------
# STATUS: TESTING ONLY. Exploratory checks to decide the paper's
# direction. Not the paper's final models.
#
# QUESTION
#   When a country starts a V-Dem autocratization episode, does its
#   statistical capacity change, compared with similar years in countries
#   that are not autocratizing? And WHICH part of capacity moves?
#
# RUN FROM: the repository root (the folder with stat_casualties_study.Rproj)
#   source("testing_functions_n_mods/optionB_testing/02_run_tests.R")
#   Takes about 5-10 minutes (many fixed-effect regressions).
#
# INPUTS  (all in data/, built by 00_prepare_inputs.R)
#   vdem_ert_episodes.csv, sci_dimensions.csv, spi_index.csv,
#   hrv_v1.csv, hrv2.csv
#
# OUTPUTS
#   outputs/figures/figNN_*.png   one figure per test (titled, captioned)
#   outputs/tables/*.csv          estimates, pre-trend tests, sample notes
#
# TESTS (each section below states its question, data, design, outputs)
#   Step 1  Overall indices: SCI (2004-2020), SPI (2016-2024)
#   Step 2  Which SPI pillar moves? All five pillars
#   Step 3  Headline contrast: Pillar 2 (dissemination) vs Pillar 5 (infrastructure)
#   Step 4  Parallel-trends checks for the pillars
#   Step 5  Robustness: drop the 2016 SPI wave
#   Step 6  What inside Pillar 2 moves? Sub-indicators 2.1, 2.2, 2.4
#   Step 7  Longer run with SCI dimensions (2004-2020)
#   Step 8  Pre-2016 dissemination: HRV transparency index, original version
#   Step 9  HRV2 (2024 update), full sample
#   Step 10 HRV2 split by period and by country set
#   Step 11 Why would post-2008 onsets look different? Erosion speed and length
#   Step 12 Sample accounting: who is in the SPI tests and who is not
#   Final   Overview table of all tests
# =====================================================================

source("testing_functions_n_mods/optionB_testing/01_helpers.R")

fits <- list()   # every fitted test, used for the overview table at the end

# ---------------------------------------------------------------------
# Load outcome data (all from data/)
# ---------------------------------------------------------------------
sci  <- read.csv(file.path(DATA_DIR, "sci_dimensions.csv"))
spi  <- read.csv(file.path(DATA_DIR, "spi_index.csv"))
hrv1 <- read.csv(file.path(DATA_DIR, "hrv_v1.csv"))
hrv2 <- read.csv(file.path(DATA_DIR, "hrv2.csv"))

# =====================================================================
# STEP 1. Overall capacity indices
# ---------------------------------------------------------------------
# Question : Do the overall indices move after an onset?
# Data     : SCI overall 2004-2020 (low/middle-income only; onset cohorts
#            2008-2017, 4 years before, 3 after). SPI overall 2016-2024
#            (all incomes; cohorts 2018-2021, 2 years before, 3 after).
# Outputs  : fig01_sci_overall.png, fig02_spi_overall.png
# =====================================================================
cat("\n=== STEP 1: overall indices ===\n")
fits$sci_overall <- fit_event_study(make_panel(sci, "sci_overall"), 2008:2017, 4, 3,
                                    "SCI overall")
plot_pair(fits$sci_overall, "SCI (0-100)",
          "Fig 1. Statistical Capacity Indicator (SCI) around autocratization onsets",
          "fig01_sci_overall.png",
          "Does overall statistical capacity fall after a country starts autocratizing? (SCI, low/middle-income countries)")
print_fit(fits$sci_overall)

fits$spi_overall <- fit_event_study(make_panel(spi, "spi_overall"), 2018:2021, 2, 3,
                                    "SPI overall")
plot_pair(fits$spi_overall, "SPI (0-100)",
          "Fig 2. Statistical Performance Indicators (SPI) around autocratization onsets",
          "fig02_spi_overall.png",
          "Does overall statistical performance fall after onset? (SPI, all income groups)")
print_fit(fits$spi_overall)

# =====================================================================
# STEP 2. Which part of capacity moves? All five SPI pillars
# ---------------------------------------------------------------------
# Question : Are some parts of statistical capacity more sensitive to
#            erosion than others?
# Data     : SPI pillar scores 2016-2024 (pillar scores start in 2016),
#            onset cohorts 2018-2021, 2 years before, 3 after. Same sample
#            design for every pillar so the panels are comparable.
#   Pillar 1 Data use        | Pillar 2 Data services (dissemination)
#   Pillar 3 Data products   | Pillar 4 Data sources | Pillar 5 Data infrastructure
# Outputs  : fig03_spi_pillars_event_study.png, fig04_spi_pillars_raw_trends.png
# =====================================================================
cat("\n=== STEP 2: all SPI pillars ===\n")
pillars <- c(p1_data_use      = "Pillar 1: Data use",
             p2_data_services = "Pillar 2: Data services (dissemination)",
             p3_data_products = "Pillar 3: Data products",
             p4_data_sources  = "Pillar 4: Data sources",
             p5_data_infra    = "Pillar 5: Data infrastructure")
pillar_fits <- lapply(names(pillars), function(v)
  fit_event_study(make_panel(spi, v, from = 2016), 2018:2021, 2, 3, pillars[[v]]))
names(pillar_fits) <- names(pillars)
fits <- c(fits, setNames(pillar_fits, paste0("spi_", names(pillars))))

plot_facets(pillar_fits, "Gap vs controls (pillar score, 0-100)",
            "Fig 3. All five SPI pillars around autocratization onsets (event study)",
            "fig03_spi_pillars_event_study.png",
            "Which parts of statistical capacity move after onset? Is any pillar more sensitive to erosion than the others?")
plot_raw_facets(pillar_fits, "Average pillar score (0-100)",
                "Fig 4. All five SPI pillars: raw trends (parallel-trends check)",
                "fig04_spi_pillars_raw_trends.png",
                "Before onset, do autocratizing countries and controls move together on each pillar?")
for (f in pillar_fits) print_fit(f)

# =====================================================================
# STEP 3. Headline contrast: Pillar 2 vs Pillar 5
# ---------------------------------------------------------------------
# Question : The mechanism says erosion hits DISSEMINATION (what the
#            statistical office releases and how openly), not the
#            machinery behind it. So Pillar 2 should fall, Pillar 5 not.
# Data     : same fits as Step 2, drawn on one set of axes.
# Output   : fig05_spi_p2_vs_p5.png
# =====================================================================
plot_overlay(list(pillar_fits$p2_data_services, pillar_fits$p5_data_infra),
             "Gap vs controls (pillar score, 0-100)",
             "Fig 5. Pillar 2 (dissemination) falls after onset, Pillar 5 (infrastructure) does not",
             "fig05_spi_p2_vs_p5.png",
             "Does erosion hit dissemination (Pillar 2) rather than infrastructure (Pillar 5)?")

# =====================================================================
# STEP 4. Parallel-trends checks for the pillars
# ---------------------------------------------------------------------
# Question : Were treated and control countries already on different paths
#            before onset?
# Checks   : (a) raw trends (Fig 4 above)
#            (b) Wald test that the pre-onset coefficients are jointly 0,
#                with 2 pre-years (cohorts 2018-2021) and 3 pre-years
#                (cohorts 2019-2021). Low power: only 1-2 leads.
#            (c) event study re-estimated with 3 pre-years
# Outputs  : tables/pretrend_tests.csv, fig06_spi_pillars_3yr_pre.png
# =====================================================================
cat("\n=== STEP 4: parallel-trends checks ===\n")
pillar_fits_pre3 <- lapply(names(pillars), function(v)
  fit_event_study(make_panel(spi, v, from = 2016), 2019:2021, 3, 3, pillars[[v]]))
plot_facets(pillar_fits_pre3, "Gap vs controls (pillar score, 0-100)",
            "Fig 6. SPI pillars with a longer pre-period (3 years before onset)",
            "fig06_spi_pillars_3yr_pre.png",
            "Do the pillar results hold, and do pre-trends stay flat, with one more year before onset?")

# =====================================================================
# STEP 5. Robustness: drop the 2016 SPI wave
# ---------------------------------------------------------------------
# Question : 2016 values look like a method break (Pillar 2 average:
#            47 in 2016, 56 in 2017, 68 in 2018).
#            Does the Pillar 2 vs Pillar 5 pattern survive without 2016?
# Data     : SPI 2017-2024, onset cohorts 2019-2021, 2 before, 3 after.
# Output   : fig07_spi_p2_p5_drop2016.png
# =====================================================================
cat("\n=== STEP 5: drop 2016 ===\n")
fits$spi_p2_no2016 <- fit_event_study(make_panel(spi, "p2_data_services", from = 2017),
                                      2019:2021, 2, 3, "Pillar 2: Data services, 2017+")
fits$spi_p5_no2016 <- fit_event_study(make_panel(spi, "p5_data_infra", from = 2017),
                                      2019:2021, 2, 3, "Pillar 5: Data infrastructure, 2017+")
plot_overlay(list(fits$spi_p2_no2016, fits$spi_p5_no2016),
             "Gap vs controls (pillar score, 0-100)",
             "Fig 7. Pillar 2 vs Pillar 5 without the 2016 SPI wave",
             "fig07_spi_p2_p5_drop2016.png",
             "Is the Pillar 2 drop an artifact of the 2016 method change?")
print_fit(fits$spi_p2_no2016); print_fit(fits$spi_p5_no2016)

# =====================================================================
# STEP 6. What inside Pillar 2 moves?
# ---------------------------------------------------------------------
# Question : Which dissemination channel drives the Pillar 2 drop?
# Data     : SPI 2016-2024, cohorts 2018-2021, 2 before, 3 after.
#   2.1 Data releases (IMF data standards SDDS / e-GDDS subscription, 0-1)
#   2.2 Online access (ODIN openness score, 0-1)
#   2.4 Data access services (microdata portal, NADA, 0-1)
#   (2.3 advisory services has no data in the SPI file)
#   Note: 2.1 barely varies (few countries change IMF status), so R may
#   print a harmless 'NaNs produced' warning from the clustered SEs.
# Output   : fig08_spi_p2_subindicators.png
# =====================================================================
cat("\n=== STEP 6: Pillar 2 sub-indicators ===\n")
subs <- c(p2_1_data_releases   = "2.1 Data releases (IMF standards)",
          p2_2_online_access   = "2.2 Online access (ODIN openness)",
          p2_4_data_access_svc = "2.4 Data access services (microdata portal)")
sub_fits <- lapply(names(subs), function(v)
  fit_event_study(make_panel(spi, v, from = 2016), 2018:2021, 2, 3, subs[[v]]))
fits <- c(fits, setNames(sub_fits, paste0("spi_", names(subs))))
plot_facets(sub_fits, "Gap vs controls (sub-indicator score, 0-1)",
            "Fig 8. Inside Pillar 2: which dissemination channel moves?",
            "fig08_spi_p2_subindicators.png",
            "Is the Pillar 2 drop driven by data releases, online openness, or microdata access?",
            ncol = 3, w = 11, h = 5)
for (f in sub_fits) print_fit(f)

# =====================================================================
# STEP 7. Longer run with SCI dimensions (2004-2020)
# ---------------------------------------------------------------------
# Question : Does the dissemination-vs-production contrast show up in
#            the longer SCI series?
# Data     : SCI 2004-2020 (low/middle-income only), cohorts 2008-2017,
#            4 before, 3 after.
#   Periodicity & timeliness (closest to dissemination, but measures how
#   often surveys/censuses happen, not how openly data are released)
#   Methodology (standards, classifications) | Source data (censuses, surveys)
# Outputs  : fig09_sci_dimensions_event_study.png, fig10_sci_dimensions_raw_trends.png
# =====================================================================
cat("\n=== STEP 7: SCI dimensions ===\n")
sci_dims <- c(sci_periodicity = "SCI: Periodicity & timeliness",
              sci_methodology = "SCI: Methodology",
              sci_source_data = "SCI: Source data",
              sci_overall     = "SCI: Overall")
sci_fits <- lapply(names(sci_dims), function(v)
  fit_event_study(make_panel(sci, v), 2008:2017, 4, 3, sci_dims[[v]]))
fits <- c(fits, setNames(sci_fits[1:3], paste0("sci_", c("periodicity", "methodology", "source"))))
plot_facets(sci_fits, "Gap vs controls (SCI score, 0-100)",
            "Fig 9. SCI dimensions around autocratization onsets (event study)",
            "fig09_sci_dimensions_event_study.png",
            "Does any SCI dimension fall after onset? (SCI, low/middle-income countries)",
            ncol = 2, w = 10, h = 7.5)
plot_raw_facets(sci_fits, "Average SCI score (0-100)",
                "Fig 10. SCI dimensions: raw trends (parallel-trends check)",
                "fig10_sci_dimensions_raw_trends.png",
                "Before onset, do autocratizing countries and controls move together on each SCI dimension?",
                ncol = 2, w = 10, h = 7.5)
for (f in sci_fits) print_fit(f)

# =====================================================================
# STEP 8. Pre-2016 dissemination: HRV transparency index (original)
# ---------------------------------------------------------------------
# Question : Does the dissemination drop also appear before 2016, in a
#            measure that covers all income groups?
# Data     : HRV v1 (Hollyer, Rosendorff & Vreeland), 1990-2010,
#            125 countries incl. high-income. The index measures whether
#            governments REPORT data to the World Bank (WDI), i.e. data
#            disclosure. Onset cohorts 1994-2007, 4 before, 3 after.
# Output   : fig11_hrv_v1.png
# =====================================================================
cat("\n=== STEP 8: HRV v1 ===\n")
hrv1_pan <- make_panel(hrv1, "hrv_index", from = 1990)
fits$hrv_v1 <- fit_event_study(hrv1_pan, 1994:2007, 4, 3, "HRV transparency index (v1)")
plot_pair(fits$hrv_v1, "HRV index",
          "Fig 11. HRV transparency index (original) around autocratization onsets, 1994-2007",
          "fig11_hrv_v1.png",
          "Does data disclosure fall after onset in the pre-2016 period, across all income groups?")
print_fit(fits$hrv_v1)
cat("  Typical within-country yearly change (SD):",
    round(mean(tapply(hrv1_pan$y, hrv1_pan$country_code, function(x) sd(diff(x))), na.rm = TRUE), 3), "\n")

# =====================================================================
# STEP 9. HRV2 (2024 update), full sample
# ---------------------------------------------------------------------
# Question : Does the HRV result hold in the updated, larger index?
# Data     : HRV2 (Hollyer, Ge, Rosendorff & Vreeland 2024 draft),
#            1990-2015, 149 countries (24 added countries start in 1993).
#            Onset cohorts 1994-2012, 4 before, 3 after.
# Output   : fig12_hrv2_full.png
# =====================================================================
cat("\n=== STEP 9: HRV2 full sample ===\n")
hrv2_pan <- make_panel(hrv2, "hrv_index", from = 1990)
fits$hrv2_full <- fit_event_study(hrv2_pan, 1994:2012, 4, 3, "HRV2 transparency index, full sample")
plot_pair(fits$hrv2_full, "HRV2 index",
          "Fig 12. HRV2 transparency index (2024 update) around autocratization onsets, 1994-2012",
          "fig12_hrv2_full.png",
          "Does the disclosure drop hold in the updated index with more countries and years?")
print_fit(fits$hrv2_full)
cat("  Typical within-country yearly change (SD):",
    round(mean(tapply(hrv2_pan$y, hrv2_pan$country_code, function(x) sd(diff(x))), na.rm = TRUE), 3), "\n")

# =====================================================================
# STEP 10. HRV2 split by period and by country set
# ---------------------------------------------------------------------
# Question : Why is the HRV2 full-sample drop smaller? Is it the newer
#            onsets (2008-2012) or the 24 added countries?
# Splits   : (a) HRV2, onset cohorts 1994-2007
#            (b) HRV2, onset cohorts 2008-2012
#            (c) HRV2, only the 125 original HRV v1 countries, cohorts 1994-2007
# Outputs  : fig13_hrv2_splits.png; HRV v1 vs HRV2 correlation printed
# =====================================================================
cat("\n=== STEP 10: HRV2 splits ===\n")
orig_countries <- unique(hrv1$country_code)
fits$hrv2_9407 <- fit_event_study(hrv2_pan, 1994:2007, 4, 3, "HRV2, onsets 1994-2007")
fits$hrv2_0812 <- fit_event_study(hrv2_pan, 2008:2012, 4, 3, "HRV2, onsets 2008-2012")
fits$hrv2_orig <- fit_event_study(hrv2_pan %>% filter(country_code %in% orig_countries),
                                  1994:2007, 4, 3, "HRV2, original 125 countries, onsets 1994-2007")
plot_facets(list(fits$hrv2_9407, fits$hrv2_0812, fits$hrv2_orig),
            "Gap vs controls (HRV2 index)",
            "Fig 13. HRV2 split: the drop is concentrated in earlier onsets and the original countries",
            "fig13_hrv2_splits.png",
            "Is the smaller HRV2 effect driven by post-2008 onsets or by the 24 added countries?",
            ncol = 3, w = 12, h = 5)
for (f in list(fits$hrv2_9407, fits$hrv2_0812, fits$hrv2_orig)) print_fit(f)
hrv_cor <- merge(hrv1_pan, hrv2_pan, by = c("country_code", "year"))
cat("  Correlation HRV v1 vs HRV2 (same country-years):", round(cor(hrv_cor$y.x, hrv_cor$y.y), 2), "\n")

# =====================================================================
# STEP 11. Why would post-2008 onsets look different? Erosion speed
# ---------------------------------------------------------------------
# Question : Is democratic decline after 2008 more gradual, so a 3-year
#            event window catches less of it?
# Data     : V-Dem electoral democracy index (v2x_polyarchy) around every
#            autocratization onset 1994-2012, split at 2008.
# Measures : change in polyarchy from the year before onset (-1) to
#            years +1 and +3; episode length (V-Dem episode start to end;
#            episodes still running in 2024 are cut off, so their length
#            is a lower bound); total decline over the episode; average
#            yearly decline (total / length).
# Outputs  : fig14_erosion_speed_by_period.png,
#            tables/erosion_speed_by_period.csv,
#            tables/erosion_speed_by_onset.csv (one row per onset)
# =====================================================================
cat("\n=== STEP 11: erosion speed by period ===\n")
poly_at <- function(cc, yr) {
  v <- ert$v2x_polyarchy[ert$country_code == cc & ert$year == yr]
  if (length(v) == 0) NA_real_ else v
}
ep_dates <- vdem %>% filter(!is.na(aut_ep_start_year)) %>%
  distinct(country_code, g = aut_ep_start_year, ep_end = aut_ep_end_year)
speed <- onsets %>% filter(g >= 1994, g <= 2012) %>%
  left_join(ep_dates, by = c("country_code", "g")) %>%
  mutate(period = ifelse(g <= 2007, "Onsets 1994-2007", "Onsets 2008-2012")) %>%
  rowwise() %>%
  mutate(p_before = poly_at(country_code, g - 1),
         change_y1 = poly_at(country_code, g + 1) - p_before,
         change_y3 = poly_at(country_code, g + 3) - p_before,
         episode_years = ep_end - g + 1,
         total_decline = poly_at(country_code, ep_end) - p_before,
         yearly_decline = total_decline / episode_years) %>%
  ungroup()
write.csv(speed %>% mutate(across(where(is.double), ~round(.x, 3))),
          file.path(TAB_DIR, "erosion_speed_by_onset.csv"), row.names = FALSE)

speed_tab <- speed %>% group_by(period) %>%
  summarise(n_onsets = n(),
            mean_change_year1 = round(mean(change_y1, na.rm = TRUE), 3),
            mean_change_year3 = round(mean(change_y3, na.rm = TRUE), 3),
            median_episode_years = median(episode_years),
            share_short_episodes_le2yrs = round(mean(episode_years <= 2), 2),
            share_long_episodes_ge10yrs = round(mean(episode_years >= 10), 2),
            median_yearly_decline = round(median(yearly_decline), 3),
            .groups = "drop")
print(as.data.frame(speed_tab), row.names = FALSE)
write.csv(speed_tab, file.path(TAB_DIR, "erosion_speed_by_period.csv"), row.names = FALSE)

# Panel A: average path of polyarchy around onset, by period
n_by_period <- speed %>% count(period) %>% mutate(lab = sprintf("%s (n = %d)", period, n))
speed_path <- speed %>%
  select(country_code, g, period, p_before) %>%
  crossing(event_time = -1:5) %>%
  rowwise() %>% mutate(change = poly_at(country_code, g + event_time) - p_before) %>% ungroup()
speed_avg <- speed_path %>% group_by(period, event_time) %>%
  summarise(mean_change = mean(change, na.rm = TRUE),
            se = sd(change, na.rm = TRUE) / sqrt(sum(!is.na(change))), .groups = "drop") %>%
  mutate(lo = mean_change - 1.96 * se, hi = mean_change + 1.96 * se) %>%
  left_join(n_by_period, by = "period")
period_cols <- setNames(c(COL_TREATED, COL_CONTROL), n_by_period$lab)

pA <- ggplot(speed_avg, aes(event_time, mean_change, colour = lab, fill = lab, shape = lab)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = -0.5, linetype = "dashed", colour = "grey50") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.8) + geom_point(size = 2.2) +
  scale_colour_manual(values = period_cols) + scale_fill_manual(values = period_cols) +
  scale_x_continuous(breaks = -1:5) +
  labs(x = "Years since autocratization onset",
       y = "Change in V-Dem electoral democracy index\n(vs year before onset)",
       colour = NULL, fill = NULL, shape = NULL,
       title = "A. Average path after onset (band = 95% CI)") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

# Panel B: each onset, episode length vs total decline
speed_pts <- speed %>% left_join(n_by_period, by = "period")
pB <- ggplot(speed_pts, aes(episode_years, total_decline, colour = lab, shape = lab)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_point(size = 2.4, alpha = 0.85) +
  geom_text(data = speed_pts %>% filter(period == "Onsets 2008-2012"),
            aes(label = country_code), size = 2.6, vjust = -0.9, show.legend = FALSE) +
  scale_colour_manual(values = period_cols) +
  labs(x = "Episode length in years (episodes still running in 2024 are cut off)",
       y = "Total decline in electoral democracy index\nover the episode",
       colour = NULL, shape = NULL,
       title = "B. Each onset: how long the episode lasted vs how far democracy fell") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

p_speed <- (pA | pB) +
  plot_annotation(
    title = "Fig 14. Erosion speed: post-2008 episodes last longer and fall more slowly per year (thin evidence: 15 onsets)",
    subtitle = paste0("Question: Is post-2008 backsliding more gradual, so a 3-year window catches less of it?\n",
                      "All V-Dem autocratization onsets 1994-2012 (not only the clean onsets used in the event studies)"),
    caption = paste(
      "How to read: in A, a lower line means democracy had fallen further since the year before onset.",
      "In B, each dot is one onset; points far to the right are long episodes, points low down are large declines.",
      "Post-2008 cases (labelled) include short drops of 1-2 years that partly reverse (as before 2008) and more episodes lasting 10+ years with slow yearly decline (33% vs 19%).",
      "Both types are hard to see in a 3-year window, which motivates measuring the SIZE of the democratic decline (dose-response).",
      sep = "\n"),
    theme = theme(plot.title = element_text(face = "bold", size = 13),
                  plot.caption = element_text(hjust = 0, size = 8, colour = "grey30")))
save_fig(p_speed, "fig14_erosion_speed_by_period.png", 13, 6.5)

# =====================================================================
# STEP 12. Sample accounting for the SPI tests
# ---------------------------------------------------------------------
# Question : SPI covers ~190 countries. Why do the pillar tests use ~140?
# Data     : SPI Pillar 2, 2016-2024, cohorts 2018-2021 (as in Step 2).
# Output   : tables/sample_accounting_spi.csv
# =====================================================================
cat("\n=== STEP 12: sample accounting (SPI Pillar 2) ===\n")
p2_pan  <- make_panel(spi, "p2_data_services", from = 2016)
p2_stk  <- build_stack(p2_pan, 2018:2021, 2, 3)
spi_c   <- unique(p2_pan$country_code)
vdem_c  <- unique(in_ep$country_code)
used    <- unique(p2_stk$country_code)
notused <- setdiff(spi_c[spi_c %in% vdem_c], used)
ep_nu   <- in_ep %>% filter(year >= 2016, country_code %in% notused) %>%
  group_by(country_code) %>% summarise(any_ep = any(aut_ep == 1), .groups = "drop")
late_or_early <- onsets %>% filter(country_code %in% notused, g >= 2016) %>% pull(country_code) %>% unique()
acct <- data.frame(
  group = c("Countries with SPI Pillar 2 scores",
            "  not in V-Dem (mostly small states), so no episode data",
            "  in V-Dem",
            "    used in the test (treated or clean control)",
            "      of which treated (clean onset 2018-2021)",
            "    not used: onset in 2016-2017 or 2022+ (outside usable cohorts)",
            "    not used: in an episode that began before 2016 (ongoing backsliders)",
            "    not used: other (not clean in any window)"),
  n = c(length(spi_c), sum(!spi_c %in% vdem_c), sum(spi_c %in% vdem_c), length(used),
        n_distinct(p2_stk$country_code[p2_stk$treated == 1]),
        length(late_or_early),
        sum(ep_nu$any_ep & !ep_nu$country_code %in% late_or_early),
        length(notused) - length(late_or_early) -
          sum(ep_nu$any_ep & !ep_nu$country_code %in% late_or_early)))
print(acct, row.names = FALSE)
cat("  Not in V-Dem:", paste(sort(setdiff(spi_c, vdem_c)), collapse = ", "), "\n")
write.csv(acct, file.path(TAB_DIR, "sample_accounting_spi.csv"), row.names = FALSE)

# =====================================================================
# FINAL. Overview tables
# ---------------------------------------------------------------------
# tables/summary_all_tests.csv    one row per test: sample, pre-trend p,
#                                 estimates at t0-t3 with 95% CI (* = CI excludes 0)
# tables/event_study_estimates.csv every coefficient from every test
# tables/pretrend_tests.csv        Wald pre-trend tests (incl. 3-pre-year SPI version)
# =====================================================================
summary_tab <- bind_rows(lapply(names(fits), function(n) summary_row(fits[[n]], n)))
write.csv(summary_tab, file.path(TAB_DIR, "summary_all_tests.csv"), row.names = FALSE)
write.csv(bind_rows(lapply(names(fits), function(n) mutate(fits[[n]]$est, test = n))),
          file.path(TAB_DIR, "event_study_estimates.csv"), row.names = FALSE)
pt_tab <- bind_rows(bind_rows(lapply(fits, `[[`, "pretest")),
                    bind_rows(lapply(pillar_fits_pre3, `[[`, "pretest")))
write.csv(pt_tab, file.path(TAB_DIR, "pretrend_tests.csv"), row.names = FALSE)

cat("\n=== OVERVIEW ===\n")
print(summary_tab, row.names = FALSE)
cat("\nPre-trend tests:\n"); print(pt_tab, row.names = FALSE)
message("Done. Figures in ", FIG_DIR, "; tables in ", TAB_DIR)
