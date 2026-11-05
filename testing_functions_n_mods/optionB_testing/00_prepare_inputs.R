# =====================================================================
# 00_prepare_inputs.R  |  Option B preliminary tests: build input files
# ---------------------------------------------------------------------
# PURPOSE
#   Builds every input the Option B tests use and saves it in
#   testing_functions_n_mods/optionB_testing/data/, so the test folder
#   is self-contained. You only need to rerun this to refresh the data;
#   the finished files are already committed in data/.
#
# RUN FROM: the repository root (the folder with stat_casualties_study.Rproj)
#
# OUTPUT FILES (data/)
#   vdem_ert_episodes.csv   V-Dem Episodes of Regime Transformation (ERT):
#                           autocratization episode flag, episode dates,
#                           electoral democracy index (v2x_polyarchy)
#   sci_dimensions.csv      World Bank Statistical Capacity Indicator (SCI),
#                           overall score and its 3 dimensions, 2004-2020
#   spi_index.csv           World Bank Statistical Performance Indicators
#                           (SPI): overall, 5 pillars, Pillar 2 sub-indicators,
#                           2004-2024 (pillar scores start 2016)
#   hrv_v1.csv              HRV transparency index, original version
#                           (1980-2010, 125 countries), from QoG
#   hrv2.csv                HRV2 transparency index, 2024 update
#                           (1980-2015, 149 countries), from the authors' file
#
# SOURCES THAT ARE NOT IN THE REPO (set the paths below if refreshing)
#   SPI    : https://github.com/worldbank/SPI (03_output_data/SPI_index.csv)
#   HRV v1 : QoG Standard Time-Series dataset (Jan 2026), variable hrv_index.
#            Same values as the QoG HRV dataset page
#            (https://datafinder.qog.gu.se/dataset/hrv).
#   HRV2   : HRV2ExpandedCntriesExtracted.RData (object TransparencyData),
#            shared by Sevastian on 5 Nov 2026 (project files).
#            Hollyer, Ge, Rosendorff & Vreeland, "Measuring Transparency,
#            Again" (2024 draft).
# =====================================================================

library(dplyr)

data_dir <- "testing_functions_n_mods/optionB_testing/data"
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

# Paths to sources outside the repo (edit if refreshing) -------------
SPI_URL   <- "https://raw.githubusercontent.com/worldbank/SPI/master/03_output_data/SPI_index.csv"
QOG_FILE  <- "path/to/qog_std_ts_jan26.csv"                 # QoG standard TS (or HRV download)
HRV2_FILE <- "path/to/HRV2ExpandedCntriesExtracted.RData"   # authors' HRV2 file

# 1. V-Dem ERT episodes (from the repo) ------------------------------
read.csv("data/input_data/ert.csv") %>%
  select(country_code = country_text_id, country_name, year,
         aut_ep, aut_ep_start_year, aut_ep_end_year,
         dem_ep, v2x_polyarchy) %>%
  write.csv(file.path(data_dir, "vdem_ert_episodes.csv"), row.names = FALSE)

# 2. SCI overall + dimensions (from the repo) ------------------------
#    IQ.SCI.OVRL overall | IQ.SCI.MTHD methodology |
#    IQ.SCI.SRCE source data | IQ.SCI.PRDC periodicity and timeliness
read.csv("data/input_data/SCI_All_Dim_TS.csv", check.names = FALSE) %>%
  transmute(country_code, country_name, year = Year,
            sci_overall     = suppressWarnings(as.numeric(`IQ.SCI.OVRL`)),
            sci_methodology = suppressWarnings(as.numeric(`IQ.SCI.MTHD`)),
            sci_source_data = suppressWarnings(as.numeric(`IQ.SCI.SRCE`)),
            sci_periodicity = suppressWarnings(as.numeric(`IQ.SCI.PRDC`))) %>%
  write.csv(file.path(data_dir, "sci_dimensions.csv"), row.names = FALSE)

# 3. SPI (World Bank GitHub) -----------------------------------------
spi_try <- try(read.csv(SPI_URL, fileEncoding = "UTF-8-BOM"), silent = TRUE)
if (!inherits(spi_try, "try-error")) {
  spi_try %>%
    transmute(country_code = iso3c, country_name = country, year = date,
              spi_overall      = SPI.INDEX,
              p1_data_use      = SPI.INDEX.PIL1,
              p2_data_services = SPI.INDEX.PIL2,
              p3_data_products = SPI.INDEX.PIL3,
              p4_data_sources  = SPI.INDEX.PIL4,
              p5_data_infra    = SPI.INDEX.PIL5,
              p2_1_data_releases   = SPI.DIM2.1.INDEX,   # IMF data standards (SDDS/e-GDDS)
              p2_2_online_access   = SPI.DIM2.2.INDEX,   # ODIN openness
              p2_4_data_access_svc = SPI.DIM2.4.INDEX) %>% # microdata portal (NADA)
    write.csv(file.path(data_dir, "spi_index.csv"), row.names = FALSE)
} else message("SPI download failed; keeping existing data/spi_index.csv")

# 4. HRV v1 (QoG) ----------------------------------------------------
if (file.exists(QOG_FILE)) {
  qog <- read.csv(QOG_FILE)
  if (!"hrv_index" %in% names(qog)) qog <- read.csv(QOG_FILE, sep = ";", dec = ",")
  qog %>% filter(!is.na(hrv_index)) %>%
    transmute(country_code = ccodealp, country_name = cname, year,
              hrv_index, hrv_lb, hrv_ub, hrv_sd) %>%
    write.csv(file.path(data_dir, "hrv_v1.csv"), row.names = FALSE)
} else message("QOG_FILE not set; keeping existing data/hrv_v1.csv")

# 5. HRV2 (authors' RData) -------------------------------------------
if (file.exists(HRV2_FILE)) {
  e <- new.env(); suppressWarnings(load(HRV2_FILE, envir = e))
  with(e$TransparencyData,
       data.frame(country_code = ccode, country_name = Country, year,
                  hrv_index = transparencyindex, hrv_lb = transparencyindexlb,
                  hrv_ub = transparencyindexub, hrv_sd = transparencyindexsd)) %>%
    write.csv(file.path(data_dir, "hrv2.csv"), row.names = FALSE)
} else message("HRV2_FILE not set; keeping existing data/hrv2.csv")

message("Inputs written to ", data_dir)
