# Sanity check: does the panel built by Comp2_panel_wrangling.R reproduce the
# sample sizes reported in the submitted paper? Run from the repo root.
# Stops with an error if any count differs.
#
# Paper (submitted to Governance): 155 countries, 2016-2023.
#   Table 1 (Stage I, SPI ~ DI) observations: FD 1,079 | FD lag 772 | TWFE 1,234 | TWFE lag 927
#   Table 2 (Stage II direct, SDG ~ DI):      FD 1,082 | FD lag 772 | TWFE 1,237 | TWFE lag 927
# Model formulas copied from panel_regression/Stage1_article.Rmd and Stage2_article.Rmd.
# Last confirmed: 2026-11-18 (all eight counts and the Table 1 coefficients match).

source("Comp2_panel_wrangling.R")
library(plm)

d  <- panel_data %>%
  dplyr::select(country_code, year, spi_comp, di_score, log_gdppc, sdg_overall) %>%
  dplyr::arrange(country_code, year)
pd <- pdata.frame(d, index = c("country_code", "year"))

f_cur <- function(y) as.formula(paste(y, "~ di_score + log_gdppc + factor(year)"))
f_dyn <- function(y) as.formula(paste(y, "~ di_score + plm::lag(di_score, 1) + plm::lag(di_score, 2) + log_gdppc + factor(year)"))

n_obs <- function(y) c(
  fd_cur = nobs(plm(f_cur(y), data = pd, model = "fd")),
  fd_dyn = nobs(plm(f_dyn(y), data = pd, model = "fd")),
  fe_cur = nobs(plm(f_cur(y), data = pd, model = "within")),
  fe_dyn = nobs(plm(f_dyn(y), data = pd, model = "within")))

got_t1 <- n_obs("spi_comp");    want_t1 <- c(1079, 772, 1234, 927)
got_t2 <- n_obs("sdg_overall"); want_t2 <- c(1082, 772, 1237, 927)
n_countries <- dplyr::n_distinct(d$country_code[!is.na(d$spi_comp) & !is.na(d$di_score) & !is.na(d$log_gdppc)])

print(rbind(table1_got = got_t1, table1_paper = want_t1, table2_got = got_t2, table2_paper = want_t2))
cat("countries in Table 1 sample:", n_countries, "(paper: 155)\n")

stopifnot(all(got_t1 == want_t1), all(got_t2 == want_t2), n_countries == 155)
cat("OK: panel reproduces the paper's sample sizes.\n")
