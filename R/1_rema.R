# GOA SST biomass estimation using the bottom trawl and longline survey indices

# Stock assessment in 2022 (projected to 2023) based on Sept 2022 GPT
# recommendations (endorsed by SSC at Oct 2022 SSC mtg):
# https://meetings.npfmc.org/CommentReview/DownloadFile?p=36b50cbe-9340-4170-bba0-347da1b0054d.pdf&fileName=C5%20GOA%20Groundfish%20Plan%20Team%20Minutes.pdf
# (1) The Team recommended excluding BTS data from 1984 and 1987 due to different survey
# methodology and to continue utilizing a two-survey model.
# (2) The Team recommended simplifying the model naming convention where Model 18 represents the
# status quo model, Model 18* is the corrected model in TMB with new data, and Model 22 is the
# model with additional observation error on BTS and LLS.
# (3) The Team recommended discontinuing the misspecified status quo model (Model 18) and bringing
# forward both the corrected model (Model 18*) and the mod

# Model naming conventions:
# Model 18* = m18: 1984-pres. corrected version of status quo model
# Model 18* no 1984/97 = m18s: 1990-pres. corrected version of status quo model
# Model 22 = m22: 1984-pres. xtra observation error for BTS and LLS
# Model 22 no 1984/87 = m22s: 1990-pres. xtra observation error for BTS and LLS

# Set up ----

# assessment year
YEAR <- 2022

libs <- c('readr', 'dplyr', 'tidyr', 'ggplot2', 'cowplot')
if(length(libs[which(libs %in% rownames(installed.packages()) == FALSE )]) > 0) {install.packages(libs[which(libs %in% rownames(installed.packages()) == FALSE)])}
lapply(libs, library, character.only = TRUE)

# install.packages("devtools")
# devtools::install_github("JaneSullivan-NOAA/rema", dependencies = TRUE)
library(rema)

# folder set up
dat_path <- paste0("data/", YEAR); dir.create(dat_path)
out_path <- paste0("results/", YEAR); dir.create(out_path)

ggplot2::theme_set(cowplot::theme_cowplot(font_size = 12) +
                     cowplot::background_grid() +
                     cowplot::panel_border())

# Read data ----

# bottom trawl survey
biomass_dat <- read_csv(paste0(dat_path, "/goa_sst_biomass_", YEAR, ".csv"))
biomass_dat %>%
  tidyr::expand(year = min(biomass_dat$year):(YEAR-1),
            strata) %>%
  left_join(biomass_dat %>%
                mutate(value = ifelse(is.na(biomass), NA,
                                      paste0(prettyNum(round(biomass, 0), big.mark = ','), ' (',
                                             format(round(cv, 3), nsmall = 3, trim = TRUE), ')')))) %>%
  pivot_wider(id_cols = c(year), names_from = strata, values_from = value) %>%
  arrange(year) %>%
  write_csv(paste0(out_path, '/biomass_data_wide.csv'))

# longline survey rpws
cpue_dat <- read_csv(paste0(dat_path, "/goa_sst_rpw_", YEAR, ".csv"))

# Model 18* ----
input <- prepare_rema_input(model_name = 'Model 18*',
                            multi_survey = 1,
                            biomass_dat = biomass_dat,
                            cpue_dat = cpue_dat,
                            # RPWs are summable
                            sum_cpue_index = TRUE,
                            end_year = YEAR + 1,
                            # three process error parameters (log_PE) estimated,
                            # indexed as follows for each biomass survey stratum
                            # (shared within an area across depths):
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1, 2, 2, 2, 3, 3, 3)),
                            q_options = list(
                              # LLS strata (n=3) indexed as follows for the
                              # biomass strata (n=9)
                              pointer_biomass_cpue_strata = c(1, 1, 1, 2, 2, 2, 3, 3, 3),
                              # one scaling parameters (log_q) estimated, shared
                              # over all three LLS strata
                              pointer_q_cpue = c(1, 1, 1)))
m18 <- fit_rema(input)
out18 <- tidy_rema(m18)
out18$parameter_estimates

# Model 18* no 1984/87 ----
input <- prepare_rema_input(model_name = 'Model 18* no 1984/87',
                            multi_survey = 1,
                            biomass_dat = biomass_dat,
                            cpue_dat = cpue_dat,
                            sum_cpue_index = TRUE,
                            # start at 1990 instead of 1984
                            start_year = 1990,
                            end_year = YEAR + 1,
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1, 2, 2, 2, 3, 3, 3)),
                            q_options = list(
                              pointer_biomass_cpue_strata = c(1, 1, 1, 2, 2, 2, 3, 3, 3),
                              pointer_q_cpue = c(1, 1, 1)))
m18s <- fit_rema(input)
out18s <- tidy_rema(m18s)
out18s$parameter_estimates

# Model 22 extra biomass + LLS RPW obs error -----
input <- prepare_rema_input(model_name = 'Model 22',
                            multi_survey = 1,
                            biomass_dat = biomass_dat,
                            cpue_dat = cpue_dat,
                            sum_cpue_index = TRUE,
                            end_year = YEAR + 1,
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1, 2, 2, 2, 3, 3, 3)),
                            q_options = list(
                              pointer_biomass_cpue_strata = c(1, 1, 1, 2, 2, 2, 3, 3, 3),
                              pointer_q_cpue = c(1, 1, 1)),
                            # ESTIMATE EXTRA biomass CV
                            extra_biomass_cv = list(assumption = 'extra_cv'),
                            # ESTIMATE EXTRA LLS RPW CV
                            extra_cpue_cv = list(assumption = 'extra_cv'))

m22 <- fit_rema(input)
out22 <- tidy_rema(m22)
out22$parameter_estimates

# Model 22 no 1984/87 -----
input <- prepare_rema_input(model_name = 'Model 22 no 1984/87',
                            multi_survey = 1,
                            biomass_dat = biomass_dat,
                            cpue_dat = cpue_dat,
                            sum_cpue_index = TRUE,
                            # start at 1990 instead of 1984
                            start_year = 1990,
                            end_year = YEAR + 1,
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1, 2, 2, 2, 3, 3, 3)),
                            q_options = list(
                              pointer_biomass_cpue_strata = c(1, 1, 1, 2, 2, 2, 3, 3, 3),
                              pointer_q_cpue = c(1, 1, 1)),
                            extra_biomass_cv = list(assumption = 'extra_cv'),
                            extra_cpue_cv = list(assumption = 'extra_cv'))

m22s <- fit_rema(input)
out22s <- tidy_rema(m22s)
out22s$parameter_estimates

# Compare M18* and M22 ----
compare <- compare_rema_models(rema_models = list(m18, m22))
compare$aic

cowplot::plot_grid(compare$plots$biomass_by_strata +
                     theme(legend.position = 'none') +
                     geom_line() +
                     labs(x = NULL, y = NULL, subtitle = 'Trawl survey biomass (t)',
                          fill = NULL, colour = NULL, shape = NULL, lty = NULL) +
                     ggplot2::scale_fill_viridis_d(direction = -1) +
                     ggplot2::scale_colour_viridis_d(direction = -1),
                   compare$plots$cpue_by_strata  +
                     facet_wrap(~strata, ncol = 1)  +
                     geom_line() +
                     labs(x = NULL, y = NULL, subtitle = 'Longline survey RPW',
                          fill = NULL, colour = NULL, shape = NULL, lty = NULL) +
                     ggplot2::scale_fill_viridis_d(direction = -1) +
                     ggplot2::scale_colour_viridis_d(direction = -1),
                   ncol = 2,
                   rel_widths = c(1.5, 1))

ggsave(filename = paste0(out_path, '/M18_M22_fits.png'),
       dpi = 400, bg = 'white', units = 'in', height = 9, width = 14)

compare$plots$total_predicted_biomass +
  labs(x = NULL, y = NULL, subtitle = 'Total predicted biomass (t)',
       fill = NULL, colour = NULL) +
  ggplot2::scale_fill_viridis_d(direction = -1) +
  ggplot2::scale_colour_viridis_d(direction = -1)

ggsave(filename = paste0(out_path, '/M18_M22_totalbiomass.png'),
       dpi = 400, bg = 'white', units = 'in', height = 3.5, width = 8)

# Compare M18* and M22 no 1984/87----
compare <- compare_rema_models(rema_models = list(m18s, m22s))
compare$aic %>% write_csv(paste0(out_path, '/m18s_m22s_aic.csv'))

cowplot::plot_grid(compare$plots$biomass_by_strata +
                     theme(legend.position = 'none') +
                     geom_line() +
                     labs(x = NULL, y = NULL, subtitle = 'Trawl survey biomass (t)',
                          fill = NULL, colour = NULL, shape = NULL, lty = NULL) +
                     ggplot2::scale_fill_viridis_d(direction = -1) +
                     ggplot2::scale_colour_viridis_d(direction = -1),
                   compare$plots$cpue_by_strata  +
                     facet_wrap(~strata, ncol = 1)  +
                     geom_line() +
                     labs(x = NULL, y = NULL, subtitle = 'Longline survey RPW',
                          fill = NULL, colour = NULL, shape = NULL, lty = NULL) +
                     ggplot2::scale_fill_viridis_d(direction = -1) +
                     ggplot2::scale_colour_viridis_d(direction = -1),
                   ncol = 2,
                   rel_widths = c(1.5, 1))

ggsave(filename = paste0(out_path, '/m18s_m22s_fits.png'),
       dpi = 400, bg = 'white', units = 'in', height = 9, width = 14)

compare$plots$total_predicted_biomass +
  labs(x = NULL, y = NULL, subtitle = 'Total predicted biomass (t)',
       fill = NULL, colour = NULL) +
  ggplot2::scale_fill_viridis_d(direction = -1) +
  ggplot2::scale_colour_viridis_d(direction = -1)

ggsave(filename = paste0(out_path, '/m18s_m22s_totalbiomass.png'),
       dpi = 400, bg = 'white', units = 'in', height = 3.5, width = 8)

# compare short and long time series ----
compare <- compare_rema_models(rema_models = list(m18, m18s, m22, m22s))

compare$plots$total_predicted_biomass +
  labs(x = NULL, y = NULL, subtitle = 'Total predicted biomass (t)',
       fill = NULL, colour = NULL) +
  ggplot2::scale_fill_viridis_d(direction = -1) +
  ggplot2::scale_colour_viridis_d(direction = -1)

ggsave(filename = paste0(out_path, '/M18_M18s_M22_M22s_totalbiomass.png'),
       dpi = 400, bg = 'white', units = 'in', height = 3.5, width = 8)

# param estimates
compare$output$parameter_estimates %>%
  write_csv(paste0(out_path, '/M18_M18s_M22_M22s_parameters.csv'))

# predicted biomass by strata and total for each model
biom <- compare$output$biomass_by_strata %>%
  pivot_wider(id_cols = c(strata, year), names_from = model_name, values_from = pred) %>%
  bind_rows(compare$output$total_predicted_biomass %>%
              mutate(strata = 'Total') %>%
              pivot_wider(id_cols = c(strata, year), names_from = model_name, values_from = pred)) %>%
  write_csv(paste0(out_path, '/M18_M18s_M22_M22s_biomass_pred.csv'))

# apportionment
appo <- compare$output$biomass_by_strata %>%
  mutate(strata = ifelse(grepl('CGOA', strata), 'CGOA',
                         ifelse(grepl('EGOA', strata), 'EGOA',
                                'WGOA'))) %>%
  group_by(model_name, strata, year) %>%
  summarize(stratum_biomass = sum(pred)) %>%
  group_by(model_name, year) %>%
  mutate(total_biomass = sum(stratum_biomass)) %>%
  ungroup() %>%
  mutate(proportion = stratum_biomass / total_biomass) %>%
  mutate(strata = factor(strata, labels = c('EGOA', 'CGOA', 'WGOA'), levels = c('EGOA', 'CGOA', 'WGOA'), ordered = TRUE)) %>%
  arrange(year, strata)

appo %>%
  pivot_wider(id_cols = c(strata, year), names_from = model_name, values_from = proportion) %>%
  write_csv(paste0(out_path, '/M18_M18s_M22_M22s_apportionment.csv'))

full_sumtable <- appo %>%
  filter(year == YEAR+1) %>%
  mutate(natmat = 0.03,
         OFL = natmat * total_biomass,
         maxABC = 0.75 * natmat * total_biomass,
         ABC = maxABC)

sumtable <- full_sumtable %>%
  distinct(model_name, year, biomass = total_biomass, OFL, maxABC) %>%
  select(model_name, year, biomass, OFL, maxABC) %>%
  write_csv(paste0(out_path, '/abc_ofl_summary.csv'))

# percent changes -----
biomass_dat %>% filter(year %in% c(2019,2021)) %>%
  mutate(strata = ifelse(grepl('CGOA', strata), 'CGOA',
                       ifelse(grepl('EGOA', strata), 'EGOA',
                              'WGOA'))) %>%
  group_by(year, strata) %>%
  summarise(biomass = sum(biomass, na.rm = T)) %>%
  pivot_wider(id_cols = c(strata), names_from = year, values_from = biomass) %>%
  mutate(percent_change = (`2021`-`2019`)/`2019`)

cpue_dat %>% filter(year %in% c(2020,2021)) %>%
  pivot_wider(id_cols = c(strata), names_from = year, values_from = cpue) %>%
  mutate(percent_change = (`2021`-`2020`)/`2020`)

old <- out1$total_predicted_biomass %>% filter(year == 2022) %>% pull(pred)
new <- out2$total_predicted_biomass %>% filter(year == 2022) %>% pull(pred)
(new-old)/old
