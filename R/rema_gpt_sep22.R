# GOA SST biomass estimation using the bottom trawl and longline survey indices

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

# longline survey rpws
cpue_dat <- read_csv(paste0(dat_path, "/goa_sst_rpw_", YEAR, ".csv"))

# SST rwout.rep data ----

# Multi-survey and multi-strata version of the random effects model (REMA).
# Example using GOA shortspine thornyhead, which has different strata
# definitions for the biomass and CPUE surveys.
admb_re <- read_admb_re(filename = 'admb_original/RWOUT.REP',
                        model_name = 'Model 18',
                        biomass_strata_names = c('CGOA (0-500 m)', 'CGOA (501-700 m)', 'CGOA (701-1000 m)',
                                                 'EGOA (0-500 m)', 'EGOA (501-700 m)', 'EGOA (701-1000 m)',
                                                 'WGOA (0-500 m)', 'WGOA (501-700 m)', 'WGOA (701-1000 m)'),
                        cpue_strata_names = c('CGOA', 'EGOA', 'WGOA'))

# check to make sure new and old data match up
compare_data <- biomass_dat %>%
  rename(new_biomass = biomass, new_cv = cv) %>%
  left_join(admb_re$biomass_dat) %>%
  mutate(tst_biomass = new_biomass - biomass,
         tst_cv = new_cv - cv)
write_csv(compare_data, paste0(out_path, '/changes_in_biomass_data.csv'))
compare_data %>%
  filter(abs(tst_biomass) > 0.001 | abs(tst_cv) > 0.001)
# There are 6 data points where the CVs don't match up. Note sure why!
# strata               year new_biomass new_cv  biomass cv       tst_biomass tst_cv
# 1 CGOA (501-700 m)   2019       6015. 0.177    6015. 0.167               0 0.0102
# 2 EGOA (501-700 m)   1984       3639. 0.103    3639. 0.0652              0 0.0380
# 3 EGOA (701-1000 m)  1984        814  0.1       814  0.000123            0 0.0999
# 4 EGOA (701-1000 m)  2009       4821. 0.0921   4821. 0.0187              0 0.0734
# 5 EGOA (701-1000 m)  2015       3686. 0.0906   3686. 0.00533             0 0.0853
# 6 WGOA (701-1000 m)  1999       1679  0.1      1679  0.0000596           0 0.0999
compare_data <- cpue_dat %>%
  rename(new_cpue = cpue, new_cv = cv) %>%
  left_join(admb_re$cpue_dat) %>%
  mutate(tst_cpue = new_cpue - cpue,
         tst_cv = new_cv - cv)
write_csv(compare_data, paste0(out_path, '/changes_in_cpue_data.csv'))
compare_data %>%
  filter(abs(tst_cpue) > 0.001 | abs(tst_cv) > 0.001)

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

# Model 22.1.a corrected Model 18 in TMB with 2020 data ----
input <- prepare_rema_input(model_name = 'Model 22.1.a',
                            multi_survey = 1,
                            admb_re = admb_re,
                            # RPWs are summable
                            sum_cpue_index = TRUE,
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
                              pointer_q_cpue = c(1, 1, 1)),
                            # assumption for how to treat zeros
                            zeros = list(assumption = 'small_constant',
                                         options_small_constant = c(0.0001, 1000)))
m1 <- fit_rema(input)
out1 <- tidy_rema(m1)
out1$parameter_estimates

# Model 22.1.b fixed Model 18 in TMB with updated 2021 data -----
input <- prepare_rema_input(model_name = 'Model 22.1.b',
                            multi_survey = 1,
                            biomass_dat = biomass_dat,
                            cpue_dat = cpue_dat,
                            end_year = YEAR,
                            # RPWs are summable
                            sum_cpue_index = TRUE,
                            # three process error parameters (log_PE) estimated,
                            # indexed as follows for each biomass survey stratum
                            # (shared within an area across depths):
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1, 2, 2, 2, 3, 3, 3)),
                            q_options = list(
                              # CPUE strata (n=3) indexed as follows for the
                              # biomass strata (n=9)
                              pointer_biomass_cpue_strata = c(1, 1, 1, 2, 2, 2, 3, 3, 3),
                              # one scaling parameters (log_q) estimated, shared
                              # over all three strata
                              pointer_q_cpue = c(1, 1, 1)),
                            # assumption for how to treat zeros
                            zeros = list(assumption = 'small_constant',
                                         options_small_constant = c(0.0001, 1000)))

m2 <- fit_rema(input)
out2 <- tidy_rema(m2)
out2$parameter_estimates

# Model 22.2.a = Model 22.1.b with extra biomass obs error -----
input <- prepare_rema_input(model_name = 'Model 22.2.a',
                            multi_survey = 1,
                            biomass_dat = biomass_dat,
                            cpue_dat = cpue_dat,
                            end_year = YEAR,
                            sum_cpue_index = TRUE,
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1, 2, 2, 2, 3, 3, 3)),
                            q_options = list(
                              pointer_biomass_cpue_strata = c(1, 1, 1, 2, 2, 2, 3, 3, 3),
                              pointer_q_cpue = c(1, 1, 1)),
                            zeros = list(assumption = 'small_constant',
                                         options_small_constant = c(0.0001, 1000)),
                            # ESTIMATE EXTRA BIOMASS CV
                            extra_biomass_cv = list(assumption = 'extra_cv'))

m3 <- fit_rema(input)
out3 <- tidy_rema(m3)
out3$parameter_estimates

# Model 22.2.b = Model 22.1.b with extra LLS RPW obs error -----
input <- prepare_rema_input(model_name = 'Model 22.2.b',
                            multi_survey = 1,
                            biomass_dat = biomass_dat,
                            cpue_dat = cpue_dat,
                            end_year = YEAR,
                            sum_cpue_index = TRUE,
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1, 2, 2, 2, 3, 3, 3)),
                            q_options = list(
                              pointer_biomass_cpue_strata = c(1, 1, 1, 2, 2, 2, 3, 3, 3),
                              pointer_q_cpue = c(1, 1, 1)),
                            zeros = list(assumption = 'small_constant',
                                         options_small_constant = c(0.0001, 1000)),
                            # ESTIMATE EXTRA LLS RPW CV
                            extra_cpue_cv = list(assumption = 'extra_cv'))

m4 <- fit_rema(input)
out4 <- tidy_rema(m4)
out4$parameter_estimates

# Model 22.2.c = Model 22.1.b with extra biomass + LLS RPW obs error -----
input <- prepare_rema_input(model_name = 'Model 22.2.c',
                            multi_survey = 1,
                            biomass_dat = biomass_dat,
                            cpue_dat = cpue_dat,
                            end_year = YEAR,
                            sum_cpue_index = TRUE,
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1, 2, 2, 2, 3, 3, 3)),
                            q_options = list(
                              pointer_biomass_cpue_strata = c(1, 1, 1, 2, 2, 2, 3, 3, 3),
                              pointer_q_cpue = c(1, 1, 1)),
                            zeros = list(assumption = 'small_constant',
                                         options_small_constant = c(0.0001, 1000)),
                            # ESTIMATE EXTRA biomass CV
                            extra_biomass_cv = list(assumption = 'extra_cv'),
                            # ESTIMATE EXTRA LLS RPW CV
                            extra_cpue_cv = list(assumption = 'extra_cv'))

m5 <- fit_rema(input)
out5 <- tidy_rema(m5)
out5$parameter_estimates

# Model 22.3 no LLS ----
input <- prepare_rema_input(model_name = 'Model 22.3',
                            multi_survey = 0,
                            end_year = YEAR,
                            biomass_dat = biomass_dat,
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1, 2, 2, 2, 3, 3, 3)),
                            zeros = list(assumption = 'small_constant',
                                         options_small_constant = c(0.0001, 1000)))

m6 <- fit_rema(input)
out6 <- tidy_rema(m6)
out6$parameter_estimates

# Compare data changes ----
compare <- compare_rema_models(rema_models = list(m1, m2))

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

ggsave(filename = paste0(out_path, '/M22.1.a_M22.1.b_fits.png'),
       dpi = 400, bg = 'white', units = 'in', height = 9, width = 14)

compare$plots$total_predicted_biomass +
  labs(x = NULL, y = NULL, subtitle = 'Total predicted biomass (t)',
       fill = NULL, colour = NULL) +
  ggplot2::scale_fill_viridis_d(direction = -1) +
  ggplot2::scale_colour_viridis_d(direction = -1)

ggsave(filename = paste0(out_path, '/M22.1.a_M22.1.b_totalbiomass.png'),
       dpi = 400, bg = 'white', units = 'in', height = 3.5, width = 8)

# compare new ----

compare <- compare_rema_models(rema_models = list(m2, m3, m4, m5),
                               biomass_ylab = 'Biomass (t)',
                               cpue_ylab = 'Relative Population Weights')
compare$aic %>% write_csv(paste0(out_path, '/aic.csv'))

cowplot::plot_grid(compare$plots$biomass_by_strata +
                     theme(legend.position = 'none') +
                     geom_line(size = 0.8) +
                     labs(x = NULL, y = NULL, subtitle = 'Trawl survey biomass (t)',
                          fill = NULL, colour = NULL, shape = NULL, lty = NULL),
                   compare$plots$cpue_by_strata  +
                     facet_wrap(~strata, ncol = 1)  +
                     geom_line(size = 0.8) +
                     labs(x = NULL, y = NULL, subtitle = 'Longline survey RPW',
                          fill = NULL, colour = NULL, shape = NULL, lty = NULL),
                   ncol = 2,
                   rel_widths = c(1.5, 1))

ggsave(filename = paste0(out_path, '/M22.1.b_M22.2.abc_fits.png'),
       dpi = 400, bg = 'white', units = 'in', height = 9, width = 14)

compare$plots$total_predicted_biomass +
  labs(x = NULL, y = NULL, subtitle = 'Total predicted biomass (t)',
       fill = NULL, colour = NULL)

ggsave(filename = paste0(out_path, '/M22.1.b_M22.2.abc_totalbiomass.png'),
       dpi = 400, bg = 'white', units = 'in', height = 3.5, width = 8)

# Top models ----

compare <- compare_rema_models(rema_models = list(m2, m5, m6),
                               biomass_ylab = 'Biomass (t)',
                               cpue_ylab = 'Relative Population Weights')
compare$plots$biomass_by_strata +
  geom_line(size = 0.8) +
  labs(x = NULL, y = NULL, subtitle = 'Trawl survey biomass (t)',
       fill = NULL, colour = NULL, shape = NULL, lty = NULL)

ggsave(filename = paste0(out_path, '/M22.1.b_M22.2.c_M22.3_fits.png'),
       dpi = 400, bg = 'white', units = 'in', height = 9, width = 14)

compare$plots$total_predicted_biomass +
  labs(x = NULL, y = NULL, subtitle = 'Total predicted biomass (t)',
       fill = NULL, colour = NULL)
ggsave(filename = paste0(out_path, '/M22.1.b_M22.2.c_M22.3_totalbiomass.png'),
       dpi = 400, bg = 'white', units = 'in', height = 3.5, width = 8)

# top models with the old ADMB model ----
compare <- compare_rema_models(rema_models = list(m2, m5, m6),
                               admb_re = admb_re)
compare$plots$total_predicted_biomass +
  labs(x = NULL, y = NULL, subtitle = 'Total predicted biomass (t)',
       fill = NULL, colour = NULL) +
  geom_line(size = 0.8)
ggsave(filename = paste0(out_path, '/M18_M22.1.b_M22.2.c_M22.3_totalbiomass.png'),
       dpi = 400, bg = 'white', units = 'in', height = 3.5, width = 8)

# save output ----

compare <- compare_rema_models(list(m1, m2, m3, m4, m5, m6))

# model 20 parameter estimates (pulled from re.std in original_admb)
data.frame(model_name = 'Model 18',
           parameter = c('log_q', 'log_pe_cgoa', 'log_pe_egoa', 'log_pe_wgoa'),
           est = c(-4.8615e-01, -2.6580e+00, -1.7364e+00, -2.2476e+00),
           se = c(1.1597e-02, 2.1528e-01, 1.7461e-01, 1.9868e-01)) %>%
  mutate(estimate = exp(est),
         std_err = exp(est) * se,
         lci = exp(est - qnorm(1 - 0.05/2) * se),
         uci = exp(est + qnorm(1 - 0.05/2) * se)) %>%
  select(-est, - se) %>%
  write_csv(paste0(out_path, '/M18_parameters.csv'))

compare$output$parameter_estimates %>%
  write_csv(paste0(out_path, '/M22.1.ab_M22.2.abc_M22.3_paramaters.csv'))

compare <- compare_rema_models(list(m1, m2, m3, m4, m5, m6),
                               admb_re = admb_re)

biom <- compare$output$biomass_by_strata %>%
  pivot_wider(id_cols = c(strata, year), names_from = model_name, values_from = pred) %>%
  bind_rows(compare$output$total_predicted_biomass %>%
              mutate(strata = 'Total') %>%
              pivot_wider(id_cols = c(strata, year), names_from = model_name, values_from = pred)) %>%
  write_csv(paste0(out_path, '/M18_M22.1.ab_M22.2.abc_M22.3_biomass_pred.csv'))

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
  write_csv(paste0(out_path, '/M18_M22.1.ab_M22.2.abc_M22.3_apportionment.csv'))

full_sumtable <- appo %>%
  filter(year == YEAR) %>%
  mutate(natmat = 0.03,
         OFL = natmat * total_biomass,
         maxABC = 0.75 * natmat * total_biomass,
         ABC = maxABC)

sumtable <- full_sumtable %>%
  distinct(model_name, year, biomass = total_biomass, OFL, maxABC) %>%
  select(model_name, year, biomass, OFL, maxABC)

statquo_biomass <- sumtable %>% filter(model_name == 'Model 18') %>% pull(biomass)
# statquo_OFL <- sumtable %>% filter(model_name == 'Model 18') %>% pull(OFL)
# statquo_maxABC <- sumtable %>% filter(model_name == 'Model 18') %>% pull(maxABC)

sumtable %>%
  filter(model_name %in% c('Model 18', 'Model 22.1.a')) %>%
  mutate(year = 2020,
         percent_change = (biomass - statquo_biomass)/statquo_biomass * 100) %>%
  write_csv(paste0(out_path, '/abc_ofl_summary_2020.csv'))

statquo_biomass <- sumtable %>% filter(model_name == 'Model 22.1.b') %>% pull(biomass)

sumtable %>%
  filter(!model_name %in% c('Model 18', 'Model 22.1.a')) %>%
  mutate(percent_change = (biomass - statquo_biomass)/statquo_biomass * 100) %>%
  write_csv(paste0(out_path, '/abc_ofl_summary_2022.csv'))

full_sumtable %>%
  mutate(year = ifelse(model_name %in% c('Model 18', 'Model 22.1.a'), 2020, year)) %>%
  pivot_wider(id_cols = c(model_name, year), names_from = strata, values_from = proportion) %>%
  write_csv(paste0(out_path, '/apportionment_summary.csv'))

compare <- compare_rema_models(list(m2, m3, m4, m5))

compare$output$cpue_by_strata %>%
  mutate(pred = prettyNum(pred, big.mark = ',', digits = 1, trim = TRUE)) %>%
  pivot_wider(id_cols = c(strata, year), names_from = model_name, values_from = pred) %>%
  mutate(`Model 22.3` = NA) %>%
  bind_rows(compare$output$total_predicted_cpue %>%
              mutate(strata = 'Total') %>%
              mutate(pred = prettyNum(pred, big.mark = ',', digits = 1, trim = TRUE)) %>%
              pivot_wider(id_cols = c(strata, year), names_from = model_name, values_from = pred) %>%
              mutate(`Model 22.3` = NA) ) %>%
  write_csv(paste0(out_path, '/M22.1.b_M22.2.abc_M22.3_cpue_pred.csv'))

compare$aic %>%
  mutate(objective_function = prettyNum(objective_function, digits = 3)) %>%
  write_csv(paste0(out_path, '/aic.csv'))

# Output for best model by AIC ----
out5$biomass_by_strata %>%
  mutate(pred = prettyNum(pred, big.mark = ',', digits = 1, trim = TRUE),
         pred_lci = prettyNum(pred_lci, big.mark = ',', digits = 1, trim = TRUE),
         pred_uci = prettyNum(pred_uci, big.mark = ',', digits = 1, trim = TRUE),
         pred = paste0(pred, " (", pred_lci, ", ", pred_uci, ")"),
         obs = ifelse(is.na(obs), NA,
                      paste0(prettyNum(obs, big.mark = ',', digits = 1, trim = TRUE),
                      " (", round(obs_cv), ")"))) %>%
  select(model_name, strata, year, pred, obs) %>%
  bind_rows(compare$output$total_predicted_biomass %>%
              mutate(strata = 'Total',red = prettyNum(pred, big.mark = ',', digits = 1, trim = TRUE),
                     pred_lci = prettyNum(pred_lci, big.mark = ',', digits = 1, trim = TRUE),
                     pred_uci = prettyNum(pred_uci, big.mark = ',', digits = 1, trim = TRUE),
                     pred = paste0(pred, " (", pred_lci, ", ", pred_uci, ")"),
                     obs = NA) %>%
              select(model_name, strata, year, pred, obs)) %>%
  write_csv(paste0(out_path, '/topmod_M22.2.c_biomass.csv'))

out5$cpue_by_strata %>%
  mutate(pred = prettyNum(pred, big.mark = ',', digits = 1, trim = TRUE),
         pred_lci = prettyNum(pred_lci, big.mark = ',', digits = 1, trim = TRUE),
         pred_uci = prettyNum(pred_uci, big.mark = ',', digits = 1, trim = TRUE),
         pred = paste0(pred, " (", pred_lci, ", ", pred_uci, ")"),
         obs = ifelse(is.na(obs), NA,
                      paste0(prettyNum(obs, big.mark = ',', digits = 1, trim = TRUE),
                             " (", round(obs_cv), ")"))) %>%
  select(model_name, strata, year, pred, obs) %>%
  bind_rows(compare$output$total_predicted_cpue %>%
              mutate(strata = 'Total',red = prettyNum(pred, big.mark = ',', digits = 1, trim = TRUE),
                     pred_lci = prettyNum(pred_lci, big.mark = ',', digits = 1, trim = TRUE),
                     pred_uci = prettyNum(pred_uci, big.mark = ',', digits = 1, trim = TRUE),
                     pred = paste0(pred, " (", pred_lci, ", ", pred_uci, ")"),
                     obs = NA) %>%
              select(model_name, strata, year, pred, obs)) %>%
  write_csv(paste0(out_path, '/topmod_M22.2.c_cpue.csv'))

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

out2$total_predicted_biomass %>%
  filter(year %in% c(2020,2021)) %>%
  pivot_wider(id_cols = model_name, names_from = year, values_from = pred) %>%
  mutate(percent_change = (`2021`-`2020`)/`2020`)

# sensitivity of 1984/1987 surveys ----

input <- prepare_rema_input(model_name = 'Model 22.2.c no 1984/87 BTS',
                            multi_survey = 1,
                            biomass_dat = biomass_dat %>%
                              filter(!year %in% c(1984, 1987)),
                            cpue_dat = cpue_dat,
                            # start_year = 1984,
                            end_year = YEAR,
                            sum_cpue_index = TRUE,
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1, 2, 2, 2, 3, 3, 3)),
                            q_options = list(
                              pointer_biomass_cpue_strata = c(1, 1, 1, 2, 2, 2, 3, 3, 3),
                              pointer_q_cpue = c(1, 1, 1)),
                            zeros = list(assumption = 'small_constant',
                                         options_small_constant = c(0.0001, 1000)),
                            # ESTIMATE EXTRA biomass CV
                            extra_biomass_cv = list(assumption = 'extra_cv'),
                            # ESTIMATE EXTRA LLS RPW CV
                            extra_cpue_cv = list(assumption = 'extra_cv'))

m5_no1980s <- fit_rema(input)
outm5_no1980s <- tidy_rema(m5_no1980s)
outm5_no1980s$parameter_estimates

compare <- compare_rema_models(list(m5, m5_no1980s))
compare$plots$biomass_by_strata
compare$plots$total_predicted_biomass
compare$output$parameter_estimates
