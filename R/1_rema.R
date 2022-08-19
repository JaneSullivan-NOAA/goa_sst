# GOA SST biomass estimation using the bottom trawl and longline survey indices

# Set up ----

# assessment year
YEAR <- 2022

libs <- c('readr', 'dplyr', 'tidyr', 'ggplot2', 'cowplot')
if(length(libs[which(libs %in% rownames(installed.packages()) == FALSE )]) > 0) {install.packages(libs[which(libs %in% rownames(installed.packages()) == FALSE)])}
lapply(libs, library, character.only = TRUE)

# install.packages("devtools")
devtools::install_github("JaneSullivan-NOAA/rema", dependencies = TRUE)
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
admb_re <- read_admb_re(filename = paste0(dat_path, '/goasst_rwout_2020.rep'),
                        biomass_strata_names = c('CGOA (0-500 m)', 'CGOA (501-700 m)', 'CGOA (701-1000 m)',
                                                 'EGOA (0-500 m)', 'EGOA (501-700 m)', 'EGOA (701-1000 m)',
                                                 'WGOA (0-500 m)', 'WGOA (501-700 m)', 'WGOA (701-1000 m)'),
                        cpue_strata_names = c('CGOA', 'EGOA', 'WGOA'),
                        model_name = 'ADMB M1: 2020 single q, q.wt=0.5*')

# check new and old data
biomass_dat %>%
  rename(new_biomass = biomass, new_cv = cv) %>%
  left_join(admb_re$biomass_dat) %>%
  mutate(tst_biomass = new_biomass - biomass,
         tst_cv = new_cv - cv) %>%
  filter(abs(tst_biomass) > 0.001 | abs(tst_cv) > 0.001)
# strata               year new_biomass new_cv  biomass cv       tst_biomass tst_cv
# 1 CGOA (501-700 m)   2019       6015. 0.177    6015. 0.167               0 0.0102
# 2 EGOA (501-700 m)   1984       3639. 0.103    3639. 0.0652              0 0.0380
# 3 EGOA (701-1000 m)  1984        814  0.1       814  0.000123            0 0.0999
# 4 EGOA (701-1000 m)  2009       4821. 0.0921   4821. 0.0187              0 0.0734
# 5 EGOA (701-1000 m)  2015       3686. 0.0906   3686. 0.00533             0 0.0853
# 6 WGOA (701-1000 m)  1999       1679  0.1      1679  0.0000596           0 0.0999
# bts %>% filter(year == 2019 & strata == 'CGOA (501-700 m)') %>%
#   group_by(regulatory_area_name) %>%
#   summarize(biomass = sum(area_biomass),
#             sd = sqrt(sum(var)))

# M1: TMB GOA SST single/shared q, q_wt = 0.5 -----
input <- prepare_rema_input(model_name = 'M1: single q, q_wt=0.5',
                            multi_survey = 1,
                            # admb_re = admb_re,
                            biomass_dat = biomass_dat,
                            cpue_dat = cpue_dat,
                            # likelihood weight
                            wt_cpue = 0.5,
                            # is the CPUE index summable? RPNs and RPWs are
                            # summable, but raw CPUE is not
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
                             zeros = list(assumption = 'NA'))

m1 <- fit_rema(input)
tidy_rema(rema_model = m1)$parameter_estimates

# M2: TMB GOA SST single q, q_wt = 1 -----
input <- prepare_rema_input(model_name = 'M2: single q, q_wt=1',
                            multi_survey = 1,
                            # admb_re = admb_re,
                            biomass_dat = biomass_dat,
                            cpue_dat = cpue_dat,
                            wt_cpue = 1, # CHANGE LIKELIHOOD WEIGHT
                            sum_cpue_index = TRUE,
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1, 2, 2, 2, 3, 3, 3)),
                            q_options = list(
                              pointer_biomass_cpue_strata = c(1, 1, 1, 2, 2, 2, 3, 3, 3),
                              pointer_q_cpue = c(1, 1, 1)),
                            zeros = list(assumption = 'NA'))

m2 <- fit_rema(input)
tidy_rema(m2)$parameter_estimates

# M3: TMB GOA SST strata q, q_wt = 0.5 -----
input <- prepare_rema_input(model_name = 'M3: strata q, q_wt=0.5',
                            multi_survey = 1,
                            # admb_re = admb_re,
                            biomass_dat = biomass_dat,
                            cpue_dat = cpue_dat,
                            wt_cpue = 0.5,
                            sum_cpue_index = TRUE,
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1, 2, 2, 2, 3, 3, 3)),
                            q_options = list(
                              pointer_biomass_cpue_strata = c(1, 1, 1, 2, 2, 2, 3, 3, 3),
                              pointer_q_cpue = c(1, 2, 3)), # STRATA-SPECIFIC q
                            zeros = list(assumption = 'NA'))

m3 <- fit_rema(input)
tidy_rema(m3)$parameter_estimates

# M4: TMB GOA SST strata q, q_wt = 1 -----

input <- prepare_rema_input(model_name = 'M4: strata q, q_wt=1',
                            multi_survey = 1,
                            # admb_re = admb_re,
                            biomass_dat = biomass_dat,
                            cpue_dat = cpue_dat,
                            wt_cpue = 1, # CHANGE LIKELIHOOD WEIGHT
                            sum_cpue_index = TRUE,
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1, 2, 2, 2, 3, 3, 3)),
                            q_options = list(
                              pointer_biomass_cpue_strata = c(1, 1, 1, 2, 2, 2, 3, 3, 3),
                              pointer_q_cpue = c(1, 2, 3)), # STRATA-SPECIFIC q
                            zeros = list(assumption = 'NA'))

m4 <- fit_rema(input)
tidy_rema(m4)$parameter_estimates

# compare weights ----
compare <- compare_rema_models(rema_models = list(m1, m2, m3, m4),
                               biomass_ylab = 'Trawl survey biomass (t)',
                               cpue_ylab = 'Longline survey RPW')
compare$aic
cowplot::plot_grid(compare$plots$biomass_by_strata +
                     theme(legend.position = 'left'),
                   compare$plots$cpue_by_strata  +
                     facet_wrap(~strata, ncol = 1)  +
                     theme(legend.position = 'none'),
                   ncol = 2,
                   rel_widths = c(2, 1))

# M5: Estimate extra variance for biomass obs ----

input <- prepare_rema_input(model_name = 'M5: M2 + xtra biomass CV',
                            multi_survey = 1,
                            # admb_re = admb_re,
                            biomass_dat = biomass_dat,
                            cpue_dat = cpue_dat,
                            wt_cpue = 1,
                            sum_cpue_index = TRUE,
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1, 2, 2, 2, 3, 3, 3)),
                            q_options = list(
                              pointer_biomass_cpue_strata = c(1, 1, 1, 2, 2, 2, 3, 3, 3),
                              pointer_q_cpue = c(1, 1, 1)),
                            zeros = list(assumption = 'NA'),
                            # ESTIMATE EXTRA BIOMASS CV
                            extra_biomass_cv = list(assumption = 'extra_cv'))

m5 <- fit_rema(input)
tidy_rema(m5)$parameter_estimates

# M6: Estimate extra variance for CPUE obs ----

input <- prepare_rema_input(model_name = 'M6: M2 + xtra RPW CV',
                            multi_survey = 1,
                            # admb_re = admb_re,
                            biomass_dat = biomass_dat,
                            cpue_dat = cpue_dat,
                            wt_cpue = 1,
                            sum_cpue_index = TRUE,
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1, 2, 2, 2, 3, 3, 3)),
                            q_options = list(
                              pointer_biomass_cpue_strata = c(1, 1, 1, 2, 2, 2, 3, 3, 3),
                              pointer_q_cpue = c(1, 1, 1)),
                            zeros = list(assumption = 'NA'),
                            # ESTIMATE EXTRA RPW CV
                            extra_cpue_cv = list(assumption = 'extra_cv'))

m6 <- fit_rema(input)
tidy_rema(m6)$parameter_estimates

# M7: Estimate extra variance for biomass + CPUE obs ----

input <- prepare_rema_input(model_name = 'M7: M2 + xtra biomass & RPW CV',
                            multi_survey = 1,
                            # admb_re = admb_re,
                            biomass_dat = biomass_dat,
                            cpue_dat = cpue_dat,
                            wt_cpue = 1,
                            sum_cpue_index = TRUE,
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1, 2, 2, 2, 3, 3, 3)),
                            q_options = list(
                              pointer_biomass_cpue_strata = c(1, 1, 1, 2, 2, 2, 3, 3, 3),
                              pointer_q_cpue = c(1, 1, 1)),
                            zeros = list(assumption = 'NA'),
                            # ESTIMATE EXTRA biomass CV
                            extra_biomass_cv = list(assumption = 'extra_cv'),
                            # ESTIMATE EXTRA RPW CV
                            extra_cpue_cv = list(assumption = 'extra_cv'))

m7 <- fit_rema(input)
tidy_rema(m7)$parameter_estimates

# Compare new ----
compare <- compare_rema_models(rema_models = list(m2, m5, m6, m7),
                               biomass_ylab = 'Biomass (t)',
                               cpue_ylab = 'Relative Population Weights')

compare$aic
cowplot::plot_grid(compare$plots$biomass_by_strata +
                     theme(legend.position = 'left'),
                   compare$plots$cpue_by_strata  +
                     facet_wrap(~strata, ncol = 1)  +
                     theme(legend.position = 'none'),
                   ncol = 2,
                   rel_widths = c(2, 1))

compare$output$parameter_estimates
compare$plots$total_predicted_biomass
compare$plots$total_predicted_cpue

# compare old ----
compare <- compare_rema_models(rema_models = list(m2, m7),
                               admb_re = admb_re,
                               biomass_ylab = 'Biomass (t)',
                               cpue_ylab = 'Relative Population Weights')

cowplot::plot_grid(compare$plots$biomass_by_strata +
                     theme(legend.position = 'left'),
                   compare$plots$cpue_by_strata  +
                     facet_wrap(~strata, ncol = 1)  +
                     theme(legend.position = 'none'),
                   ncol = 2,
                   rel_widths = c(2, 1))

compare$plots$total_predicted_biomass
compare$plots$total_predicted_cpue
