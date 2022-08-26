# This script compares the following models using the data from the 2020 GOA
# thornyheads assessment (using data as presented in 2020):
# 1. Status quo (Model 20): this is the incorrect version of ADMB where
# predicted CPUE is defined outside of the SEPARABLE_FUNCTION.
# 2. Model 22.1.a:  the corrected version of Model 20 in TMB

# Set up ----

libs <- c('readr', 'dplyr', 'tidyr', 'ggplot2', 'cowplot')
if(length(libs[which(libs %in% rownames(installed.packages()) == FALSE )]) > 0) {install.packages(libs[which(libs %in% rownames(installed.packages()) == FALSE)])}
lapply(libs, library, character.only = TRUE)

# install.packages("devtools")
# devtools::install_github("afsc-assessments/rema", dependencies = TRUE, build_vignettes = TRUE)
library(rema)

# folder set up
YEAR <- 2022 # assessment year
dat_path <- paste0("data/", YEAR); dir.create(dat_path)
out_path <- paste0("results/", YEAR); dir.create(out_path)

ggplot2::theme_set(cowplot::theme_cowplot(font_size = 12) +
                     cowplot::background_grid() +
                     cowplot::panel_border())

# ADMB models ----

# Incorrect version
oldadmb <- read_admb_re(filename = 'admb_original/RWOUT.REP',
                        model_name = 'Model 20',
                        biomass_strata_names = c('CGOA (0-500 m)', 'CGOA (501-700 m)', 'CGOA (701-1000 m)',
                                                 'EGOA (0-500 m)', 'EGOA (501-700 m)', 'EGOA (701-1000 m)',
                                                 'WGOA (0-500 m)', 'WGOA (501-700 m)', 'WGOA (701-1000 m)'),
                        cpue_strata_names = c('CGOA', 'EGOA', 'WGOA'))
oldadmb$biomass_dat %>% filter(year == 1987)

# Corrected version - this model structure was too complicated to easily fix the
# error.
# fixedadmb <- read_admb_re(filename = 'admb_fixed/RWOUT.REP',
#                         model_name = 'Model 22.1.a',
#                         biomass_strata_names = c('CGOA (0-500 m)', 'CGOA (501-700 m)', 'CGOA (701-1000 m)',
#                                                  'EGOA (0-500 m)', 'EGOA (501-700 m)', 'EGOA (701-1000 m)',
#                                                  'WGOA (0-500 m)', 'WGOA (501-700 m)', 'WGOA (701-1000 m)'),
#                         cpue_strata_names = c('CGOA', 'EGOA', 'WGOA'))

fixedadmb$biomass_dat %>% print(n=Inf)

# TMB model with zeros as 0.0001 CV = 1000 ----
input <- prepare_rema_input(model_name = 'Model 22.1.a',
                            multi_survey = 1,
                            admb_re = oldadmb,
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

input$biomass_dat %>% filter(year == 1987)

m1 <- fit_rema(input)
out1 <- tidy_rema(m1)
out1$parameter_estimates

# Compare models ----

compare <- compare_rema_models(rema_models = list(m1),
                               admb_re = oldadmb)
compare2 <- compare_rema_models(rema_models = list(m1))

biomass_by_strata <- compare$output$biomass_by_strata
total_predicted_biomass <- compare$output$total_predicted_biomass

p1 <- ggplot(data = biomass_by_strata,
             aes(x = year, y = pred,
                 col = model_name)) +
  geom_ribbon(aes(ymin = pred_lci, ymax = pred_uci,
                  fill = model_name), col = NA,
              alpha = 0.25) +
  geom_line(aes(lty = model_name)) +
  facet_wrap(~strata, nrow = NULL) +
  geom_point(data = out1$biomass_by_strata, aes(x = year, y = obs), col = 'black') +
  geom_errorbar(data = out1$biomass_by_strata,
                aes(x = year, ymin = obs_lci, ymax = obs_uci),
                col = 'black', lty = 1) +
  scale_y_continuous(labels = scales::comma, expand = c(0.01, 0), limits = c(0, NA)) +
  labs(x = NULL, y = NULL, subtitle = 'Trawl survey biomass (t)',
       fill = NULL, colour = NULL, shape = NULL, lty = NULL) +
  ggplot2::scale_fill_viridis_d(direction = -1) +
  ggplot2::scale_colour_viridis_d(direction = -1) +
  theme(legend.position = 'left')

p1

p2 <- compare2$plots$cpue_by_strata +
  facet_wrap(~strata, ncol = 1) +
  theme(legend.position = 'none') +
  labs(y = NULL, subtitle = 'Relative Population Weights')

plot_grid(p1, p2, ncol = 2, rel_widths = c(0.7, 0.3))

ggsave(filename = paste0(out_path, '/M20_M22.1.a_fits.png'),
       dpi = 400, bg = 'white', units = 'in', height = 9, width = 14)

ggplot(data = total_predicted_biomass,
       aes(x = year, y = pred,
           col = model_name)) +
  geom_ribbon(aes(ymin = pred_lci, ymax = pred_uci,
                  fill = model_name), col = NA,
              alpha = 0.25) +
  geom_line() +
  scale_y_continuous(labels = scales::comma) + #, expand = c(0, 0), limits = c(0, NA)) +
  labs(x = NULL, y = NULL, subtitle = 'Total predicted biomass (t)',
       fill = NULL, colour = NULL) +
  ggplot2::scale_fill_viridis_d(direction = -1) +
  ggplot2::scale_colour_viridis_d(direction = -1)

ggsave(filename = paste0(out_path, '/M20_M22.1.a_totalbiomass.png'),
       dpi = 400, bg = 'white', units = 'in', height = 4, width = 7)
