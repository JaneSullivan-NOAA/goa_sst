# GOA shortspine thornyhead queries
# RACE species code 30020

# set up ----

# assessment year
YEAR <- 2022

libs <- c('readr', 'dplyr', 'tidyr', 'RODBC', 'ggplot2')
if(length(libs[which(libs %in% rownames(installed.packages()) == FALSE )]) > 0) {install.packages(libs[which(libs %in% rownames(installed.packages()) == FALSE)])}
lapply(libs, library, character.only = TRUE)

# folder set up
dat_path <- paste0("data/", YEAR); dir.create(dat_path)
raw_path <- paste0(dat_path, "/raw"); dir.create(raw_path)
out_path <- paste0("results/", YEAR); dir.create(out_path)

# database connection ----

# Enter your username and password for the AKFIN database. Note that these
# credentials are different than what you may use to access AKFIN Answer.
# Contact AKFIN for more information.
username_akfin = 'my_username'
password_akfin = 'my_password'
username_akfin = 'jsullivan'
password_akfin = 'sculja22'
channel_akfin <- odbcConnect("akfin", uid = username_akfin, pwd = password_akfin, believeNRows=FALSE)

# INPFC and depth area look up ----

# INPFC = International North Pacific Fisheries Commission

query <- "select distinct   survey, regulatory_area_name, inpfc_area,
                            summary_area_depth, min_depth, max_depth, area
          from              afsc.race_goastrataaigoa
          where             survey = 'GOA'
          order by          regulatory_area_name asc, min_depth asc"

areadepth <- sqlQuery(channel_akfin, query) %>% rename_all(tolower)
areadepth

# need to remove area estimate dups
areadepth <- areadepth %>%
  select(-area) %>%
  distinct() %>%
  arrange(regulatory_area_name, min_depth)

# GOA bottom trawl survey (BTS) biomass ----

query <- "select   survey, year, summary_area_depth, species_code,
                   area_biomass as biomass, biomass_var as var,
                   haul_count, catch_count
          from     afsc.race_biomassinpfcdepthaigoa
          where    species_code in ('30020') and
                   survey = 'GOA'"

bts <- sqlQuery(channel_akfin, query) %>%
  write_csv(paste0(raw_path, "/goabts_sst_raw_", YEAR, ".csv"))

bts <- bts %>%
  rename_all(tolower) %>%
  left_join(areadepth) %>%
  rename(area_biomass = biomass) %>%
  mutate(strata = case_when(regulatory_area_name == 'EASTERN GOA' & max_depth <= 500 ~ 'EGOA (0-500 m)',
                            regulatory_area_name == 'EASTERN GOA' & min_depth >= 501 & max_depth <= 700 ~ 'EGOA (501-700 m)',
                            regulatory_area_name == 'EASTERN GOA' & min_depth >= 701 ~ 'EGOA (701-1000 m)',
                            regulatory_area_name == 'CENTRAL GOA' & max_depth <= 500 ~ 'CGOA (0-500 m)',
                            regulatory_area_name == 'CENTRAL GOA' & min_depth >= 501 & max_depth <= 700 ~ 'CGOA (501-700 m)',
                            regulatory_area_name == 'CENTRAL GOA' & min_depth >= 701 ~ 'CGOA (701-1000 m)',
                            regulatory_area_name == 'WESTERN GOA' & max_depth <= 500 ~ 'WGOA (0-500 m)',
                            regulatory_area_name == 'WESTERN GOA' & min_depth >= 501 & max_depth <= 700 ~ 'WGOA (501-700 m)',
                            regulatory_area_name == 'WESTERN GOA' & min_depth >= 701 ~ 'WGOA (701-1000 m)'
                            )) %>%
  arrange(strata, year)

# 2020 status quo: assume strata that have a biomass estimate but no variance
# (i.e. there was a catch = 1) has CV=0.1.
bts %>% filter(area_biomass > 0 & is.na(var))
bts <- bts %>% mutate(var = ifelse(area_biomass > 0 & is.na(var), (0.1 * area_biomass)^2, var))

biomass_dat <- bts %>%
  group_by(strata, year) %>%
  dplyr::summarise(biomass = sum(area_biomass, na.rm = TRUE),
                   cv = sqrt(sum(var, na.rm = TRUE)) / sum(area_biomass, na.rm = TRUE)) %>%
 full_join(expand.grid(strata = unique(bts$strata),
                        year = unique(bts$year))) %>%
  arrange(strata, year) %>%
  write_csv(paste0(dat_path, "/goa_sst_biomass_", YEAR, ".csv"))

biomass_dat %>% pivot_wider(id_cols = year, names_from = strata, values_from = biomass)

# LLS Relative Population Weights ----

query <- "select   species_code, year, council_management_area, rpw, rpw_var
         from      afsc.lls_fmp_subarea_all_strata
         where     species_code = '30020' and
                   council_management_area in ('Central Gulf of Alaska', 'Eastern Gulf of Alaska', 'Western Gulf of Alaska') and
                   country in ('United States')
         order by  year asc"

lls <- sqlQuery(channel_akfin, query) %>%
  rename_all(tolower) %>%
  write_csv(paste0(raw_path, "/goalls_sst_raw_", YEAR, ".csv"))

cpue_dat <- lls %>%
  mutate(strata = case_when(council_management_area == 'Central Gulf of Alaska' ~ 'CGOA',
                            council_management_area == 'Eastern Gulf of Alaska' ~ 'EGOA',
                            council_management_area == 'Western Gulf of Alaska' ~ 'WGOA')) %>%
  group_by(strata, year) %>%
  dplyr::summarise(cpue = sum(rpw, na.rm = TRUE),
                   cv = sqrt(sum(rpw_var, na.rm = TRUE)) / cpue) %>%
  write_csv(paste0(dat_path, "/goa_sst_rpw_", YEAR, ".csv"))
