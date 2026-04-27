#MBA Community data processing
#Sabrina N. Grant 
#April 22, 2026

rm(list=ls())
librarian::shelf("tidyverse", "vegan", "ggplot2", "dplyr", "googledrive", "janitor", "here")

################################################################################
#load data and set directories
drive_deauth()
drive_auth(scopes = "https://www.googleapis.com/auth/drive")


file <- drive_find("kelp_recovery_data.rda",
                   shared_drive = "MBA_kelp_recovery_research")
drive_download(
  as_id(file$id),
  path = "kelp_recovery_data.rda",
  overwrite = TRUE
)

load("kelp_recovery_data.rda")
ls()

output <- here::here("data_files","mba_data","processed")

################################################################################
#processing quad data 

###averaging counts and % cover per transect

#helper function to get the most common value (mode)
get_mode <- function(x) {
  ux <- na.omit(unique(x))
  ux[which.max(tabulate(match(x, ux)))]
}

#averaging across quads
quad_transect_avgs <- quad_data %>%
  group_by(site, site_type, zone, survey_date, transect) %>%
  summarise(
    across(where(is.numeric) & !any_of("quadrat"), ~ mean(.x, na.rm = TRUE)),
    substrate = get_mode(substrate),
    .groups = "drop")

#now averaging per survey 
quad_avgs <- quad_transect_avgs %>%
  group_by(site, site_type, zone, survey_date) %>%
  #creating new column that has proportion of concealed urchins
  mutate(
    prop_purp_concealedm2 = if_else(
      purple_urchin_densitym2 > 0,
      purple_urchin_conceiledm2 / purple_urchin_densitym2,
      0),
    prop_red_concealedm2 = if_else(
      red_urchin_densitym2 > 0,
      red_urchin_conceiledm2 / red_urchin_densitym2,
      0)) %>%
  #dropping the old concealed density columns
  dplyr::select(-purple_urchin_conceiledm2, -red_urchin_conceiledm2) %>%
  #summarise by group
  summarise(
    across(where(is.numeric) & !any_of("transect"), ~mean(.x, na.rm = TRUE)),
    substrate = get_mode(substrate),
    .groups = "drop" )

################################################################################
#processing kelp data
#averaging kelp transects
kelp_avgs <- kelp_data %>%
  group_by(site, site_type, zone, survey_date) %>%
  summarise( 
    across(where(is.numeric) & !any_of("transect"), ~mean(.x, na.rm = TRUE)), 
    .groups = "drop") 

#combining quad and kelp data frames 
merged_data <- kelp_avgs %>%
  #keeping coordinates from kelp_avgs
  left_join(
    quad_avgs %>% dplyr::select(-latitude, -longitude),  # drop coords from quad_avgs before join
    by = c("site", "site_type", "zone", "survey_date")
  ) %>%
  mutate(year = as.numeric(format(survey_date, "%Y")))

#creating a new column that just has years for future analysis  
merged_data$year <- as.numeric(format(merged_data$survey_date, "%Y"))

#standardizing kelp counts to quad data 

#identify columns that contain "20m2"
cols_to_scale <- grep("20m2", names(merged_data), value = TRUE)

#divide those columns by 20 to get density / m2
merged_data[cols_to_scale] <- lapply(merged_data[cols_to_scale], function(x) x / 20)

#rename those columns by removing "20m2_"
names(merged_data)[names(merged_data) %in% cols_to_scale] <-
  gsub("20m2","", cols_to_scale)

#replaces NAs with zero
merged_data[sapply(merged_data, is.numeric)] <-
  lapply(merged_data[sapply(merged_data, is.numeric)],
         function(x) { x[is.na(x)] <- 0; x })  

################################################################################
#normalizing (z-score) data to compare UPC and Swath algae

norm_data <- merged_data %>%
  mutate(year = as.factor(year)) %>% # making year a factor instead of numeric
  mutate(site_type = as.factor(site_type)) %>%
  dplyr::select(site, site_type, zone, year, where(is.numeric)) %>%
  dplyr::select(-any_of(c("longitude", "latitude", "relief_cm", "risk_index"))) %>% # getting rid of columns I don't care about
  dplyr::mutate(across(where(is.numeric), scale)) # z-scoring data


################################################################################
#export


#Export
write.csv(norm_data,file.path(output, "normalized_kelp_quad_data_CC.csv"), row.names = FALSE)
write.csv(merged_data,file.path(output, "kelp_quad_data_CC.csv"), row.names = FALSE)

################################################################################
#last write April 23, 2026
