#Statistical analysis of kelp density data
#Sabrina N. Grant
#April 22, 2026

rm(list=ls())
librarian::shelf(tidyverse, here, glmmTMB, emmeans, DHARMa)

##require:
#kelp_swath_counts_CC.csv

################################################################################
#set directories and load data
basedir2 <- here::here("data_files","monitoring_data","processed")
output2 <- here::here("Macro_Manuscript")

kelp_swath_processed <- read.csv(file.path(basedir2, "kelp_swath_counts_CC.csv")) %>%
  janitor::clean_names()

################################################################################
#clean and prepare data
#***same pipeline as visualization script***

kelp_swath_clean1 <- kelp_swath_processed %>%
  dplyr::select('year', 'mhw', 'site', 'affiliated_mpa', 'zone', 'transect',
                'pterygophora_californica', 'eisenia_arborea', 'macrocystis_pyrifera',
                'macrocystis_stipes', 'nereocystis_luetkeana', 'laminaria_setchellii') %>%
  mutate(macrocystis_stipes = replace_na(macrocystis_stipes, 0))

kelp_swath_clean2 <- kelp_swath_clean1 %>%
  pivot_longer(pterygophora_californica:laminaria_setchellii,
               names_to = "species", values_to = "total_count") %>%
  group_by(year, mhw, site, affiliated_mpa, zone, transect, species) %>%
  summarise(raw_density = total_count/60, .groups = "drop")

kelp_swath_clean3 <- kelp_swath_clean2 %>%
  group_by(year, mhw, site, affiliated_mpa, zone, species) %>%
  summarise(mean_density_transect = mean(raw_density), .groups = "drop")

kelp_swath_clean_4 <- kelp_swath_clean3 %>%
  group_by(year, mhw, site, affiliated_mpa, species) %>%
  summarise(mean_density = mean(mean_density_transect), .groups = "drop")

kelp_swath_clean5 <- kelp_swath_clean_4 %>%
  mutate(region = ifelse(site %in% c("BLUEFISH_DC", "BLUEFISH_UC", "BUTTERFLY_DC",
                                     "BUTTERFLY_UC", "MONASTERY_UC", "LONE_TREE",
                                     "MONASTERY_DC", "PESCADERO_DC", "PESCADERO_UC",
                                     "STILLWATER_DC", "STILLWATER_UC", "WESTON_DC",
                                     "WESTON_UC"), "Carmel", "Monterey"))

kelp_swath_clean6 <- kelp_swath_clean5 %>%
  mutate(kelp_type = case_when(
    species %in% c("macrocystis_pyrifera", "macrocystis_stipes") ~ "macro",
    species == "nereocystis_luetkeana" ~ "nereo",
    species %in% c("eisenia_arborea", "laminaria_setchellii",
                   "pterygophora_californica") ~ "understory"
  ))

################################################################################
#prepare data for models
#set mhw as factor with before as reference level
#no transformation needed — glmmTMB with Tweedie handles zeros and skew directly

kelp_model_data <- kelp_swath_clean6 %>%
  mutate(mhw = factor(mhw, levels = c("before", "during", "after")),
         site = as.factor(site),
         region = as.factor(region),
         kelp_type = as.factor(kelp_type))

#understory data — averaged across understory species before modeling
understory_model_data <- kelp_model_data %>%
  filter(kelp_type == "understory") %>%
  group_by(year, mhw, site, region) %>%
  summarise(mean_density = mean(mean_density), .groups = "drop") %>%
  mutate(site = as.factor(site))

################################################################################
#fit glmmTMB models with Tweedie distribution
#fixed effect: mhw (before/during/after)
#random effect: site (repeated sampling at same locations)
#Tweedie family handles continuous non-negative data with excess zeros

################################################################################
#macro models - individuals (macrocystis_pyrifera)

macro_ind_carmel <- glmmTMB(mean_density ~ mhw + (1|site),
                            data = kelp_model_data %>%
                              filter(region == "Carmel",
                                     species == "macrocystis_pyrifera"),
                            family = tweedie(link = "log"),
                            REML = TRUE)

macro_ind_monterey <- glmmTMB(mean_density ~ mhw + (1|site),
                              data = kelp_model_data %>%
                                filter(region == "Monterey",
                                       species == "macrocystis_pyrifera"),
                              family = tweedie(link = "log"),
                              REML = TRUE)

################################################################################
#macro models - stipes (macrocystis_stipes)

macro_stipes_carmel <- glmmTMB(mean_density ~ mhw + (1|site),
                               data = kelp_model_data %>%
                                 filter(region == "Carmel",
                                        species == "macrocystis_stipes"),
                               family = tweedie(link = "log"),
                               REML = TRUE)

macro_stipes_monterey <- glmmTMB(mean_density ~ mhw + (1|site),
                                 data = kelp_model_data %>%
                                   filter(region == "Monterey",
                                          species == "macrocystis_stipes"),
                                 family = tweedie(link = "log"),
                                 REML = TRUE)

################################################################################
#nereo models

nereo_carmel <- glmmTMB(mean_density ~ mhw + (1|site),
                        data = kelp_model_data %>%
                          filter(region == "Carmel",
                                 species == "nereocystis_luetkeana"),
                        family = tweedie(link = "log"),
                        REML = TRUE)

nereo_monterey <- glmmTMB(mean_density ~ mhw + (1|site),
                          data = kelp_model_data %>%
                            filter(region == "Monterey",
                                   species == "nereocystis_luetkeana"),
                          family = tweedie(link = "log"),
                          REML = TRUE)

################################################################################
#understory models

understory_carmel <- glmmTMB(mean_density ~ mhw + (1|site),
                             data = understory_model_data %>%
                               filter(region == "Carmel"),
                             family = tweedie(link = "log"),
                             REML = TRUE)

understory_monterey <- glmmTMB(mean_density ~ mhw + (1|site),
                               data = understory_model_data %>%
                                 filter(region == "Monterey"),
                               family = tweedie(link = "log"),
                               REML = TRUE)

################################################################################
#check model assumptions using DHARMa simulated residuals
#DHARMa simulates residuals from the fitted model and compares them to observed
#values — this works correctly for non-normal distributions like Tweedie

check_assumptions_dharma <- function(model, model_name) {
  sim_res <- simulateResiduals(fittedModel = model, n = 1000)
  plot(sim_res, main = model_name)
  #also test for overdispersion and zero inflation explicitly
  cat("\n---", model_name, "---\n")
  print(testDispersion(sim_res))
  print(testZeroInflation(sim_res))
}

check_assumptions_dharma(macro_ind_carmel,    "Macro individuals - Carmel")
check_assumptions_dharma(macro_ind_monterey,  "Macro individuals - Monterey")
check_assumptions_dharma(macro_stipes_carmel,   "Macro stipes - Carmel")
check_assumptions_dharma(macro_stipes_monterey, "Macro stipes - Monterey")
check_assumptions_dharma(nereo_carmel,        "Nereo - Carmel")
check_assumptions_dharma(nereo_monterey,      "Nereo - Monterey")
check_assumptions_dharma(understory_carmel,   "Understory - Carmel")
check_assumptions_dharma(understory_monterey, "Understory - Monterey")

################################################################################
#extract results
#using emmeans for pairwise comparisons between MHW periods
#type = "response" back-transforms estimates to the original density scale

get_results <- function(model, model_name) {
  #overall likelihood ratio test for mhw effect
  anova_res <- car::Anova(model, type = "III")
  
  #pairwise comparisons between before/during/after
  emm <- emmeans(model, ~ mhw, type = "response")
  pairs_res <- pairs(emm, adjust = "tukey") %>%
    as.data.frame() %>%
    mutate(model = model_name)
  
  #estimated marginal means on response scale
  emm_df <- as.data.frame(emm) %>%
    mutate(model = model_name)
  
  list(anova = anova_res,
       emmeans = emm_df,
       pairwise = pairs_res)
}

results_macro_ind_carmel      <- get_results(macro_ind_carmel,      "Macro individuals - Carmel")
results_macro_ind_monterey    <- get_results(macro_ind_monterey,    "Macro individuals - Monterey")
results_macro_stipes_carmel   <- get_results(macro_stipes_carmel,   "Macro stipes - Carmel")
results_macro_stipes_monterey <- get_results(macro_stipes_monterey, "Macro stipes - Monterey")
results_nereo_carmel          <- get_results(nereo_carmel,          "Nereo - Carmel")
results_nereo_monterey        <- get_results(nereo_monterey,        "Nereo - Monterey")
results_understory_carmel     <- get_results(understory_carmel,     "Understory - Carmel")
results_understory_monterey   <- get_results(understory_monterey,   "Understory - Monterey")

################################################################################
#compile pairwise comparison results into one summary table

pairwise_summary <- bind_rows(
  results_macro_ind_carmel$pairwise,
  results_macro_ind_monterey$pairwise,
  results_macro_stipes_carmel$pairwise,
  results_macro_stipes_monterey$pairwise,
  results_nereo_carmel$pairwise,
  results_nereo_monterey$pairwise,
  results_understory_carmel$pairwise,
  results_understory_monterey$pairwise
) %>%
  dplyr::select(model, contrast, ratio, SE, df, z.ratio, p.value) %>%
  #note: emmeans returns ratio (not estimate) when type = "response" for log link
  mutate(p.value = round(p.value, 4),
         ratio = round(ratio, 3),
         SE = round(SE, 3),
         significance = case_when(
           p.value < 0.001 ~ "***",
           p.value < 0.01  ~ "**",
           p.value < 0.05  ~ "*",
           TRUE            ~ "ns"
         ))

print(pairwise_summary)

################################################################################
#compile emmeans summary table

emmeans_summary <- bind_rows(
  results_macro_ind_carmel$emmeans,
  results_macro_ind_monterey$emmeans,
  results_macro_stipes_carmel$emmeans,
  results_macro_stipes_monterey$emmeans,
  results_nereo_carmel$emmeans,
  results_nereo_monterey$emmeans,
  results_understory_carmel$emmeans,
  results_understory_monterey$emmeans
)

print(emmeans_summary)

################################################################################
#last write April 22, 2026