#visualizations from MPLA data
#Sabrina N. Grant
#April 22, 2026

rm(list=ls())
librarian::shelf(tidyverse, here, ggplot2, patchwork, mgcv)

##require:
#kelp_swath_counts_CC.csv

################################################################################
#set directories and load data
basedir2 <- here::here("data_files","monitoring_data","processed")           
output1 <- "/Users/sabrinagrant/Desktop/code_repositories/macro_manuscript/macro_manuscript/figures/MPLA"

kelp_swath_processed <- read.csv(file.path(basedir2, "kelp_swath_counts_CC.csv")) %>%
  janitor::clean_names()

kelp_upc_processed <- read.csv(file.path(basedir2, "kelp_upc_cov_CC.csv")) %>%
  janitor::clean_names()


################################################################################
#clean kelp swath

kelp_swath_clean1 <- kelp_swath_processed %>%
  #selecting only for species I'm interested in (macro, nereo, ptery, lsetch, eisarb)
  dplyr::select('year',
                'mhw',
                'site', 
                'affiliated_mpa', 
                'zone', 
                'transect', 
                'pterygophora_californica', 
                'eisenia_arborea', 
                'macrocystis_pyrifera',
                'macrocystis_stipes',
                'nereocystis_luetkeana', 
                'laminaria_setchellii') %>%
  #converting macro stipe NAs to zeroes 
  mutate(macrocystis_stipes = replace_na(macrocystis_stipes, 0))
  

#convert total count to density/m2
#***operating on the assumption that transects are 30 x 2 m2 therefore total surface area is 60m2 
kelp_swath_clean2 <- kelp_swath_clean1 %>%
  #first need to pivot back to long format
  pivot_longer(pterygophora_californica:laminaria_setchellii, names_to = "species", values_to = "total_count") %>%
  #then dividing total count by 60 to get individuals / m2
  group_by(year, mhw, site, affiliated_mpa, zone, transect, species) %>%
  summarise(raw_density = total_count/60)
 
#averaging density across transects
kelp_swath_clean3 <- kelp_swath_clean2 %>%
  group_by(year, mhw, site, affiliated_mpa, zone, species) %>%
  summarise(mean_density_transect = mean(raw_density))

#averaging density across zones 
kelp_swath_clean_4 <- kelp_swath_clean3 %>%
  group_by(year, mhw, site, affiliated_mpa, species) %>%
  summarise(mean_density = mean(mean_density_transect), .groups = "drop")

#organizing sites by Carmel vs Monterey
kelp_swath_clean5 <- kelp_swath_clean_4 %>%
  mutate(region = ifelse(site %in% c("BLUEFISH_DC",
                                     "BLUEFISH_UC",
                                     "BUTTERFLY_DC",   
                                     "BUTTERFLY_UC",
                                     "MONASTERY_UC",
                                     "LONE_TREE",
                                     "MONASTERY_DC",
                                     "PESCADERO_DC",
                                     "PESCADERO_UC",
                                     "STILLWATER_DC",
                                     "STILLWATER_UC",
                                     "WESTON_DC",
                                     "WESTON_UC"), "Carmel", "Monterey"))

#creating a column that has species-type specification
kelp_swath_clean6 <- kelp_swath_clean5 %>%
  mutate(kelp_type = case_when(
    species %in% c("macrocystis_pyrifera", "macrocystis_stipes") ~ "macro",
    species == "nereocystis_luetkeana" ~ "nereo",
    species %in% c("eisenia_arborea", "laminaria_setchellii", "pterygophora_californica") ~ "understory"
  ))

################################################################################

#visualizations of raw density

################################################################################

#prepare data
#overall (combined regions) - macro individuals
kelp_swath_overall_ind <- kelp_swath_clean6 %>%
  filter(species != "macrocystis_stipes") %>%
  group_by(year, mhw, kelp_type) %>%
  summarise(mean_density = mean(mean_density), .groups = "drop")

#overall (combined regions) - macro stipes
kelp_swath_overall_stipes <- kelp_swath_clean6 %>%
  filter(species != "macrocystis_pyrifera") %>%
  group_by(year, mhw, kelp_type) %>%
  summarise(mean_density = mean(mean_density), .groups = "drop")

#by region - macro individuals
kelp_swath_region_ind <- kelp_swath_clean6 %>%
  filter(species != "macrocystis_stipes") %>%
  group_by(year, mhw, region, kelp_type) %>%
  summarise(mean_density = mean(mean_density), .groups = "drop")

#by region - macro stipes
kelp_swath_region_stipes <- kelp_swath_clean6 %>%
  filter(species != "macrocystis_pyrifera") %>%
  group_by(year, mhw, region, kelp_type) %>%
  summarise(mean_density = mean(mean_density), .groups = "drop")

################################################################################
# define shared elements

mhw_rect <- annotate("rect", xmin = 2014, xmax = 2016,
                     ymin = -Inf, ymax = Inf,
                     alpha = 0.2, fill = "tomato")

ssw_line <- geom_vline(xintercept = 2013,
                             linetype = "dotted",
                             color = "darkgrey",
                             linewidth = 0.8)

kelp_colors <- c("macro" = "#2a9d8f",
                 "nereo" = "#e9c46a",
                 "understory" = "#264653")

kelp_labels <- c("macro" = expression(italic("Macrocystis pyrifera")),
                 "nereo" = expression(italic("Nereocystis luetkeana")),
                 "understory" = "Understory species")

shared_theme <- theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "white"),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 11, face = "bold"),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 9),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.position = "bottom")

################################################################################
#Fig. 1: Mean density of kelp species in Central California over time (2007-2022). A: Macrocystis pyrifera as # of individuals/m-2, B: Macrocystis pyrifera as # of stipes/m-2

a1 <- ggplot(data = kelp_swath_overall_ind,
             aes(x = year, y = mean_density, color = kelp_type)) +
  mhw_rect +
  ssw_line+
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_color_manual(values = kelp_colors, labels = kelp_labels) +
  scale_x_continuous(breaks = seq(2007, 2024, by = 3)) +
  labs(x = "Year",
       y = expression("Mean density (individuals m"^-2*")"),
       title = "A") +
  shared_theme

b1 <- ggplot(data = kelp_swath_overall_stipes,
             aes(x = year, y = mean_density, color = kelp_type)) +
  mhw_rect +
  ssw_line+
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_color_manual(values = kelp_colors, labels = kelp_labels) +
  scale_x_continuous(breaks = seq(2007, 2024, by = 3)) +
  labs(x = "Year",
       y = expression("Mean density (stipes m"^-2*")"),
       title = "B") +
  shared_theme

fig1 <- a1 | b1 +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

################################################################################
#Fig. 2: Mean density of kelp species in Carmel vs. Monterey over time (2007-2022). A: Macrocystis pyrifera as # of individuals/m-2, B: Macrocystis pyrifera as # of stipes/m-2

a2 <- ggplot(data = kelp_swath_region_ind %>% filter(region == "Carmel"),
             aes(x = year, y = mean_density, color = kelp_type)) +
  mhw_rect +
  ssw_line+
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_color_manual(values = kelp_colors, labels = kelp_labels) +
  scale_x_continuous(breaks = seq(2007, 2024, by = 3)) +
  labs(x = "Year",
       y = expression("Mean density (individuals m"^-2*")"),
       title = "A. Carmel ") +
  shared_theme

b2 <- ggplot(data = kelp_swath_region_ind %>% filter(region == "Monterey"),
             aes(x = year, y = mean_density, color = kelp_type)) +
  mhw_rect +
  ssw_line+
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_color_manual(values = kelp_colors, labels = kelp_labels) +
  scale_x_continuous(breaks = seq(2007, 2024, by = 3)) +
  labs(x = "Year",
       y = "",
       title = "B.  Monterey") +
  shared_theme

c2 <- ggplot(data = kelp_swath_region_stipes %>% filter(region == "Carmel"),
             aes(x = year, y = mean_density, color = kelp_type)) +
  mhw_rect +
  ssw_line+
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_color_manual(values = kelp_colors, labels = kelp_labels) +
  scale_x_continuous(breaks = seq(2007, 2024, by = 3)) +
  labs(x = "Year",
       y = expression("Mean density (stipes m"^-2*")"),
       title = "C. Carmel ") +
  shared_theme

d2 <- ggplot(data = kelp_swath_region_stipes %>% filter(region == "Monterey"),
             aes(x = year, y = mean_density, color = kelp_type)) +
  mhw_rect +
  ssw_line+
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_color_manual(values = kelp_colors, labels = kelp_labels) +
  scale_x_continuous(breaks = seq(2007, 2024, by = 3)) +
  labs(x = "Year",
       y = "",
       title = "D.  Monterey") +
  shared_theme

fig2 <- (a2 | b2) / (c2 | d2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

################################################################################
# print figures

fig1
fig2

################################################################################

#Now getting averages at the site level 
# overall (combined regions) - individuals, site level
kelp_swath_overall_ind_site <- kelp_swath_clean6 %>%
  filter(species != "macrocystis_stipes") %>%
  group_by(year, mhw, site, kelp_type) %>%
  summarise(mean_density = mean(mean_density), .groups = "drop")

# overall (combined regions) - stipes, site level
kelp_swath_overall_stipes_site <- kelp_swath_clean6 %>%
  filter(species != "macrocystis_pyrifera") %>%
  group_by(year, mhw, site, kelp_type) %>%
  summarise(mean_density = mean(mean_density), .groups = "drop")

# by region - individuals, site level
kelp_swath_region_ind_site <- kelp_swath_clean6 %>%
  filter(species != "macrocystis_stipes") %>%
  group_by(year, mhw, region, site, kelp_type) %>%
  summarise(mean_density = mean(mean_density), .groups = "drop")

# by region - stipes, site level
kelp_swath_region_stipes_site <- kelp_swath_clean6 %>%
  filter(species != "macrocystis_pyrifera") %>%
  group_by(year, mhw, region, site, kelp_type) %>%
  summarise(mean_density = mean(mean_density), .groups = "drop")

################################################################################
#Fig. 3: Mean density of kelp species in Central California fitted with a Generalize Additive Model over time (2007-2022). A: Macrocystis pyrifera as # of individuals/m-2, B: Macrocystis pyrifera as # of stipes/m-2

a3 <- ggplot(data = kelp_swath_overall_ind_site, 
       aes(x = year, y = mean_density, color = kelp_type, fill = kelp_type)) +
  mhw_rect +
  ssw_line+
  geom_point(size = 1.5, alpha = 0.4) +  # individual site points, slightly transparent
  geom_smooth(method = "gam",
              formula = y ~ s(x, bs = "cr", sp = 0.05),  # cubic regression spline, lambda = 0.05
              se = TRUE,                                   # 95% CI shading
              linewidth = 0.8) + 
  scale_color_manual(values = kelp_colors, labels = kelp_labels) +
  scale_fill_manual(values = kelp_colors, labels = kelp_labels) +
  scale_x_continuous(breaks = seq(2007, 2024, by = 3)) +
  labs(x = "Year",
       y = expression("Mean density (individuals m"^-2*")"),
       title = "A. Number of Individuals") +
  shared_theme

b3 <- ggplot(data = kelp_swath_overall_stipes_site,
             aes(x = year, y = mean_density, color = kelp_type, fill = kelp_type)) +
  mhw_rect +
  ssw_line+
  geom_point(size = 1.5, alpha = 0.4) +  # individual site points, slightly transparent
  geom_smooth(method = "gam",
              formula = y ~ s(x, bs = "cr", sp = 0.05),  # cubic regression spline, lambda = 0.05
              se = TRUE,                                   # 95% CI shading
              linewidth = 0.8) + 
  scale_color_manual(values = kelp_colors, labels = kelp_labels) +
  scale_fill_manual(values = kelp_colors, labels = kelp_labels) +
  scale_x_continuous(breaks = seq(2007, 2024, by = 3)) +
  labs(x = "Year",
       y = expression("Mean density (stipes m"^-2*")"),
       title = "B. Number of Stipes") +
  shared_theme

fig3 <- a3 | b3 +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

################################################################################
#Fig. 4: Mean density of kelp species across regions fitted with a Generalize Additive Model over time (2007-2022). A: Carmel Macrocystis pyrifera as # of individuals/m-2, B: Monterey Macrocystis pyrifera as # of individuals/m-2, C: Carmel macpyr as #of stipes, D: Monterey Macppyr as # of stipes

a4 <- ggplot(data = kelp_swath_region_ind_site %>% filter(region == "Carmel"), 
          aes(x = year, y = mean_density, color = kelp_type, fill = kelp_type)) +
  mhw_rect +
  ssw_line+
  geom_point(size = 1.5, alpha = 0.4) +  # individual site points, slightly transparent
  geom_smooth(method = "gam",
              formula = y ~ s(x, bs = "cr", sp = 0.05),  # cubic regression spline, lambda = 0.05
              se = TRUE,                                   # 95% CI shading
              linewidth = 0.8) + 
  scale_color_manual(values = kelp_colors, labels = kelp_labels) +
  scale_fill_manual(values = kelp_colors, labels = kelp_labels) +
  scale_x_continuous(breaks = seq(2007, 2024, by = 3)) +
  labs(x = "Year",
       y = expression("Mean density (individuals m"^-2*")"),
       title = "A. Carmel") +
  shared_theme

b4 <- ggplot(data = kelp_swath_region_ind_site %>% filter(region == "Monterey"), 
             aes(x = year, y = mean_density, color = kelp_type, fill = kelp_type)) +
  mhw_rect +
  ssw_line+
  geom_point(size = 1.5, alpha = 0.4) +  # individual site points, slightly transparent
  geom_smooth(method = "gam",
              formula = y ~ s(x, bs = "cr", sp = 0.05),  # cubic regression spline, lambda = 0.05
              se = TRUE,                                   # 95% CI shading
              linewidth = 0.8) + 
  scale_color_manual(values = kelp_colors, labels = kelp_labels) +
  scale_fill_manual(values = kelp_colors, labels = kelp_labels) +
  scale_x_continuous(breaks = seq(2007, 2024, by = 3)) +
  labs(x = "Year",
       y = expression("Mean density (individuals m"^-2*")"),
       title = "B. Monterey") +
  shared_theme

c4 <- ggplot(data = kelp_swath_region_stipes_site %>% filter(region == "Carmel"), 
             aes(x = year, y = mean_density, color = kelp_type, fill = kelp_type)) +
  mhw_rect +
  ssw_line+
  geom_point(size = 1.5, alpha = 0.4) +  # individual site points, slightly transparent
  geom_smooth(method = "gam",
              formula = y ~ s(x, bs = "cr", sp = 0.05),  # cubic regression spline, lambda = 0.05
              se = TRUE,                                   # 95% CI shading
              linewidth = 0.8) + 
  scale_color_manual(values = kelp_colors, labels = kelp_labels) +
  scale_fill_manual(values = kelp_colors, labels = kelp_labels) +
  scale_x_continuous(breaks = seq(2007, 2024, by = 3)) +
  labs(x = "Year",
       y = expression("Mean density (stipes m"^-2*")"),
       title = "C") +
  shared_theme

d4 <- ggplot(data = kelp_swath_region_stipes_site %>% filter(region == "Monterey"), 
             aes(x = year, y = mean_density, color = kelp_type, fill = kelp_type)) +
  mhw_rect +
  ssw_line+
  geom_point(size = 1.5, alpha = 0.4) +  # individual site points, slightly transparent
  geom_smooth(method = "gam",
              formula = y ~ s(x, bs = "cr", sp = 0.05),  # cubic regression spline, lambda = 0.05
              se = TRUE,                                   # 95% CI shading
              linewidth = 0.8) + 
  scale_color_manual(values = kelp_colors, labels = kelp_labels) +
  scale_fill_manual(values = kelp_colors, labels = kelp_labels) +
  scale_x_continuous(breaks = seq(2007, 2024, by = 3)) +
  labs(x = "Year",
       y = expression("Mean density (stipes m"^-2*")"),
       title = "D") +
  shared_theme

fig4 <- (a4 | b4) / (c4 | d4) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

################################################################################
# print figures

fig3
fig4

################################################################################

#visualizations of proportional abundance

################################################################################

#preparing dfs

kelp_swath_overall_ind_prop <- kelp_swath_overall_ind %>%
  group_by(year, mhw) %>%
  mutate(kelp_type = factor(kelp_type,
                            levels = c("understory","nereo", "macro"))) %>%
  mutate(prop_abundance = mean_density / sum(mean_density, na.rm = TRUE)) %>%
  ungroup()


kelp_swath_overall_stipes_prop <- kelp_swath_overall_stipes %>%
  group_by(year, mhw) %>%
  mutate(kelp_type = factor(kelp_type,
                            levels = c("understory","nereo", "macro"))) %>%
  mutate(prop_abundance = mean_density / sum(mean_density, na.rm = TRUE)) %>%
  ungroup()

kelp_swath_region_ind_prop <- kelp_swath_region_ind %>%
  group_by(year, mhw, region) %>%
  mutate(kelp_type = factor(kelp_type,
                            levels = c("understory","nereo", "macro"))) %>%
  mutate(prop_abundance = mean_density / sum(mean_density, na.rm = TRUE)) %>%
  ungroup()

kelp_swath_region_stipes_prop <- kelp_swath_region_stipes %>%
  group_by(year, mhw, region) %>%
  mutate(kelp_type = factor(kelp_type,
                            levels = c("understory","nereo", "macro"))) %>%
  mutate(prop_abundance = mean_density / sum(mean_density, na.rm = TRUE)) %>%
  ungroup()

################################################################################
#Fig 5: Proportional abundance of kelp type in Central California over time
#A: macro # of individuals
#B: macro # of stipes

a5 <- ggplot(kelp_swath_overall_ind_prop,
       aes(x = year, y = prop_abundance, fill = kelp_type)) +
  geom_col(width = 0.8, color = "white") +
  scale_fill_manual(values = kelp_colors, labels = kelp_labels) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = c(0, 0)) +
  labs(x = "Year",
       y = "Proportional abundance",
       title = "A. (number of individuals)", 
       fill = NULL) +
  shared_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

b5 <- ggplot(kelp_swath_overall_stipes_prop,
             aes(x = year, y = prop_abundance, fill = kelp_type)) +
  geom_col(width = 0.8, color = "white") +
  scale_fill_manual(values = kelp_colors, labels = kelp_labels) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = c(0, 0)) +
  labs(x = "Year",
       y = "Proportional abundance",
       title = "B. (number of stipes)", 
       fill = NULL) +
  shared_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

fig5 <- a5 + b5

################################################################################
#Fig 6: Proportional abundance of kelp type in Central California across regions over time
#A: Carmel macro # of individuals
#B: Monterey macro # of individuals
#C: Carmel macro # of stipes 
#D: Monterey macro # of stipes

a6 <- ggplot(data = kelp_swath_region_ind_prop %>% filter(region == "Carmel"), 
             aes(x = year, y = prop_abundance, fill = kelp_type)) +
  geom_col(width = 0.8, color = "white") +
  scale_fill_manual(values = kelp_colors, labels = kelp_labels) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = c(0, 0)) +
  labs(x = "Year",
       y = "Proportional abundance",
       title = "A. Carmel (number of individuals)",
       fill = NULL) +
  shared_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

b6 <- ggplot(data = kelp_swath_region_ind_prop %>% filter(region == "Monterey"), 
           aes(x = year, y = prop_abundance, fill = kelp_type)) +
  geom_col(width = 0.8, color = "white") +
  scale_fill_manual(values = kelp_colors, labels = kelp_labels) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = c(0, 0)) +
  labs(x = "Year",
       y = "Proportional abundance",
       title = "B. Monterey (number of individuals)",
       fill = NULL) +
  shared_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

c6 <- ggplot(data = kelp_swath_region_stipes_prop %>% filter(region == "Carmel"), 
             aes(x = year, y = prop_abundance, fill = kelp_type)) +
  geom_col(width = 0.8, color = "white") +
  scale_fill_manual(values = kelp_colors, labels = kelp_labels) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = c(0, 0)) +
  labs(x = "Year",
       y = "Proportional abundance",
       title = "C. (number of stipes)",
       fill = NULL) +
  shared_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

d6 <- ggplot(data = kelp_swath_region_stipes_prop %>% filter(region == "Monterey"), 
             aes(x = year, y = prop_abundance, fill = kelp_type)) +
  geom_col(width = 0.8, color = "white") +
  scale_fill_manual(values = kelp_colors, labels = kelp_labels) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = c(0, 0)) +
  labs(x = "Year",
       y = "Proportional abundance",
       title = "D. (number of stipes)",
       fill = NULL) +
  shared_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

fig6 <- (a6 | b6) / (c6 | d6) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

################################################################################
# print figures

fig5
fig6

################################################################################
#SAVING ALL FIGURES I WANT 

ggsave(file.path(output1, "Fig_density_#plants_GAM.png"),
       plot = a3, width = 18, height = 10, dpi = 300, bg = "white")

ggsave(file.path(output1, "Fig_density_#stipes_GAM.png"),
       plot = b3, width = 18, height = 10, dpi = 300, bg = "white")

ggsave(file.path(output1, "Fig_Carmel_#plants_GAM.png"),
       plot = a4, width = 18, height = 10, dpi = 300, bg = "white")

ggsave(file.path(output1, "Fig_Monterey_#plants_GAM.png"),
       plot = b4, width = 18, height = 10, dpi = 300, bg = "white")

ggsave(file.path(output1, "Fig_Carmel_#stipes_GAM.png"),
       plot = c4, width = 18, height = 10, dpi = 300, bg = "white")

ggsave(file.path(output1, "Fig_Monterey_#stipes_GAM.png"),
       plot = d4, width = 18, height = 10, dpi = 300, bg = "white")


ggsave(file.path(output1, "Fig_prop_abundance_#plants_GAM.png"),
       plot = a5, width = 18, height = 10, dpi = 300, bg = "white")

ggsave(file.path(output1, "Fig_prop_abundance_#stipes_GAM.png"),
       plot = b5, width = 18, height = 10, dpi = 300, bg = "white")

ggsave(file.path(output1, "Fig_Carmel_prop_abundance_#plants_GAM.png"),
       plot = a6, width = 18, height = 10, dpi = 300, bg = "white")

ggsave(file.path(output1, "Fig_Carmel_prop_abundance_#stipes_GAM.png"),
       plot = c6, width = 18, height = 10, dpi = 300, bg = "white")

ggsave(file.path(output1, "Fig_Monterey_prop_abundance_#plants_GAM.png"),
       plot = b6, width = 18, height = 10, dpi = 300, bg = "white")

ggsave(file.path(output1, "Fig_Monterey_prop_abundance_#stipes_GAM.png"),
       plot = d6, width = 18, height = 10, dpi = 300, bg = "white")


################################################################################
#last write April 22, 2026

