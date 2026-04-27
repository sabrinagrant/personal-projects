#MBA Community data visualizations/ analyses
#Sabrina N. Grant 
#April 23, 2026

rm(list=ls())
librarian::shelf(tidyverse, here, ggplot2, MASS, caret, patchwork, grid, RColorBrewer, reshape2, vegan, ggtext)


##require:
#normalized_kelp_quad_data_CC.csv
#kelp_quad_data_CC.csv

################################################################################
#set directories and load data
basedir2 <- here::here("data_files","mba_data","processed")           
output2 <- "/Users/sabrinagrant/Desktop/code_repositories/macro_manuscript/macro_manuscript/figures/MBA"



normalized_kelp_quad_data_CC <- read.csv(file.path(basedir2, "normalized_kelp_quad_data_CC.csv")) %>%
  janitor::clean_names()

kelp_quad_data_CC <- read.csv(file.path(basedir2, "kelp_quad_data_CC.csv")) %>%
  janitor::clean_names()

  
################################################################################
#renaming dfs
  
merged_data <- kelp_quad_data_CC
  
norm_data <- normalized_kelp_quad_data_CC  


#filtering for 2024 data 
  year_one_data <- merged_data %>% 
    filter(year == 2024) %>%
    dplyr::select(year, site, zone, site_type, macro_stipe_density, everything()) %>%
    mutate(site_type = as.factor(site_type))
  
  
#getting mean of numeric values across grouping categories
  year_one_data %>%
    group_by(site, site_type, zone, year) %>%
    summarise(
      across(where(is.numeric), \(x) mean(x, na.rm = TRUE)),
      .groups = "drop")
  
  #filtering for 2024 z-scored data 
  z_year_one_data <- norm_data %>% 
    filter(year == 2024) %>%
    dplyr::select(year, site, zone, site_type, macro_stipe_density, everything()) %>%
    mutate(site_type = as.factor(site_type)) 
  
  #getting mean of numeric values across grouping categories
  z_year_one_data_means <-  z_year_one_data %>%
    group_by(site, site_type, zone, year) %>%
    summarise(
      across(where(is.numeric), \(x) mean(x, na.rm = TRUE)),
      .groups = "drop")

################################################################################
#LDA on Z-Scored Data

lda_model_z <- lda(site_type ~ macro_stipe_density + 
                       density_ptecal + density_nerlue + 
                       density_lamset + density_eisarb + 
                       density_lamstump + 
                       density_macstump + 
                       lamr + 
                       macr +
                       macj + 
                       nerj + 
                       ptej + 
                       lsetj + 
                       cov_articulated_coralline + 
                       cov_crustose_coralline + 
                       cov_encrusting_red + 
                       cov_fleshy_red + 
                       cov_stephanocystis + 
                       cov_dictyoneurum_spp  + 
                       cov_desmarestia_spp + 
                       cov_lam_holdfast_live + 
                       cov_mac_holdfast_live, 
                     data = z_year_one_data)
  
#lda summary
lda_model_z
  
#getting LDA predictions
pred_z <- predict(lda_model_z)
pred_class_z <- pred_z$class

#getting % accuracy
confusionMatrix(pred_class_z, z_year_one_data$site_type)
#predicts at 81% accuracy which is the same as non z scored


################################################################################
#Fig. 7 A: LDA model B: Coefficients 

#A: LDA Model Fig
#getting LDA scores 
lda_scores_z <- as.data.frame(predict(lda_model_z)$x) 
lda_scores_z$site_type <- z_year_one_data$site_type 
  
#plotting LD1 vs LD2 
a7 <- ggplot(lda_scores_z, aes(x = LD1, y = LD2, color = site_type, shape = site_type)) + geom_point(size = 3.5, alpha = 0.9) + 
    stat_ellipse(level = 0.95, size = 1, alpha = 0.6) + 
    scale_color_brewer(palette = "Set1") + 
    theme_classic(base_size = 16) + 
    labs( title = "", subtitle = "", x = paste0("LD1 (", round(lda_model_z$svd[1]^2 / sum(lda_model_z$svd^2) * 100, 1), "%)"), y = paste0("LD2 (", round(lda_model_z$svd[2]^2 / sum(lda_model_z$svd^2) * 100, 1), "%)"), color = "Site Type", shape = "Site Type" ) + 
    theme( plot.title = element_text(face = "bold", size = 20, hjust = 0.5, margin = margin(b = 10)), plot.subtitle = element_text(size = 14, hjust = 0.5, margin = margin(b = 15)), axis.title = element_text(face = "bold", size = 15), axis.text = element_text(size = 13), legend.title = element_text(face = "bold", size = 14), legend.text = element_text(size = 12), panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "white", linewidth = 0.3) )

#B: Coefficients of LDA fig
# Get coefficients into tidy format
  coef_z <- as.data.frame(lda_model_z$scaling)
  coef_z$variable <- rownames(coef_z)
  coef_long <- pivot_longer(coef_z, cols = c(LD1, LD2), names_to = "LD", values_to = "coefficient")
  
  # clean up variable names for display
  coef_long <- coef_long %>%
    mutate(variable = str_replace_all(variable, "_", " ") %>%
             str_to_title())
  
  b7 <- ggplot(coef_long, aes(x = coefficient, 
                              y = reorder(variable, coefficient), 
                              fill = LD)) +
    geom_col(position = position_dodge(width = 0.7), 
             width = 0.6, 
             alpha = 0.85) +
    geom_vline(xintercept = 0, linewidth = 0.7, color = "black", linetype = "solid") +
    scale_fill_brewer(palette = "Set1") +
    theme_classic(base_size = 16) +
    labs(
      x = "Discriminant Coefficient",
      y = NULL,
      fill = "Linear\nDiscriminant"
    ) +
    theme(
      axis.title = element_text(face = "bold", size = 15),
      axis.text.y = element_text(size = 11),
      axis.text.x = element_text(size = 13),
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(size = 12),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3)
    )
  
  a7 
  b7
  fig7 <- a7 + b7 +
    plot_layout(ncol = 2, widths = c(1, 1.4)) +
    plot_annotation(
      tag_levels = "A",
      theme = theme(
        plot.tag = element_text(face = "bold", size = 18)
      )
    )
  fig7
  
################################################################################

#LDA model on 2025 data 

  #filtering for 2025 data 
  z_year_two_data <- norm_data %>% 
    filter(year == 2025) %>%
    dplyr::select(year, site, zone, site_type, macro_stipe_density, everything())
  
  #getting mean of numeric values across grouping categories
  z_year_two_data %>%
    group_by(site, site_type, zone, year) %>%
    summarise(
      across(where(is.numeric), \(x) mean(x, na.rm = TRUE)),
      .groups = "drop")
  
  #predicting 2025 patch types based on 2024 LDA model
  z_year_two_data <- z_year_two_data %>%
    mutate(site_type_pred = predict(lda_model_z, newdata = .)$class) %>%
    dplyr::select(year, site, zone, site_type, site_type_pred, everything())
  
  #setting up the 2025 df with the predictions column 
  z_year_one_data <- z_year_one_data %>%
    mutate(site_type_pred = predict(lda_model_z, newdata = .)$class) %>%
    dplyr::select(year, site, zone, site_type, site_type_pred, everything())
  
  ##REJOINING 2024 AND 2025 DFS
  mergedz <- z_year_one_data %>%
    distinct(site, zone, site_type, site_type_pred) %>%
    rename(site_type_pred_2024 = site_type_pred) %>%
    left_join(
      z_year_two_data %>%
        distinct(site, zone, site_type, site_type_pred) %>%
        rename(site_type_pred_2025 = site_type_pred),
      by = c("site", "zone", "site_type"))

################################################################################ 

#Fig. 8: Heatmap of patch-type transitions from 2024 to 2025 (on z-scored data)

  # build cross tab
  cross_tab <- table(
    `2024` = mergedz$site_type_pred_2024,
    `2025` = mergedz$site_type_pred_2025
  )
  
  # melt for ggplot
  cross_tab_melt <- reshape2::melt(cross_tab) %>%
    rename(Var1 = `2024`, Var2 = `2025`)
  
  # set factor order
  desired_order <- c("BAR", "INCIP", "FOR")
  
  cross_tab_melt <- cross_tab_melt %>%
    mutate(
      Var1 = factor(Var1, levels = desired_order),
      Var2 = factor(Var2, levels = desired_order)
    )
  
  # heatmap
  fig8 <- ggplot(cross_tab_melt, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile(color = "grey90", linewidth = 0.3) +
    geom_text(aes(label = value,
                  color = value > max(value) / 2),  # cleaner than ifelse in aes
              size = 5, fontface = "bold", show.legend = FALSE) +
    scale_color_manual(values = c("TRUE" = "white", "FALSE" = "black")) +
    scale_fill_gradientn(colors = brewer.pal(9, "Blues")) +  # 9 is max for Blues
    labs(
      x = "Predicted Patch Type — 2024",
      y = "Predicted Patch Type — 2025",
      title = "Transition of Predicted Patch Types from 2024 to 2025",
      fill = "Count"
    ) +
    theme_minimal(base_size = 15) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 13),
      axis.text.y = element_text(face = "bold", size = 13),
      axis.title = element_text(face = "bold", size = 14),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18, 
                                margin = margin(b = 10)),
      legend.title = element_text(face = "bold", size = 13),
      legend.text = element_text(size = 12),
      panel.grid = element_blank()
    )
  
  fig8
  
################################################################################
#Creating a new df that has LDA patch transitions for my SIMPER analysis (raw non-normalized data)

  predictions_df <- mergedz %>%
    left_join(merged_data, by = c("site", "zone", "site_type")) %>%
    dplyr::select(c(site,zone,site_type_pred_2024,site_type_pred_2025, macro_stipe_density, density_ptecal, density_nerlue, density_lamset, density_eisarb, density_lamstump, density_macstump, lamr , macr , macj , nerj , ptej , lsetj , cov_articulated_coralline , cov_crustose_coralline , cov_encrusting_red , cov_fleshy_red , cov_stephanocystis , cov_dictyoneurum_spp  , cov_desmarestia_spp , cov_lam_holdfast_live ,cov_mac_holdfast_live,year))
  
  patch_transitions <- predictions_df %>%
    mutate(
      #assigning the correct prediction for each year
      predicted_patch_type = case_when(
        year == 2024 ~ as.character(site_type_pred_2024),
        year == 2025 ~ as.character(site_type_pred_2025),
        TRUE ~ NA_character_),
      #creating the transition column
      patch_transition = paste(site_type_pred_2024, "→", site_type_pred_2025)
    ) %>%
    #dropping the old wide prediction columns, keep everything else (including species)
    dplyr::select(site, zone, year, predicted_patch_type, patch_transition, everything(), -site_type_pred_2024, -site_type_pred_2025)
  
  

  #creating groupings that I want for future analyses 
  patch_groupings <- patch_transitions %>%
    mutate(patches_groups = case_when(
      year == 2024 & predicted_patch_type == "BAR"   ~ "2024 Barrens",
      year == 2024 & predicted_patch_type == "INCIP" ~ "2024 Incipients",
      year == 2024 & predicted_patch_type == "FOR"   ~ "2024 Forests", 
      year == 2025 & patch_transition == "BAR → INCIP"  ~ "2025 Recovered Incipient Patches",
      year == 2025 & patch_transition == "BAR → FOR"    ~ "2025 BAR Recovered Forest Patches", #forests from barrens
      year == 2025 & patch_transition == "INCIP → FOR"  ~ "2025 INCIP Recovered Forest Patches", #forests from incipients
      year == 2025 & patch_transition == "BAR → BAR"    ~ "2025 Persistant Barren Patches",
      year == 2025 & patch_transition == "FOR → FOR"    ~ "2025 Persistant Forest Patches",
      year == 2025 & patch_transition == "INCIP → INCIP"~ "2025 Persistant Incipient Patches",
      year == 2025 & patch_transition == "INCIP → BAR"  ~ "2025 Regressed Barren Patches",
      year == 2025 & patch_transition == "FOR → BAR"    ~ "2025 Regressed Barren Patches",
      year == 2025 & patch_transition == "FOR → INCIP"  ~ "2025 Regressed Incipient Patches",
      TRUE ~ NA_character_))
  
################################################################################
#SIMPER Analysis on raw data. Separating it by % cover species vs count species 

#new df 
groups <- patch_groupings %>%
    dplyr::select(c(site, zone, year, patches_groups, macro_stipe_density:cov_mac_holdfast_live)) %>%
    filter(!is.na(patches_groups)) %>% #just keeping rows with my patch transitions a.k.a kelp recovery
    droplevels()
  
  #just percent cover species 
  raw_cover <- groups %>%
    dplyr::select(-c(macro_stipe_density:lsetj))
  cover_meta <- raw_cover %>%
    dplyr::select(c(site, zone, year, patches_groups))
  cover_comm<- raw_cover %>%
    dplyr::select(-c(site, zone, year, patches_groups))
  
  sim_cover <- simper(cover_comm, cover_meta$patches_groups, permutations = 999)
  
  summary(sim_cover)
  
  
  #just density species 
  raw_density <- groups %>%
    dplyr::select(c(site:lsetj))
  
  
  den_meta <- raw_density %>%
    dplyr::select(c(site, zone, year, patches_groups))
  den_comm<- raw_density %>%
    dplyr::select(-c(site, zone, year, patches_groups))
  den_comm_eps <- den_comm + 1e-6 #adding in a tiny value to account for zero inflation, dont want to get rid of zeros 
  d <- vegdist(den_comm_eps, method="bray")
  sim_den <- simper(den_comm_eps, den_meta$patches_groups, permutations = 999)
  
  summary(sim_den)

  
################################################################################
#SIMPER
 
#labeling species
   #species label lookups using html italic tags for ggtext
  den_labels <- c(
    "density_nerlue"      = "<i>N. luetkeana</i>",
    "ptej"                = "<i>P. californica</i> (juvenile)",
    "density_ptecal"      = "<i>P. californica</i>",
    "macro_stipe_density" = "<i>M. pyrifera</i> (# of stipes)",
    "macr"                = "<i>M. pyrifera</i> (recruit)",
    "lsetj"               = "<i>L. setchellii</i> (juvenile)",
    "density_lamstump"    = "Laminariales spp. (stump)",
    "macj"                = "<i>M. pyrifera</i> (juvenile)",
    "density_lamset"      = "<i>L. setchellii</i>",
    "nerj"                = "<i>N. luetkeana</i> (juvenile)",
    "lamr"                = "Laminariales spp. (recruit)")
  
  cov_labels <- c(
    "cov_crustose_coralline"    = "Crustose Coralline",
    "cov_encrusting_red"        = "Encrusting Red",
    "cov_fleshy_red"            = "Fleshy Red",
    "cov_articulated_coralline" = "Articulated Coralline",
    "cov_desmarestia_spp"       = "<i>Desmarestia</i> spp.",
    "cov_stephanocystis"        = "<i>S. osmundea</i>",
    "cov_dictyoneurum_spp"      = "<i>Dictyoneurum</i> spp.",
    "cov_mac_holdfast_live"     = "<i>M. pyrifera</i> holdfast",
    "cov_lam_holdfast_live"     = "Laminariales spp. holdfast")

################################################################################
#Fig. 9: SIMPER - Top contributors to patch transitions (density + cover)

  #bars colored by direction of change (increase = blue, decrease = red)
  #panel A = density spp, panel B = cover spp per transition
  plot_simper <- function(den_obj, 
                          cov_obj,
                          comparison,        
                          title_label,       
                          top_n = 10,
                          increase_color = "#2166AC",
                          decrease_color = "#D6604D") {
    
    #internal helper to tidy a simper object
    tidy_simper <- function(simper_obj, label_lookup) {
      
      if (!comparison %in% names(simper_obj)) {
        stop(paste0("Comparison '", comparison, "' not found in simper object.\n",
                    "Available comparisons:\n",
                    paste(names(simper_obj), collapse = "\n")))
      }
      
      as.data.frame(simper_obj[[comparison]]) %>%
        filter(!is.na(average), average > 0) %>%
        arrange(desc(average)) %>%
        mutate(
          percent        = 100 * average / sum(average),
          cumsum_percent = cumsum(percent),
          change         = ifelse(avb > ava, "Increase", "Decrease"),
          species        = dplyr::recode(species, !!!label_lookup),
          species        = factor(species, levels = rev(unique(species)))
        ) %>%
        slice_head(n = top_n)
    }
    
    den_df <- tidy_simper(den_obj, den_labels)
    cov_df <- tidy_simper(cov_obj, cov_labels)
    
    #internal helper to build a single panel
    build_panel <- function(df, panel_title, label_lookup) {
      
      #recode species to html labels for rendering
      df <- df %>%
          mutate(species_label = dplyr::recode(as.character(species), 
                                               !!!label_lookup,
                                               .default = as.character(species)),
                 species_label = factor(species_label, 
                                        levels = rev(unique(species_label))))
      
      ggplot(df, aes(x = species_label, y = percent, fill = change)) +
        geom_col(width = 0.7, alpha = 0.9, color = "black") +
        geom_line(aes(x = species_label, y = cumsum_percent, group = 1),
                  color = "black", linewidth = 0.9) +
        geom_point(aes(x = species_label, y = cumsum_percent),
                   color = "black", size = 2.5, shape = 21, fill = "white") +
        coord_flip() +
        scale_y_continuous(
          name = "Contribution (%)",
          sec.axis = sec_axis(~ ., name = "Cumulative (%)")
        ) +
        scale_fill_manual(values = c("Increase" = increase_color,
                                     "Decrease" = decrease_color),
                          name = "Direction") +
        labs(title = panel_title, x = NULL) +
        theme_minimal(base_size = 14) +
        theme(
          panel.grid.major.y = element_blank(),
          panel.grid.minor   = element_blank(),
          panel.grid.major.x = element_line(color = "grey80", linetype = "dotted"),
          axis.text.y        = ggtext::element_markdown(size = 12),  #renders html italic
          axis.text.x        = element_text(size = 12),
          axis.title.x       = element_text(face = "bold", size = 13),
          legend.title       = element_text(face = "bold", size = 13),
          legend.text        = element_text(size = 12),
          plot.title         = element_text(size = 14, face = "bold", hjust = 0.5)
        )
    }
    
    pA <- build_panel(den_df, paste0(title_label, " — Density Species"), den_labels)
    pB <- build_panel(cov_df, paste0(title_label, " — Percent Cover Species"), cov_labels)
    
    combined <- pA / pB +
      plot_annotation(
        tag_levels = "A",
        theme = theme(plot.tag = element_text(face = "bold", size = 16))
      )
    
    return(list(plot = combined, den_data = den_df, cov_data = cov_df))}
    
################################################################################  
  
  #BAR → FOR
  fig9a <- plot_simper(sim_den, sim_cover,
                       comparison  = "2024 Barrens_2025 BAR Recovered Forest Patches",
                       title_label = "BAR → FOR")
  
  #BAR → INCIP
  fig9b <- plot_simper(sim_den, sim_cover,
                       comparison  = "2024 Barrens_2025 Recovered Incipient Patches",
                       title_label = "BAR → INCIP")
  
  #INCIP → FOR
  fig9c <- plot_simper(sim_den, sim_cover,
                       comparison  = "2024 Incipients_2025 INCIP Recovered Forest Patches",
                       title_label = "INCIP → FOR")
  
  fig9a
  fig9b
  fig9c
  
################################################################################
#Looking at regional differences in patch transitions 
  
  #First, I need to create a new column that has Carmel or Monterey region designation 
  #adding region column based on site number
  regional_data <- mergedz %>%
    mutate(region = case_when(
      site %in% c("REC_01", "REC_02", "REC_03", "REC_04", "REC_05", "REC_10") ~ "Carmel",
      site %in% c("REC_06", "REC_07", "REC_11", "REC_12") ~ "Monterey",
      TRUE ~ NA_character_))

################################################################################
#HEATMAPS by Region
  
  #monterey
  monterey_transitions <- regional_data %>%
    filter(region == "Monterey")
  
  cross_tab_monterey <- table(
    `2024` = monterey_transitions$site_type_pred_2024,
    `2025` = monterey_transitions$site_type_pred_2025)
  
  #melt for ggplot
  cross_tab_melt_monterey <- reshape2::melt(cross_tab_monterey) %>%
    rename(Var1 = `2024`, Var2 = `2025`)
  
  #set factor order
  desired_order <- c("BAR", "INCIP", "FOR")
  
  cross_tab_melt_monterey <- cross_tab_melt_monterey %>%
    mutate(
      Var1 = factor(Var1, levels = desired_order),
      Var2 = factor(Var2, levels = desired_order))
  
  #heatmap
monterey_heatmap <- ggplot(cross_tab_melt_monterey, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile(color = "grey90", linewidth = 0.3) +
    geom_text(aes(label = value,
                  color = value > max(value) / 2),  # cleaner than ifelse in aes
              size = 5, fontface = "bold", show.legend = FALSE) +
    scale_color_manual(values = c("TRUE" = "white", "FALSE" = "black")) +
    scale_fill_gradientn(colors = brewer.pal(9, "Blues")) +  # 9 is max for Blues
    labs(
      x = "Predicted Patch Type — 2024",
      y = "Predicted Patch Type — 2025",
      title = "Monterey Patch Transitions",
      fill = "Count") +
    theme_minimal(base_size = 15) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 13),
      axis.text.y = element_text(face = "bold", size = 13),
      axis.title = element_text(face = "bold", size = 14),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18, 
                                margin = margin(b = 10)),
      legend.title = element_text(face = "bold", size = 13),
      legend.text = element_text(size = 12),
      panel.grid = element_blank())

  #carmel
carmel_transitions <- regional_data %>%
  filter(region == "Carmel")

cross_tab_carmel <- table(
  `2024` = carmel_transitions$site_type_pred_2024,
  `2025` = carmel_transitions$site_type_pred_2025)

  #melt for ggplot
cross_tab_melt_carmel <- reshape2::melt(cross_tab_carmel) %>%
  rename(Var1 = `2024`, Var2 = `2025`)

  #set factor order
desired_order <- c("BAR", "INCIP", "FOR")

cross_tab_melt_carmel <- cross_tab_melt_carmel %>%
  mutate(
    Var1 = factor(Var1, levels = desired_order),
    Var2 = factor(Var2, levels = desired_order))

  #heatmap
carmel_heatmap <- ggplot(cross_tab_melt_carmel, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "grey90", linewidth = 0.3) +
  geom_text(aes(label = value,
                color = value > max(value) / 2),  # cleaner than ifelse in aes
            size = 5, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = c("TRUE" = "white", "FALSE" = "black")) +
  scale_fill_gradientn(colors = brewer.pal(9, "Blues")) +  # 9 is max for Blues
  labs(
    x = "Predicted Patch Type — 2024",
    y = "Predicted Patch Type — 2025",
    title = "Carmel Patch Transitions",
    fill = "Count"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 13),
    axis.text.y = element_text(face = "bold", size = 13),
    axis.title = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18, 
                              margin = margin(b = 10)),
    legend.title = element_text(face = "bold", size = 13),
    legend.text = element_text(size = 12),
    panel.grid = element_blank())

################################################################################
#SIMPER by Region 

#adding region column to raw data based on site number
merged_data <- merged_data %>%
  mutate(region = case_when(
    site %in% c("REC_01", "REC_02", "REC_03", "REC_04", "REC_05", "REC_10") ~ "Carmel",
    site %in% c("REC_06", "REC_07", "REC_11", "REC_12") ~ "Monterey",
    TRUE ~ NA_character_))  

predictions_df <- mergedz %>%
  left_join(merged_data, by = c("site", "zone", "site_type")) %>%
  dplyr::select(c(site,zone,region, site_type_pred_2024,site_type_pred_2025, macro_stipe_density, density_ptecal, density_nerlue, density_lamset, density_eisarb, density_lamstump, density_macstump, lamr , macr , macj , nerj , ptej , lsetj , cov_articulated_coralline , cov_crustose_coralline , cov_encrusting_red , cov_fleshy_red , cov_stephanocystis , cov_dictyoneurum_spp  , cov_desmarestia_spp , cov_lam_holdfast_live ,cov_mac_holdfast_live,year))

patch_transitions <- predictions_df %>%
  mutate(
    #assigning the correct prediction for each year
    predicted_patch_type = case_when(
      year == 2024 ~ as.character(site_type_pred_2024),
      year == 2025 ~ as.character(site_type_pred_2025),
      TRUE ~ NA_character_),
    #creating the transition column
    patch_transition = paste(site_type_pred_2024, "→", site_type_pred_2025)
  ) %>%
  #dropping the old wide prediction columns, keep everything else (including species)
  dplyr::select(site, zone, region, year, predicted_patch_type, patch_transition, everything(), -site_type_pred_2024, -site_type_pred_2025)



#creating groupings that I want for future analyses 
patch_groupings <- patch_transitions %>%
  mutate(patches_groups = case_when(
    year == 2024 & predicted_patch_type == "BAR"   ~ "2024 Barrens",
    year == 2024 & predicted_patch_type == "INCIP" ~ "2024 Incipients",
    year == 2024 & predicted_patch_type == "FOR"   ~ "2024 Forests", 
    year == 2025 & patch_transition == "BAR → INCIP"  ~ "2025 Recovered Incipient Patches",
    year == 2025 & patch_transition == "BAR → FOR"    ~ "2025 BAR Recovered Forest Patches", #forests from barrens
    year == 2025 & patch_transition == "INCIP → FOR"  ~ "2025 INCIP Recovered Forest Patches", #forests from incipients
    year == 2025 & patch_transition == "BAR → BAR"    ~ "2025 Persistant Barren Patches",
    year == 2025 & patch_transition == "FOR → FOR"    ~ "2025 Persistant Forest Patches",
    year == 2025 & patch_transition == "INCIP → INCIP"~ "2025 Persistant Incipient Patches",
    year == 2025 & patch_transition == "INCIP → BAR"  ~ "2025 Regressed Barren Patches",
    year == 2025 & patch_transition == "FOR → BAR"    ~ "2025 Regressed Barren Patches",
    year == 2025 & patch_transition == "FOR → INCIP"  ~ "2025 Regressed Incipient Patches",
    TRUE ~ NA_character_))

################################################################################
#Regional SIMPER analysis - Carmel vs Monterey

#filtering patch_groupings by region then rebuilding comm matrices
build_regional_simper <- function(region_label) {
  
  groups_regional <- patch_groupings %>%
    filter(region == region_label) %>%
    dplyr::select(site, zone, year, patches_groups, 
                  macro_stipe_density:cov_mac_holdfast_live) %>%
    filter(!is.na(patches_groups)) %>%
    droplevels()
  
  #percent cover
  raw_cover_r <- groups_regional %>%
    dplyr::select(-c(macro_stipe_density:lsetj))
  cover_meta_r <- raw_cover_r %>%
    dplyr::select(site, zone, year, patches_groups)
  cover_comm_r <- raw_cover_r %>%
    dplyr::select(-c(site, zone, year, patches_groups))
  
  #density
  raw_density_r <- groups_regional %>%
    dplyr::select(site:lsetj)
  den_meta_r <- raw_density_r %>%
    dplyr::select(site, zone, year, patches_groups)
  den_comm_r <- raw_density_r %>%
    dplyr::select(-c(site, zone, year, patches_groups))
  den_comm_eps_r <- den_comm_r + 1e-6
  
  #check if there are enough groups to run simper (need at least 2)
  if (length(unique(cover_meta_r$patches_groups)) < 2) {
    warning(paste0("Not enough groups in ", region_label, " to run SIMPER"))
    return(NULL)
  }
  
  sim_cover_r <- simper(cover_comm_r, cover_meta_r$patches_groups, permutations = 999)
  sim_den_r   <- simper(den_comm_eps_r, den_meta_r$patches_groups, permutations = 999)
  
  return(list(sim_cover = sim_cover_r, sim_den = sim_den_r))
}

#running for each region
simper_carmel   <- build_regional_simper("Carmel")
simper_monterey <- build_regional_simper("Monterey")

#pulling out into global environment so summary() works on them
sim_den_carmel   <- simper_carmel$sim_den
sim_cover_carmel <- simper_carmel$sim_cover

sim_den_monterey   <- simper_monterey$sim_den
sim_cover_monterey <- simper_monterey$sim_cover

#SIMPER summary stats
summary(sim_den_carmel)
summary(sim_cover_carmel)
summary(sim_den_monterey)
summary(sim_cover_monterey)

################################################################################
#Fig. 11: Regional SIMPER - Carmel vs Monterey patch transitions (density + cover)

#Carmel figures

#BAR → FOR (Carmel)
fig11a_carmel <- plot_simper(simper_carmel$sim_den, simper_carmel$sim_cover,
                             comparison  = "2024 Barrens_2025 BAR Recovered Forest Patches",
                             title_label = "Carmel: BAR → FOR")
fig11a_carmel$plot

#BAR → INCIP (Carmel)
fig11b_carmel <- plot_simper(simper_carmel$sim_den, simper_carmel$sim_cover,
                             comparison  = "2024 Barrens_2025 Recovered Incipient Patches",
                             title_label = "Carmel: BAR → INCIP")
fig11b_carmel$plot

#INCIP → FOR (Carmel)
fig11c_carmel <- plot_simper(simper_carmel$sim_den, simper_carmel$sim_cover,
                             comparison  = "2024 Incipients_2025 INCIP Recovered Forest Patches",
                             title_label = "Carmel: INCIP → FOR")
fig11c_carmel$plot


#Monterey figures

#BAR → FOR (Monterey)
fig11a_monterey <- plot_simper(simper_monterey$sim_den, simper_monterey$sim_cover,
                               comparison  = "2024 Barrens_2025 BAR Recovered Forest Patches",
                               title_label = "Monterey: BAR → FOR")
fig11a_monterey$plot

#BAR → INCIP (Monterey)
fig11b_monterey <- plot_simper(simper_monterey$sim_den, simper_monterey$sim_cover,
                               comparison  = "2024 Barrens_2025 Recovered Incipient Patches",
                               title_label = "Monterey: BAR → INCIP")
fig11b_monterey$plot

#INCIP → FOR (Monterey)
fig11c_monterey <- plot_simper(simper_monterey$sim_den, simper_monterey$sim_cover,
                               comparison  = "2024 Incipients_2025 INCIP Recovered Forest Patches",
                               title_label = "Monterey: INCIP → FOR")
fig11c_monterey$plot


#combining Carmel and Monterey side by side per transition for easy comparison
fig11_bar_for <- fig11a_carmel$plot + fig11a_monterey$plot +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = "BAR → FOR: Carmel vs Monterey",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
      plot.tag   = element_text(face = "bold", size = 16)
    )
  )

fig11_bar_incip <- fig11b_carmel$plot + fig11b_monterey$plot +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = "BAR → INCIP: Carmel vs Monterey",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
      plot.tag   = element_text(face = "bold", size = 16)
    )
  )

fig11_incip_for <- fig11c_carmel$plot + fig11c_monterey$plot +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = "INCIP → FOR: Carmel vs Monterey",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
      plot.tag   = element_text(face = "bold", size = 16)
    )
  )

fig11_bar_for
fig11_bar_incip
fig11_incip_for



################################################################################
#SAVING ALL FIGURES THAT I WANT

ggsave(file.path(output2, "Fig_LDA_Model.png"),
       plot = a7, width = 10, height = 6, dpi = 300, bg = "white")

ggsave(file.path(output2, "Fig_LDA_Coefficients.png"),
       plot = b7, width = 10, height = 6, dpi = 300, bg = "white")

ggsave(file.path(output2, "Fig_Heatmap.png"),
       plot = fig8, width = 10, height = 6,  dpi = 300, bg = "white")

ggsave(file.path(output2, "Fig_SIMPER_BAR_INCIP.png"),
       plot = fig9b$plot, width = 9, height = 7,  dpi = 300, bg = "white")

ggsave(file.path(output2, "Fig_SIMPER_BAR_FOR.png"),
       plot = fig9a$plot, width = 9, height = 7,  dpi = 300, bg = "white")

ggsave(file.path(output2, "Fig_SIMPER_INCIP_FOR.png"),
       plot = fig9c$plot, width = 9, height = 7,  dpi = 300, bg = "white")

ggsave(file.path(output2, "Fig_Monterey_Heatmap.png"),
       plot = monterey_heatmap, width = 10, height = 6,  dpi = 300, bg = "white")

ggsave(file.path(output2, "Fig_Carmel_Heatmap.png"),
       plot = carmel_heatmap, width = 10, height = 6,  dpi = 300, bg = "white")

ggsave(file.path(output2, "Fig_Carmel_SIMPER_BAR_FOR.png"),
       plot = fig11a_carmel$plot, width = 9, height = 7, dpi = 300, bg = "white")

ggsave(file.path(output2, "Fig_Carmel_SIMPER_BAR_INCIP.png"),
       plot = fig11b_carmel$plot, width = 9, height = 7, dpi = 300, bg = "white")

ggsave(file.path(output2, "Fig_Carmel_SIMPER_INCIP_FOR.png"),
       plot = fig11c_carmel$plot, width = 9, height = 7, dpi = 300, bg = "white")

ggsave(file.path(output2, "Fig_Monterey_SIMPER_BAR_INCIP.png"),
       plot = fig11a_monterey$plot, width = 9, height = 7,  dpi = 300, bg = "white")

ggsave(file.path(output2, "Fig_Monterey_SIMPER_BAR_FOR.png"),
       plot = fig11b_monterey$plot, width = 9, height = 7,  dpi = 300, bg = "white")

ggsave(file.path(output2, "Fig_Monterey_SIMPER_INCIP_FOR.png"),
       plot = fig11c_monterey$plot, width = 9, height = 7,  dpi = 300, bg = "white")

################################################################################
#last write: April 25, 2026
