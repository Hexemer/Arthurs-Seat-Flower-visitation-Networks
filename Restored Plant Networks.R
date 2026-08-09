##################################
#Restored Plant Networks
##################################

# Aim 1:
# Determine which insects visit the restored Dianthus deltoides
# and Silene viscaria populations and the extent to which their
# visitors are shared with co-flowering plant species.
#
# Network structure is characterised using connectance,
# weighted NODF and modularity.


source("Overall Survey and Data.R")


#Prepare species-level interaction data

network_data <- pollinators %>%
  filter(
    Group != "None",
    Visits > 0,
    !is.na(Plant.ID),
    !is.na(Site),
    !is.na(Time.Period)
  ) %>%
  mutate(
    Species = trimws(
      as.character(Species)
    ),
    
    Genus = trimws(
      as.character(Genus)
    ),
    
    Group = trimws(
      as.character(Group)
    ),
    
    Plant.ID = factor(
      Plant.ID,
      levels = levels(
        plant_names$Plant.ID
      )
    )
  ) %>%
  
  left_join(
    plant_names,
    by = "Plant.ID"
  ) %>%
  
  mutate(
    Pollinator.Label = case_when(
      !is.na(Species) &
        Species != "" &
        Species != "NA" &
        Species != "Unknown" &
        Species != "None" ~
        
        sub(
          "^([A-Za-z])[A-Za-z]+ ",
          "\\1. ",
          Species
        ),
      
      !is.na(Genus) &
        Genus != "" ~
        
        paste0(
          Genus,
          " sp."
        ),
      
      TRUE ~
        Group
    ),
    
    Plant.Label = sub(
      "^([A-Za-z])[A-Za-z]+ ",
      "\\1. ",
      Plant.Species
    ),
    
    Site = factor(
      Site,
      levels = c(
        "1",
        "2"
      ),
      labels = c(
        "Site 1",
        "Site 2"
      )
    ),
    
    Time.Period = factor(
      Time.Period
    )
  ) %>%
  
  filter(
    !is.na(Plant.Species),
    !is.na(Pollinator.Label),
    Pollinator.Label != ""
  )


#Restored-plant visitor assemblages

focal_pollinators <- pollinators %>%
  filter(
    as.character(Plant.ID) %in% c(
      "2",
      "5",
      "10"
    ),
    Group != "None",
    Visits > 0
  ) %>%
  
  mutate(
    Focal_Plant = case_when(
      as.character(Plant.ID) == "2" ~
        "Dianthus deltoides, Site 1",
      
      as.character(Plant.ID) == "5" ~
        "Silene viscaria, Site 1",
      
      as.character(Plant.ID) == "10" ~
        "Silene viscaria, Site 2"
    ),
    
    Pollinator_Taxon = case_when(
      !is.na(Species) &
        Species != "" &
        Species != "None" &
        Species != "Unknown" ~
        
        Species,
      
      !is.na(Genus) &
        Genus != "" ~
        
        paste0(
          Genus,
          " sp."
        ),
      
      TRUE ~
        as.character(Group)
    )
  ) %>%
  
  group_by(
    Focal_Plant,
    Pollinator_Taxon,
    Group
  ) %>%
  
  summarise(
    Visits = sum(
      Visits,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  
  arrange(
    Focal_Plant,
    desc(Visits)
  )

focal_pollinators


#Restored-plant visitor richness and total visitation

focal_pollinator_summary <- focal_pollinators %>%
  group_by(
    Focal_Plant
  ) %>%
  
  summarise(
    Total.Visits = sum(
      Visits,
      na.rm = TRUE
    ),
    
    Pollinator.Taxa = n_distinct(
      Pollinator_Taxon
    ),
    
    .groups = "drop"
  )

focal_pollinator_summary


#Create plant x pollinator network matrices

make_species_network <- function(data) {
  
  if (
    nrow(data) == 0
  ) {
    
    return(
      matrix(
        numeric(0),
        nrow = 0,
        ncol = 0
      )
    )
  }
  
  web_table <- xtabs(
    Visits ~ Plant.Label + Pollinator.Label,
    data = data
  )
  
  web <- matrix(
    data = as.numeric(
      web_table
    ),
    nrow = nrow(
      web_table
    ),
    ncol = ncol(
      web_table
    ),
    dimnames = dimnames(
      web_table
    )
  )
  
  web <- web[
    rowSums(web) > 0,
    colSums(web) > 0,
    drop = FALSE
  ]
  
  web
}


overall_species_web <- network_data %>%
  make_species_network()


site1_species_web <- network_data %>%
  filter(
    Site == "Site 1"
  ) %>%
  make_species_network()


site2_species_web <- network_data %>%
  filter(
    Site == "Site 2"
  ) %>%
  make_species_network()


#Inspect site networks

site1_species_web
site2_species_web


#Network metric functions

clean_network_matrix <- function(web) {
  
  web <- as.matrix(
    web
  )
  
  storage.mode(
    web
  ) <- "numeric"
  
  web <- web[
    rowSums(web) > 0,
    colSums(web) > 0,
    drop = FALSE
  ]
  
  if (
    nrow(web) < 2 ||
    ncol(web) < 2
  ) {
    
    stop(
      "The network must contain at least two non-empty rows and columns."
    )
  }
  
  web
}


#Weighted NODF

calculate_weighted_nodf <- function(web) {
  
  web <- clean_network_matrix(
    web
  )
  
  as.numeric(
    bipartite::networklevel(
      web,
      index = "weighted NODF"
    )
  )
}


#Connectance

calculate_connectance <- function(web) {
  
  web <- clean_network_matrix(
    web
  )
  
  as.numeric(
    bipartite::networklevel(
      web,
      index = "connectance"
    )
  )
}


#Modularity

calculate_modularity <- function(
    web,
    module_reps = 100
) {
  
  web <- clean_network_matrix(
    web
  )
  
  result <- bipartite::metaComputeModules(
    web,
    N = module_reps,
    method = "Beckett"
  )
  
  as.numeric(
    result@likelihood
  )
}


#Site-level network metrics

set.seed(
  123
)

site1_modularity <- calculate_modularity(
  site1_species_web,
  module_reps = 100
)

set.seed(
  456
)

site2_modularity <- calculate_modularity(
  site2_species_web,
  module_reps = 100
)


site_metrics <- tibble(
  Site = c(
    "Site 1",
    "Site 2"
  ),
  
  Plant_Nodes = c(
    nrow(
      site1_species_web
    ),
    nrow(
      site2_species_web
    )
  ),
  
  Pollinator_Nodes = c(
    ncol(
      site1_species_web
    ),
    ncol(
      site2_species_web
    )
  ),
  
  Total_Visits = c(
    sum(
      site1_species_web
    ),
    sum(
      site2_species_web
    )
  ),
  
  Number_of_Links = c(
    sum(
      site1_species_web > 0
    ),
    sum(
      site2_species_web > 0
    )
  ),
  
  Connectance = c(
    calculate_connectance(
      site1_species_web
    ),
    calculate_connectance(
      site2_species_web
    )
  ),
  
  Weighted_NODF = c(
    calculate_weighted_nodf(
      site1_species_web
    ),
    calculate_weighted_nodf(
      site2_species_web
    )
  ),
  
  Modularity_Q = c(
    site1_modularity,
    site2_modularity
  )
) %>%
  
  mutate(
    Connectance = round(
      Connectance,
      3
    ),
    
    Weighted_NODF = round(
      Weighted_NODF,
      2
    ),
    
    Modularity_Q = round(
      Modularity_Q,
      3
    )
  )

site_metrics


#Null-model helper functions

calculate_ses <- function(
    observed,
    null_values
) {
  
  null_values <- null_values[
    is.finite(
      null_values
    )
  ]
  
  null_sd <- sd(
    null_values
  )
  
  if (
    length(null_values) < 2 ||
    !is.finite(null_sd) ||
    null_sd == 0
  ) {
    
    return(
      NA_real_
    )
  }
  
  (
    observed -
      mean(
        null_values
      )
  ) /
    null_sd
}


calculate_two_tailed_p <- function(
    observed,
    null_values
) {
  
  null_values <- null_values[
    is.finite(
      null_values
    )
  ]
  
  if (
    length(null_values) == 0
  ) {
    
    return(
      NA_real_
    )
  }
  
  lower_p <- (
    sum(
      null_values <= observed
    ) + 1
  ) /
    (
      length(null_values) + 1
    )
  
  upper_p <- (
    sum(
      null_values >= observed
    ) + 1
  ) /
    (
      length(null_values) + 1
    )
  
  min(
    1,
    2 * min(
      lower_p,
      upper_p
    )
  )
}


#Weighted NODF null-model standardisation

standardise_weighted_nodf <- function(
    web,
    null_reps = 500,
    seed = 123
) {
  
  web <- clean_network_matrix(
    web
  )
  
  observed <- calculate_weighted_nodf(
    web
  )
  
  set.seed(
    seed
  )
  
  r2d_nulls <- bipartite::nullmodel(
    web,
    N = null_reps,
    method = "r2dtable"
  )
  
  r2d_values <- vapply(
    r2d_nulls,
    calculate_weighted_nodf,
    numeric(1)
  )
  
  set.seed(
    seed
  )
  
  vaz_nulls <- bipartite::nullmodel(
    web,
    N = null_reps,
    method = "vaznull"
  )
  
  vaz_values <- vapply(
    vaz_nulls,
    calculate_weighted_nodf,
    numeric(1)
  )
  
  summary_table <- tibble(
    Metric = "Weighted NODF",
    Observed = observed,
    
    Null_Mean_R2d = mean(
      r2d_values,
      na.rm = TRUE
    ),
    
    Null_SD_R2d = sd(
      r2d_values,
      na.rm = TRUE
    ),
    
    z_R2d = calculate_ses(
      observed,
      r2d_values
    ),
    
    p_R2d = calculate_two_tailed_p(
      observed,
      r2d_values
    ),
    
    Null_Mean_VazNull = mean(
      vaz_values,
      na.rm = TRUE
    ),
    
    Null_SD_VazNull = sd(
      vaz_values,
      na.rm = TRUE
    ),
    
    z_VazNull = calculate_ses(
      observed,
      vaz_values
    ),
    
    p_VazNull = calculate_two_tailed_p(
      observed,
      vaz_values
    )
  )
  
  list(
    summary = summary_table,
    r2d_values = r2d_values,
    vaz_values = vaz_values
  )
}


#Modularity null-model standardisation

standardise_modularity <- function(
    web,
    null_reps = 500,
    module_reps = 500,
    seed = 123
) {
  
  web <- clean_network_matrix(
    web
  )
  
  observed <- calculate_modularity(
    web,
    module_reps = module_reps
  )
  
  set.seed(
    seed
  )
  
  r2d_nulls <- bipartite::nullmodel(
    web,
    N = null_reps,
    method = "r2dtable"
  )
  
  r2d_values <- vapply(
    r2d_nulls,
    
    function(null_web) {
      
      calculate_modularity(
        null_web,
        module_reps = module_reps
      )
    },
    
    numeric(1)
  )
  
  set.seed(
    seed
  )
  
  vaz_nulls <- bipartite::nullmodel(
    web,
    N = null_reps,
    method = "vaznull"
  )
  
  vaz_values <- vapply(
    vaz_nulls,
    
    function(null_web) {
      
      calculate_modularity(
        null_web,
        module_reps = module_reps
      )
    },
    
    numeric(1)
  )
  
  summary_table <- tibble(
    Metric = "Modularity Q",
    Observed = observed,
    
    Null_Mean_R2d = mean(
      r2d_values,
      na.rm = TRUE
    ),
    
    Null_SD_R2d = sd(
      r2d_values,
      na.rm = TRUE
    ),
    
    z_R2d = calculate_ses(
      observed,
      r2d_values
    ),
    
    p_R2d = calculate_two_tailed_p(
      observed,
      r2d_values
    ),
    
    Null_Mean_VazNull = mean(
      vaz_values,
      na.rm = TRUE
    ),
    
    Null_SD_VazNull = sd(
      vaz_values,
      na.rm = TRUE
    ),
    
    z_VazNull = calculate_ses(
      observed,
      vaz_values
    ),
    
    p_VazNull = calculate_two_tailed_p(
      observed,
      vaz_values
    )
  )
  
  list(
    summary = summary_table,
    r2d_values = r2d_values,
    vaz_values = vaz_values
  )
}


#Connectance null-model standardisation

# VazNull preserves connectance, so only R2d is used for SES.

standardise_connectance <- function(
    web,
    null_reps = 500,
    seed = 123
) {
  
  web <- clean_network_matrix(
    web
  )
  
  observed <- calculate_connectance(
    web
  )
  
  set.seed(
    seed
  )
  
  r2d_nulls <- bipartite::nullmodel(
    web,
    N = null_reps,
    method = "r2dtable"
  )
  
  r2d_values <- vapply(
    r2d_nulls,
    calculate_connectance,
    numeric(1)
  )
  
  summary_table <- tibble(
    Metric = "Connectance",
    Observed = observed,
    
    Null_Mean_R2d = mean(
      r2d_values,
      na.rm = TRUE
    ),
    
    Null_SD_R2d = sd(
      r2d_values,
      na.rm = TRUE
    ),
    
    z_R2d = calculate_ses(
      observed,
      r2d_values
    ),
    
    p_R2d = calculate_two_tailed_p(
      observed,
      r2d_values
    ),
    
    Null_Mean_VazNull = NA_real_,
    Null_SD_VazNull = NA_real_,
    z_VazNull = NA_real_,
    p_VazNull = NA_real_
  )
  
  list(
    summary = summary_table,
    r2d_values = r2d_values
  )
}


#Run null models

# 500 null networks are generated for each null model.
# Modularity uses repeated module detection.


#Site 1

site1_nodf_results <- standardise_weighted_nodf(
  site1_species_web,
  null_reps = 500,
  seed = 123
)

site1_modularity_results <- standardise_modularity(
  site1_species_web,
  null_reps = 500,
  module_reps = 500,
  seed = 123
)

site1_connectance_results <- standardise_connectance(
  site1_species_web,
  null_reps = 500,
  seed = 123
)


#Site 2

site2_nodf_results <- standardise_weighted_nodf(
  site2_species_web,
  null_reps = 500,
  seed = 456
)

site2_modularity_results <- standardise_modularity(
  site2_species_web,
  null_reps = 500,
  module_reps = 500,
  seed = 456
)

site2_connectance_results <- standardise_connectance(
  site2_species_web,
  null_reps = 500,
  seed = 456
)


#Combine null-model results

site1_results <- bind_rows(
  site1_modularity_results$summary,
  site1_connectance_results$summary,
  site1_nodf_results$summary
) %>%
  
  mutate(
    Site = "Site 1",
    .before = 1
  )


site2_results <- bind_rows(
  site2_modularity_results$summary,
  site2_connectance_results$summary,
  site2_nodf_results$summary
) %>%
  
  mutate(
    Site = "Site 2",
    .before = 1
  )


network_null_results <- bind_rows(
  site1_results,
  site2_results
)

network_null_results


#Rounded null-model reporting table

network_null_results_report <- network_null_results %>%
  
  mutate(
    Observed = round(
      Observed,
      3
    ),
    
    Null_Mean_R2d = round(
      Null_Mean_R2d,
      3
    ),
    
    Null_SD_R2d = round(
      Null_SD_R2d,
      3
    ),
    
    z_R2d = round(
      z_R2d,
      2
    ),
    
    p_R2d = round(
      p_R2d,
      3
    ),
    
    Null_Mean_VazNull = round(
      Null_Mean_VazNull,
      3
    ),
    
    Null_SD_VazNull = round(
      Null_SD_VazNull,
      3
    ),
    
    z_VazNull = round(
      z_VazNull,
      2
    ),
    
    p_VazNull = round(
      p_VazNull,
      3
    )
  )

network_null_results_report

#Plant identity and attractiveness

# Test whether differences in visitation among plant species
# remain after accounting for the number of flowers available.

plant_analysis_data <- survey_data %>%
  filter(
    !(Plant.ID %in% c(
      "6",
      "13"
    )),
    !is.na(Number.of.Flowers),
    Number.of.Flowers > 0,
    !is.na(Plant.ID),
    !is.na(Site),
    !is.na(Time.Period)
  ) %>%
  
  mutate(
    Plant.ID = factor(
      Plant.ID,
      levels = levels(
        plant_names$Plant.ID
      )
    )
  ) %>%
  
  left_join(
    plant_names,
    by = "Plant.ID"
  ) %>%
  
  filter(
    !is.na(Plant.Species)
  ) %>%
  
  mutate(
    Plant.Species = droplevels(
      factor(
        Plant.Species
      )
    ),
    
    Site = droplevels(
      factor(
        Site
      )
    ),
    
    Time.Period = droplevels(
      factor(
        Time.Period
      )
    ),
    
    Plant.Common.Name = case_when(
      Plant.Species ==
        "Hypochaeris radicata" ~
        "Yellow flatweed",
      
      Plant.Species ==
        "Dianthus deltoides" ~
        "Maiden pink",
      
      Plant.Species ==
        "Thymus drucei" ~
        "Thyme",
      
      Plant.Species ==
        "Lotus corniculatus" ~
        "Bird's-foot trefoil",
      
      Plant.Species ==
        "Silene viscaria" ~
        "Sticky catchfly",
      
      Plant.Species ==
        "Teucrium scorodonia" ~
        "Woodland germander",
      
      Plant.Species ==
        "Galium verum" ~
        "Lady's bedstraw",
      
      Plant.Species ==
        "Cirsium arvense" ~
        "Creeping thistle",
      
      TRUE ~
        as.character(
          Plant.Species
        )
    )
  )


#Check plant representation

table(
  plant_analysis_data$Plant.Species
)


#Descriptive visitation by plant species

plant_descriptive_summary <- plant_analysis_data %>%
  group_by(
    Plant.Species
  ) %>%
  
  summarise(
    Surveys = n(),
    
    Total.Visits = sum(
      Visits,
      na.rm = TRUE
    ),
    
    Mean.Flowers = mean(
      Number.of.Flowers,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

plant_descriptive_summary


#Fit abundance-only model

# This model assumes visitation is explained by floral availability,
# site and time period, without an additional plant-species effect.

ModelAbundance_nb <- MASS::glm.nb(
  Visits ~
    Site +
    Time.Period +
    offset(
      log(
        Number.of.Flowers
      )
    ),
  
  data =
    plant_analysis_data
)


#Fit plant-species model

# The offset accounts for the number of flowers available.
# Plant.Species therefore tests whether plant identity explains
# visitation beyond differences in floral abundance.

ModelSpecies_nb <- MASS::glm.nb(
  Visits ~
    Plant.Species +
    Site +
    Time.Period +
    offset(
      log(
        Number.of.Flowers
      )
    ),
  
  data =
    plant_analysis_data
)


summary(
  ModelSpecies_nb
)


#Test overall plant-species effect

Plant_Anova <- car::Anova(
  ModelSpecies_nb,
  type = 2
)

Plant_Anova


#Does adding plant species improve model fit?

Plant_model_comparison <- anova(
  ModelAbundance_nb,
  ModelSpecies_nb,
  test = "Chisq"
)

Plant_model_comparison


#Compare AIC

AIC(
  ModelAbundance_nb,
  ModelSpecies_nb
)


#Estimated attractiveness of each plant species

# offset = 0 corresponds to log(1 flower).
# Estimates are therefore predicted visits per flower
# per 20-minute survey.

Species_emmeans <- emmeans(
  ModelSpecies_nb,
  ~ Plant.Species,
  type = "response",
  offset = 0
)

Species_emmeans


#Pairwise comparisons among plant species

Species_pairs <- pairs(
  Species_emmeans,
  adjust = "tukey"
)

Species_pairs


#Create plant-attractiveness results table

Species_table <- as.data.frame(
  Species_emmeans
) %>%
  
  left_join(
    plant_analysis_data %>%
      dplyr::distinct(
        Plant.Species,
        Plant.Common.Name
      ),
    
    by = "Plant.Species"
  ) %>%
  
  mutate(
    Type = case_when(
      Plant.Species %in% c(
        "Dianthus deltoides",
        "Silene viscaria"
      ) ~
        "Restored",
      
      TRUE ~
        "Co-flowering"
    )
  ) %>%
  
  dplyr::select(
    Plant.Common.Name,
    Plant.Species,
    Type,
    response,
    SE,
    asymp.LCL,
    asymp.UCL
  ) %>%
  
  rename(
    `Common name` =
      Plant.Common.Name,
    
    `Scientific name` =
      Plant.Species,
    
    `Plant type` =
      Type,
    
    `Predicted visits per flower per 20-minute survey` =
      response,
    
    `Standard error` =
      SE,
    
    `Lower 95% CI` =
      asymp.LCL,
    
    `Upper 95% CI` =
      asymp.UCL
  ) %>%
  
  arrange(
    desc(
      `Predicted visits per flower per 20-minute survey`
    )
  ) %>%
  
  mutate(
    across(
      where(
        is.numeric
      ),
      ~ round(
        .x,
        3
      )
    )
  )

Species_table


#Create plant-attractiveness figure

Species_plot_data <- as.data.frame(
  Species_emmeans
) %>%
  
  left_join(
    plant_analysis_data %>%
      dplyr::distinct(
        Plant.Species,
        Plant.Common.Name
      ),
    
    by = "Plant.Species"
  ) %>%
  
  mutate(
    Type = case_when(
      Plant.Species %in% c(
        "Dianthus deltoides",
        "Silene viscaria"
      ) ~
        "Restored",
      
      TRUE ~
        "Co-flowering"
    )
  ) %>%
  
  arrange(
    response
  ) %>%
  
  mutate(
    Plant.Label = paste0(
      Plant.Common.Name,
      "\n",
      Plant.Species
    ),
    
    Plant.Label = factor(
      Plant.Label,
      levels = Plant.Label
    )
  )


Plant_attractiveness_plot <- ggplot(
  Species_plot_data,
  
  aes(
    x = Plant.Label,
    y = response,
    fill = Type
  )
) +
  
  geom_col(
    width = 0.7
  ) +
  
  geom_errorbar(
    aes(
      ymin = asymp.LCL,
      ymax = asymp.UCL
    ),
    
    width = 0.2,
    linewidth = 0.7
  ) +
  
  scale_fill_manual(
    values = c(
      "Restored" = "#E78AC3",
      "Co-flowering" = "#66C2A5"
    )
  ) +
  
  scale_y_log10() +
  
  labs(
    x = "Plant species",
    y = "Predicted visits per flower per 20-minute survey",
    fill = NULL
  ) +
  
  theme_classic(
    base_size = 14
  ) +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    
    legend.position = "top"
  )


Plant_attractiveness_plot


#Pollinator sharing with co-flowering species

# Pollinator sharing is defined as the same pollinator taxon
# visiting a restored and co-flowering species at the same site
# during the same time period.

restored_species <- c(
  "Dianthus deltoides",
  "Silene viscaria"
)


#Create plant-occurrence dataset

plant_occasion <- pollinators %>%
  
  mutate(
    Plant.ID = factor(
      Plant.ID,
      levels = levels(
        plant_names$Plant.ID
      )
    )
  ) %>%
  
  left_join(
    plant_names,
    by = "Plant.ID"
  ) %>%
  
  mutate(
    Site = case_when(
      as.character(Site) %in% c(
        "1",
        "Site 1"
      ) ~ "Site 1",
      
      as.character(Site) %in% c(
        "2",
        "Site 2"
      ) ~ "Site 2",
      
      TRUE ~
        NA_character_
    ),
    
    Time.Period = trimws(
      as.character(
        Time.Period
      )
    )
  ) %>%
  
  filter(
    !is.na(Site),
    !is.na(Time.Period),
    Time.Period != "",
    !is.na(Plant.Species)
  ) %>%
  
  dplyr::select(
    Site,
    Time.Period,
    Plant.Species
  ) %>%
  
  distinct()


#Create positive-interaction dataset

interaction_data <- pollinators %>%
  
  mutate(
    Plant.ID = factor(
      Plant.ID,
      levels = levels(
        plant_names$Plant.ID
      )
    )
  ) %>%
  
  left_join(
    plant_names,
    by = "Plant.ID"
  ) %>%
  
  mutate(
    Site = case_when(
      as.character(Site) %in% c(
        "1",
        "Site 1"
      ) ~ "Site 1",
      
      as.character(Site) %in% c(
        "2",
        "Site 2"
      ) ~ "Site 2",
      
      TRUE ~
        NA_character_
    ),
    
    Time.Period = trimws(
      as.character(
        Time.Period
      )
    ),
    
    Pollinator.Taxon = case_when(
      !is.na(Species) &
        Species != "" &
        Species != "NA" &
        Species != "Unknown" &
        Species != "None" ~
        
        trimws(
          as.character(
            Species
          )
        ),
      
      !is.na(Genus) &
        Genus != "" ~
        
        paste0(
          trimws(
            as.character(
              Genus
            )
          ),
          " sp."
        ),
      
      TRUE ~
        trimws(
          as.character(
            Group
          )
        )
    )
  ) %>%
  
  filter(
    !is.na(Site),
    !is.na(Time.Period),
    Time.Period != "",
    !is.na(Plant.Species),
    !is.na(Pollinator.Taxon),
    Pollinator.Taxon != "",
    Group != "None",
    !is.na(Visits),
    Visits > 0
  ) %>%
  
  dplyr::select(
    Site,
    Time.Period,
    Plant.Species,
    Pollinator.Taxon,
    Visits
  )


#Identify periods when restored species were present

restored_periods <- plant_occasion %>%
  filter(
    Plant.Species %in%
      restored_species
  ) %>%
  
  transmute(
    Site,
    Time.Period,
    Restored.Species =
      Plant.Species
  ) %>%
  
  distinct()


#Identify periods when co-flowering species were present

native_periods <- plant_occasion %>%
  filter(
    !Plant.Species %in%
      restored_species
  ) %>%
  
  transmute(
    Site,
    Time.Period,
    Native.Species =
      Plant.Species
  ) %>%
  
  distinct()


#Calculate temporal overlap

temporal_overlap <- restored_periods %>%
  
  inner_join(
    native_periods,
    by = c(
      "Site",
      "Time.Period"
    )
  )


overlap_summary <- temporal_overlap %>%
  
  group_by(
    Site,
    Restored.Species,
    Native.Species
  ) %>%
  
  summarise(
    Overlapping.Time.Periods = n_distinct(
      Time.Period
    ),
    
    .groups = "drop"
  )


#Summarise restored-plant interactions

restored_period_interactions <- interaction_data %>%
  
  filter(
    Plant.Species %in%
      restored_species
  ) %>%
  
  group_by(
    Site,
    Time.Period,
    Plant.Species,
    Pollinator.Taxon
  ) %>%
  
  summarise(
    Restored.Visits = sum(
      Visits,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
  
  rename(
    Restored.Species =
      Plant.Species
  )


#Summarise co-flowering plant interactions

native_period_interactions <- interaction_data %>%
  
  filter(
    !Plant.Species %in%
      restored_species
  ) %>%
  
  group_by(
    Site,
    Time.Period,
    Plant.Species,
    Pollinator.Taxon
  ) %>%
  
  summarise(
    Native.Visits = sum(
      Visits,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
  
  rename(
    Native.Species =
      Plant.Species
  )


#Identify shared pollinator use

period_sharing <- restored_period_interactions %>%
  
  inner_join(
    native_period_interactions,
    by = c(
      "Site",
      "Time.Period",
      "Pollinator.Taxon"
    )
  )


#Summarise sharing by restored-co-flowering pair

shared_pollinator_summary <- period_sharing %>%
  
  group_by(
    Site,
    Restored.Species,
    Native.Species
  ) %>%
  
  summarise(
    Shared.Pollinator.Taxa = n_distinct(
      Pollinator.Taxon
    ),
    
    Time.Period.Sharing.Events = n_distinct(
      interaction(
        Time.Period,
        Pollinator.Taxon,
        drop = TRUE
      )
    ),
    
    Time.Periods.With.Shared.Pollinators = n_distinct(
      Time.Period
    ),
    
    Total.Restored.Visits = sum(
      Restored.Visits,
      na.rm = TRUE
    ),
    
    Total.Native.Visits = sum(
      Native.Visits,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


#Complete restored-co-flowering interaction table

potential_interactions <- overlap_summary %>%
  
  left_join(
    shared_pollinator_summary,
    by = c(
      "Site",
      "Restored.Species",
      "Native.Species"
    )
  ) %>%
  
  mutate(
    across(
      c(
        Shared.Pollinator.Taxa,
        Time.Period.Sharing.Events,
        Time.Periods.With.Shared.Pollinators,
        Total.Restored.Visits,
        Total.Native.Visits
      ),
      ~ replace_na(
        .x,
        0
      )
    ),
    
    Proportion.Time.Periods.With.Sharing = if_else(
      Overlapping.Time.Periods > 0,
      
      Time.Periods.With.Shared.Pollinators /
        Overlapping.Time.Periods,
      
      NA_real_
    ),
    
    Shared.Any.Pollinators =
      Shared.Pollinator.Taxa > 0
  ) %>%
  
  arrange(
    Site,
    Restored.Species,
    desc(
      Shared.Pollinator.Taxa
    ),
    desc(
      Time.Period.Sharing.Events
    ),
    Native.Species
  )


#View complete pollinator-sharing results

print(
  potential_interactions,
  n = Inf
)


#Species sharing at least one pollinator

native_species_with_sharing <- potential_interactions %>%
  filter(
    Shared.Any.Pollinators
  ) %>%
  
  distinct(
    Native.Species
  ) %>%
  
  arrange(
    Native.Species
  )

native_species_with_sharing


#Reporting table

pollinator_sharing_table <- potential_interactions %>%
  
  dplyr::select(
    Site,
    Restored.Species,
    Native.Species,
    Overlapping.Time.Periods,
    Shared.Pollinator.Taxa,
    Time.Period.Sharing.Events,
    Time.Periods.With.Shared.Pollinators,
    Proportion.Time.Periods.With.Sharing,
    Shared.Any.Pollinators
  )

print(
  pollinator_sharing_table,
  n = Inf
)


#Pollinator-sharing figure

sharing_plot_data <- potential_interactions %>%
  filter(
    Shared.Any.Pollinators
  ) %>%
  
  mutate(
    Site = factor(
      Site,
      levels = c(
        "Site 1",
        "Site 2"
      )
    ),
    
    Restored.Species = factor(
      Restored.Species,
      levels = c(
        "Dianthus deltoides",
        "Silene viscaria"
      )
    ),
    
    Native.Species = reorder(
      Native.Species,
      Shared.Pollinator.Taxa
    )
  )


Pollinator_sharing_plot <- ggplot(
  sharing_plot_data,
  aes(
    x = Shared.Pollinator.Taxa,
    y = Native.Species,
    size = Time.Period.Sharing.Events
  )
) +
  
  geom_point(
    shape = 21,
    fill = "plum1",
    colour = "grey20",
    stroke = 0.7,
    alpha = 0.9
  ) +
  
  facet_grid(
    Restored.Species ~ Site,
    scales = "free_y",
    space = "free_y"
  ) +
  
  scale_x_continuous(
    breaks = scales::pretty_breaks()
  ) +
  
  scale_size_continuous(
    name = "Sharing events",
    range = c(
      3,
      10
    )
  ) +
  
  labs(
    x = "Number of shared pollinator taxa",
    y = "Co-flowering plant species"
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    strip.text = element_text(
      face = "bold"
    ),
    
    axis.text.y = element_text(
      face = "italic"
    ),
    
    legend.title = element_text(
      face = "bold"
    )
  )

Pollinator_sharing_plot


#Final Aim 1 outputs

focal_pollinators

focal_pollinator_summary

plant_descriptive_summary

Plant_Anova

Plant_model_comparison

Species_emmeans

Species_pairs

Species_table

Plant_attractiveness_plot

site_metrics

network_null_results_report

pollinator_sharing_table

Pollinator_sharing_plot