##################################
#Aim 3 Temporal Variation and Networks
##################################

# Aim 3:
# Examine how pollinator visitation varies throughout the day,
# how plant-pollinator network structure changes among time periods,
# and what the timing of shared pollinator use reveals about
# potential interactions between restored and co-flowering plants.


source("Overall Survey and Data.R")

source("Visitation Time and Site.R")


#Time-period visitation

# Time.Period is tested in the shared Time.Period + Site GLMM
# fitted in Visitation Time and Site.R.

TimeSite_Anova


#Estimated visitation by time period

Time_emmeans$emmeans


#Pairwise comparisons among time periods

Time_emmeans$contrasts


#Descriptive visitation by time period

time_descriptive_summary <- survey_data %>%
  group_by(
    Time.Period
  ) %>%
  summarise(
    Surveys = n(),
    
    Total.Visits = sum(
      Visits,
      na.rm = TRUE
    ),
    
    Mean.Visits = mean(
      Visits,
      na.rm = TRUE
    ),
    
    Median.Visits = median(
      Visits,
      na.rm = TRUE
    ),
    
    Zero.Visit.Surveys = sum(
      Visits == 0
    ),
    
    .groups = "drop"
  )

time_descriptive_summary


#Descriptive visitation by site and time period

site_time_visitation <- survey_data %>%
  group_by(
    Site,
    Time.Period
  ) %>%
  summarise(
    Surveys = n(),
    
    Total.Visits = sum(
      Visits,
      na.rm = TRUE
    ),
    
    Mean.Visits = mean(
      Visits,
      na.rm = TRUE
    ),
    
    Median.Visits = median(
      Visits,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

site_time_visitation


#Prepare species-level network data

species_network_data <- pollinators %>%
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
      Time.Period,
      levels = c(
        "1",
        "2",
        "3"
      )
    )
  ) %>%
  
  filter(
    !is.na(Plant.Species),
    !is.na(Pollinator.Label),
    Pollinator.Label != ""
  )


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


#Site 1 temporal networks

site1_tp1_web <- species_network_data %>%
  filter(
    Site == "Site 1",
    Time.Period == "1"
  ) %>%
  make_species_network()


site1_tp2_web <- species_network_data %>%
  filter(
    Site == "Site 1",
    Time.Period == "2"
  ) %>%
  make_species_network()


site1_tp3_web <- species_network_data %>%
  filter(
    Site == "Site 1",
    Time.Period == "3"
  ) %>%
  make_species_network()


#Site 2 temporal networks

site2_tp1_web <- species_network_data %>%
  filter(
    Site == "Site 2",
    Time.Period == "1"
  ) %>%
  make_species_network()


site2_tp2_web <- species_network_data %>%
  filter(
    Site == "Site 2",
    Time.Period == "2"
  ) %>%
  make_species_network()


site2_tp3_web <- species_network_data %>%
  filter(
    Site == "Site 2",
    Time.Period == "3"
  ) %>%
  make_species_network()


#Inspect temporal matrices

site1_tp1_web
site1_tp2_web
site1_tp3_web

site2_tp1_web
site2_tp2_web
site2_tp3_web


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

calculate_modularity_meta <- function(
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


#Site x time-period network metrics

set.seed(
  123
)

site_time_metrics <- tibble(
  Site = c(
    "Site 1",
    "Site 1",
    "Site 1",
    "Site 2",
    "Site 2",
    "Site 2"
  ),
  
  Time.Period = factor(
    c(
      "1",
      "2",
      "3",
      "1",
      "2",
      "3"
    ),
    levels = c(
      "1",
      "2",
      "3"
    )
  ),
  
  Plant_Nodes = c(
    nrow(site1_tp1_web),
    nrow(site1_tp2_web),
    nrow(site1_tp3_web),
    nrow(site2_tp1_web),
    nrow(site2_tp2_web),
    nrow(site2_tp3_web)
  ),
  
  Pollinator_Nodes = c(
    ncol(site1_tp1_web),
    ncol(site1_tp2_web),
    ncol(site1_tp3_web),
    ncol(site2_tp1_web),
    ncol(site2_tp2_web),
    ncol(site2_tp3_web)
  ),
  
  Total_Visits = c(
    sum(site1_tp1_web),
    sum(site1_tp2_web),
    sum(site1_tp3_web),
    sum(site2_tp1_web),
    sum(site2_tp2_web),
    sum(site2_tp3_web)
  ),
  
  Number_of_Links = c(
    sum(site1_tp1_web > 0),
    sum(site1_tp2_web > 0),
    sum(site1_tp3_web > 0),
    sum(site2_tp1_web > 0),
    sum(site2_tp2_web > 0),
    sum(site2_tp3_web > 0)
  ),
  
  Connectance = c(
    calculate_connectance(
      site1_tp1_web
    ),
    
    calculate_connectance(
      site1_tp2_web
    ),
    
    calculate_connectance(
      site1_tp3_web
    ),
    
    calculate_connectance(
      site2_tp1_web
    ),
    
    calculate_connectance(
      site2_tp2_web
    ),
    
    calculate_connectance(
      site2_tp3_web
    )
  ),
  
  Weighted_NODF = c(
    calculate_weighted_nodf(
      site1_tp1_web
    ),
    
    calculate_weighted_nodf(
      site1_tp2_web
    ),
    
    calculate_weighted_nodf(
      site1_tp3_web
    ),
    
    calculate_weighted_nodf(
      site2_tp1_web
    ),
    
    calculate_weighted_nodf(
      site2_tp2_web
    ),
    
    calculate_weighted_nodf(
      site2_tp3_web
    )
  ),
  
  Modularity_Q = c(
    calculate_modularity_meta(
      site1_tp1_web
    ),
    
    calculate_modularity_meta(
      site1_tp2_web
    ),
    
    calculate_modularity_meta(
      site1_tp3_web
    ),
    
    calculate_modularity_meta(
      site2_tp1_web
    ),
    
    calculate_modularity_meta(
      site2_tp2_web
    ),
    
    calculate_modularity_meta(
      site2_tp3_web
    )
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


site_time_metrics


#Temporal pollinator sharing

# Sharing is only counted when the same pollinator taxon
# visited a restored and co-flowering species at the same
# site during the same time period.

restored_species <- c(
  "Dianthus deltoides",
  "Silene viscaria"
)


#Plant occurrence data

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


#Positive interactions

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


#Restored species present in each period

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


#Co-flowering species present in each period

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


#Potential temporal overlap

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


#Restored-plant interactions

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


#Co-flowering plant interactions

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


#Shared pollinator use within the same time period

period_sharing <- restored_period_interactions %>%
  
  inner_join(
    native_period_interactions,
    by = c(
      "Site",
      "Time.Period",
      "Pollinator.Taxon"
    )
  )


#Summarise temporal sharing by plant pair

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
    
    .groups = "drop"
  )


#Complete temporal-sharing table

temporal_sharing_table <- overlap_summary %>%
  
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
        Time.Periods.With.Shared.Pollinators
      ),
      
      ~ replace_na(
        .x,
        0
      )
    ),
    
    Proportion.Time.Periods.With.Sharing =
      Time.Periods.With.Shared.Pollinators /
      Overlapping.Time.Periods,
    
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
    )
  )


print(
  temporal_sharing_table,
  n = Inf
)


#Summarise sharing separately for each time period

sharing_by_time_period <- period_sharing %>%
  
  group_by(
    Site,
    Restored.Species,
    Time.Period
  ) %>%
  
  summarise(
    Shared.Pollinator.Taxa = n_distinct(
      Pollinator.Taxon
    ),
    
    Native.Plant.Species = n_distinct(
      Native.Species
    ),
    
    Sharing.Events = n_distinct(
      interaction(
        Native.Species,
        Pollinator.Taxon,
        drop = TRUE
      )
    ),
    
    .groups = "drop"
  ) %>%
  
  arrange(
    Site,
    Restored.Species,
    Time.Period
  )


sharing_by_time_period


#Final Aim 3 outputs

TimeSite_Anova

Time_emmeans$emmeans

Time_emmeans$contrasts

time_descriptive_summary

site_time_visitation

site_time_metrics

temporal_sharing_table

sharing_by_time_period
