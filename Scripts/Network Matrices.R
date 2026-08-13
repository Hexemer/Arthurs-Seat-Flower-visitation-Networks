# NETWORK METRICS
# Site-level and Site × Time Period plant-flower visitor networks


# Packages

library(dplyr)
library(tidyr)
library(tibble)
library(bipartite)


# Load data

pollinators <- read.csv(
  "Raw Data/Pollinators.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)

#Although the dataframe and some internal object names use the term
#"pollinator", the observations represent flower visitation rather than
#confirmed pollination.


# Prepare network data

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
    
    Group = trimws(
      as.character(Group)
    ),
    
    Plant.ID = factor(
      Plant.ID
    ),
    
    Site = factor(
      Site
    ),
    
    Time.Period = factor(
      Time.Period
    )
  )


# Plant scientific names

plant_names <- data.frame(
  Plant.ID = factor(
    1:14
  ),
  
  Plant.Species = c(
    "Hypochaeris radicata",   # 1
    "Dianthus deltoides",     # 2
    "Thymus drucei",          # 3
    "Lotus corniculatus",     # 4
    "Silene viscaria",        # 5
    "Galium saxatile",        # 6
    "Hypochaeris radicata",   # 7
    "Thymus drucei",          # 8
    "Lotus corniculatus",     # 9
    "Silene viscaria",        # 10
    "Teucrium scorodonia",    # 11
    "Galium verum",           # 12
    "Galium verum",           # 13
    "Cirsium arvense"         # 14
  ),
  
  stringsAsFactors = FALSE
)


# Prepare species-level interaction data

#Use species-level identity wherever a reliable species name is
#available. Otherwise, retain the broader flower-visitor group.

species_network_data <- network_data %>%
  left_join(
    plant_names,
    by = "Plant.ID"
  ) %>%
  filter(
    !is.na(Plant.Species)
  ) %>%
  mutate(
    
    Flower.Visitor.Label = case_when(
      !is.na(Species) &
        Species != "" &
        Species != "NA" &
        Species != "Unknown" ~
        
        sub(
          "^([A-Za-z])[A-Za-z]+ ",
          "\\1. ",
          Species
        ),
      
      TRUE ~
        as.character(Group)
    ),
    
    #Abbreviate plant scientific names
    
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
    )
  ) %>%
  filter(
    !is.na(Flower.Visitor.Label),
    Flower.Visitor.Label != ""
  )


# Create network matrices

#Function to create a plant × flower-visitor matrix

make_species_network <- function(data) {
  
  #Return an empty matrix if no interactions are present
  
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
    Visits ~
      Plant.Label +
      Flower.Visitor.Label,
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
  
  
  #Remove plants and flower visitors with no interactions
  
  web <- web[
    rowSums(web) > 0,
    colSums(web) > 0,
    drop = FALSE
  ]
  
  
  web
}


# Overall network

overall_species_web <- species_network_data %>%
  make_species_network()


# Site 1

site1_species_web <- species_network_data %>%
  filter(
    Site == "Site 1"
  ) %>%
  make_species_network()


# Site 2

site2_species_web <- species_network_data %>%
  filter(
    Site == "Site 2"
  ) %>%
  make_species_network()


# Inspect site-level matrices

site1_species_web
site2_species_web


# Site × Time Period networks

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


# Metric functions

#Clean a network matrix before calculating metrics

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


# Weighted NODF

#Measure of weighted nestedness

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


# Connectance

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


# Modularity

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


# Site-level network metrics

set.seed(123)

site1_modularity <- calculate_modularity_meta(
  site1_species_web,
  module_reps = 500
)


set.seed(456)

site2_modularity <- calculate_modularity_meta(
  site2_species_web,
  module_reps = 500
)


site_metrics <- tibble(
  Site = c(
    "Site 1",
    "Site 2"
  ),
  
  Plant.Nodes = c(
    nrow(
      site1_species_web
    ),
    nrow(
      site2_species_web
    )
  ),
  
  Flower.Visitor.Nodes = c(
    ncol(
      site1_species_web
    ),
    ncol(
      site2_species_web
    )
  ),
  
  Total.Visits = c(
    sum(
      site1_species_web
    ),
    sum(
      site2_species_web
    )
  ),
  
  Number.of.Links = c(
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
  
  Weighted.NODF = c(
    calculate_weighted_nodf(
      site1_species_web
    ),
    
    calculate_weighted_nodf(
      site2_species_web
    )
  ),
  
  Modularity.Q = c(
    site1_modularity,
    site2_modularity
  )
) %>%
  mutate(
    Connectance = round(
      Connectance,
      3
    ),
    
    Weighted.NODF = round(
      Weighted.NODF,
      2
    ),
    
    Modularity.Q = round(
      Modularity.Q,
      3
    )
  )


site_metrics


# Site × Time Period network metrics

#Modularity contains a stochastic component, so seeds are set
#to make the results reproducible.

set.seed(123)

site1_tp1_modularity <- calculate_modularity_meta(
  site1_tp1_web,
  module_reps = 500
)

site1_tp2_modularity <- calculate_modularity_meta(
  site1_tp2_web,
  module_reps = 500
)

site1_tp3_modularity <- calculate_modularity_meta(
  site1_tp3_web,
  module_reps = 500
)


set.seed(456)

site2_tp1_modularity <- calculate_modularity_meta(
  site2_tp1_web,
  module_reps = 500
)

site2_tp2_modularity <- calculate_modularity_meta(
  site2_tp2_web,
  module_reps = 500
)

site2_tp3_modularity <- calculate_modularity_meta(
  site2_tp3_web,
  module_reps = 500
)


site_time_metrics <- tibble(
  Network = c(
    "Site 1 - Time Period 1",
    "Site 1 - Time Period 2",
    "Site 1 - Time Period 3",
    "Site 2 - Time Period 1",
    "Site 2 - Time Period 2",
    "Site 2 - Time Period 3"
  ),
  
  Plant.Nodes = c(
    nrow(site1_tp1_web),
    nrow(site1_tp2_web),
    nrow(site1_tp3_web),
    nrow(site2_tp1_web),
    nrow(site2_tp2_web),
    nrow(site2_tp3_web)
  ),
  
  Flower.Visitor.Nodes = c(
    ncol(site1_tp1_web),
    ncol(site1_tp2_web),
    ncol(site1_tp3_web),
    ncol(site2_tp1_web),
    ncol(site2_tp2_web),
    ncol(site2_tp3_web)
  ),
  
  Total.Visits = c(
    sum(site1_tp1_web),
    sum(site1_tp2_web),
    sum(site1_tp3_web),
    sum(site2_tp1_web),
    sum(site2_tp2_web),
    sum(site2_tp3_web)
  ),
  
  Number.of.Links = c(
    sum(site1_tp1_web > 0),
    sum(site1_tp2_web > 0),
    sum(site1_tp3_web > 0),
    sum(site2_tp1_web > 0),
    sum(site2_tp2_web > 0),
    sum(site2_tp3_web > 0)
  ),
  
  Connectance = c(
    calculate_connectance(site1_tp1_web),
    calculate_connectance(site1_tp2_web),
    calculate_connectance(site1_tp3_web),
    calculate_connectance(site2_tp1_web),
    calculate_connectance(site2_tp2_web),
    calculate_connectance(site2_tp3_web)
  ),
  
  Weighted.NODF = c(
    calculate_weighted_nodf(site1_tp1_web),
    calculate_weighted_nodf(site1_tp2_web),
    calculate_weighted_nodf(site1_tp3_web),
    calculate_weighted_nodf(site2_tp1_web),
    calculate_weighted_nodf(site2_tp2_web),
    calculate_weighted_nodf(site2_tp3_web)
  ),
  
  Modularity.Q = c(
    site1_tp1_modularity,
    site1_tp2_modularity,
    site1_tp3_modularity,
    site2_tp1_modularity,
    site2_tp2_modularity,
    site2_tp3_modularity
  )
) %>%
  mutate(
    Connectance = round(
      Connectance,
      3
    ),
    
    Weighted.NODF = round(
      Weighted.NODF,
      2
    ),
    
    Modularity.Q = round(
      Modularity.Q,
      3
    )
  )


site_time_metrics


# Null-model functions

#Null-model standardisation is applied to the site-level networks.

#Metrics:
#1. Weighted NODF
#2. Modularity Q
#3. Connectance

#Null models:
#R2d = r2dtable
#VazNull = vaznull

#500 null networks are generated for each null model.


# Standardised effect size

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


# Two-tailed null-model p-value

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


# Weighted NODF null models

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
  
  
  #R2d
  
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
  
  
  #VazNull
  
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
    
    Null.Mean.R2d = mean(
      r2d_values,
      na.rm = TRUE
    ),
    
    Null.SD.R2d = sd(
      r2d_values,
      na.rm = TRUE
    ),
    
    SES.R2d = calculate_ses(
      observed,
      r2d_values
    ),
    
    p.R2d = calculate_two_tailed_p(
      observed,
      r2d_values
    ),
    
    Null.Mean.VazNull = mean(
      vaz_values,
      na.rm = TRUE
    ),
    
    Null.SD.VazNull = sd(
      vaz_values,
      na.rm = TRUE
    ),
    
    SES.VazNull = calculate_ses(
      observed,
      vaz_values
    ),
    
    p.VazNull = calculate_two_tailed_p(
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


# Modularity null models

standardise_modularity <- function(
    web,
    null_reps = 500,
    module_reps = 500,
    seed = 123
) {
  
  web <- clean_network_matrix(
    web
  )
  
  
  observed <- calculate_modularity_meta(
    web,
    module_reps = module_reps
  )
  
  
  #R2d
  
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
      
      calculate_modularity_meta(
        null_web,
        module_reps = module_reps
      )
    },
    
    numeric(1)
  )
  
  
  #VazNull
  
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
      
      calculate_modularity_meta(
        null_web,
        module_reps = module_reps
      )
    },
    
    numeric(1)
  )
  
  
  summary_table <- tibble(
    Metric = "Modularity Q",
    Observed = observed,
    
    Null.Mean.R2d = mean(
      r2d_values,
      na.rm = TRUE
    ),
    
    Null.SD.R2d = sd(
      r2d_values,
      na.rm = TRUE
    ),
    
    SES.R2d = calculate_ses(
      observed,
      r2d_values
    ),
    
    p.R2d = calculate_two_tailed_p(
      observed,
      r2d_values
    ),
    
    Null.Mean.VazNull = mean(
      vaz_values,
      na.rm = TRUE
    ),
    
    Null.SD.VazNull = sd(
      vaz_values,
      na.rm = TRUE
    ),
    
    SES.VazNull = calculate_ses(
      observed,
      vaz_values
    ),
    
    p.VazNull = calculate_two_tailed_p(
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


# Connectance null models

#VazNull preserves connectance, so only R2d is used
#for standardisation of connectance.

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
    
    Null.Mean.R2d = mean(
      r2d_values,
      na.rm = TRUE
    ),
    
    Null.SD.R2d = sd(
      r2d_values,
      na.rm = TRUE
    ),
    
    SES.R2d = calculate_ses(
      observed,
      r2d_values
    ),
    
    p.R2d = calculate_two_tailed_p(
      observed,
      r2d_values
    ),
    
    Null.Mean.VazNull = NA_real_,
    Null.SD.VazNull = NA_real_,
    SES.VazNull = NA_real_,
    p.VazNull = NA_real_
  )
  
  
  list(
    summary = summary_table,
    r2d_values = r2d_values
  )
}


# Run null models

# Site 1

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


# Site 2

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


# Combine null-model results

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


# Rounded reporting table

network_null_results_report <- network_null_results %>%
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

network_null_results_report


# Final outputs

site_metrics

site_time_metrics

network_null_results_report