#Master Script

#Import Data

pollinators <- read.csv(
  "Pollinators.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)

#Load Packages

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(glmmTMB)
library(MASS)
library(car)
library(emmeans)
library(DHARMa)
library(vegan)
library(bipartite)
library(tibble)
library(ggtext)


#Make sure all required columns are in the dataset

required_columns <- c(
  "Pollinator.ID",
  "Date",
  "Species",
  "Genus",
  "Family",
  "Order",
  "Group",
  "Visits",
  "Block",
  "Replication",
  "Time.Period",
  "Site",
  "Plant.ID",
  "Number.of.Flowers",
  "Abundance",
  "Start.Time",
  "Temperature..C.degrees.",
  "Solar.Power..W.m.2.",
  "Humidity....RH.",
  "Wind.Speed..m.s."
)

missing_columns <- setdiff(
  required_columns,
  names(pollinators)
)

if (length(missing_columns) > 0) {
  stop(
    paste(
      "These required columns are missing:",
      paste(missing_columns, collapse = ", ")
    )
  )
}


#Clean data

pollinators <- pollinators %>%
  mutate(
    Order = str_squish(Order),
    Family = str_squish(Family),
    Genus = str_squish(Genus),
    Species = str_squish(Species),
    Group = str_squish(Group),
    
    # Standardise group terminology from other fly to non-syrphid fly (true fly not hoverfly)
    Group = case_when(
      Group == "Other Fly" ~ "Non-syrphid Fly",
      TRUE ~ Group
    )
  )


#Prep variables

pollinators$Group <- factor(trimws(pollinators$Group))
pollinators$Plant.ID <- factor(pollinators$Plant.ID)
pollinators$Time.Period <- factor(pollinators$Time.Period)
pollinators$Site <- factor(pollinators$Site)
pollinators$Block <- factor(pollinators$Block)
pollinators$Replication <- factor(pollinators$Replication)

analysis_data <- droplevels(
  subset(
    pollinators,
    Group != "None"
  )
)


#Change reference to bumblebee

analysis_data$Group <- relevel(
  analysis_data$Group,
  ref = "Bumblebee"
)


#Remove zero-pollinator records
#Keep one row per unique taxon
#Sort

pollinator_list <- pollinators %>%
  filter(Group != "None") %>%
  distinct(
    Order,
    Family,
    Species,
    Group
  ) %>%
  arrange(
    Order,
    Family,
    Species
  )

pollinator_list


survey_data <- pollinators %>%
  mutate(
    # Convert Date from day/month/year format
    Date = as.Date(
      Date,
      format = "%d/%m/%Y"
    ),
    
    # Correct the two confirmed date-entry errors
    Date = case_when(
      Date == as.Date("2026-06-06") ~
        as.Date("2026-06-05"),
      
      Date == as.Date("2026-07-03") ~
        as.Date("2026-07-02"),
      
      TRUE ~ Date
    ),
    
    Time.Period = factor(Time.Period),
    Site = factor(Site),
    Plant.ID = factor(Plant.ID),
    Block = factor(Block),
    Replication = factor(Replication)
  ) %>%
  
  group_by(
    Date,
    Block,
    Replication,
    Time.Period,
    Site,
    Plant.ID,
    Start.Time
  ) %>%
  
  summarise(
    # Sum visits from all pollinator groups recorded during
    # the same 20-minute plant survey
    Visits = sum(
      Visits,
      na.rm = TRUE
    ),
    
    Number.of.Flowers = first(
      Number.of.Flowers
    ),
    
    Abundance = first(
      Abundance
    ),
    
    Minutes = first(
      Minutes
    ),
    
    .groups = "drop"
  ) %>%
  
  mutate(
    Date = droplevels(factor(Date)),
    Time.Period = droplevels(factor(Time.Period)),
    Site = droplevels(factor(Site)),
    Plant.ID = droplevels(factor(Plant.ID))
  )

#Check the survey-level data

head(survey_data)

# Number of individual surveys
nrow(survey_data)

# Distribution of visits
summary(survey_data$Visits)

# Number of zero-visit surveys
sum(survey_data$Visits == 0)

# Confirm there is one row per survey
anyDuplicated(
  survey_data[
    c(
      "Date",
      "Block",
      "Replication",
      "Time.Period",
      "Site",
      "Plant.ID",
      "Start.Time"
    )
  ]
)

# Check corrected survey dates
table(survey_data$Date)

# Check plant sample sizes
table(survey_data$Plant.ID)

# Check plant IDs within sites
table(
  survey_data$Plant.ID,
  survey_data$Site
)
###############################################################
# Questions 1 and 2
# 1. Does pollinator visitation vary throughout the day?
# 2. Does pollinator visitation vary between sites?
###############################################################

#Prepare dataset

# Plant IDs 6 and 13 received no pollinator visits in any survey.
# They are retained in descriptive summaries but excluded from
# inferential count models because they provide no within-plant
# variation in visitation.

survey_model <- survey_data %>%
  filter(
    !(Plant.ID %in% c("6", "13"))
  ) %>%
  mutate(
    Plant.ID = droplevels(Plant.ID),
    Date = droplevels(Date),
    Site = droplevels(Site),
    Time.Period = droplevels(Time.Period)
  )


# Check final analysis structure

nrow(survey_model)

table(survey_model$Date)
table(survey_model$Plant.ID)
table(survey_model$Time.Period)
table(survey_model$Site)


#Block or no block
#Block was determined as the first two weeks and the second two weeks
#This was determined using a negative binomial glm

Model_noBlock <- glm.nb(
  Visits ~
    Time.Period +
    Site +
    Plant.ID,
  data = survey_model
)

Model_withBlock <- glm.nb(
  Visits ~
    Time.Period +
    Site +
    Plant.ID +
    Block,
  data = survey_model
)

anova(
  Model_noBlock,
  Model_withBlock,
  test = "Chisq"
)

#Addition of block was not significant


#Time-period and site model

#Inclusion of date and plant ID as random effects

ModelTimeSite <- glmmTMB(
  Visits ~
    Time.Period +
    Site +
    (1 | Date) +
    (1 | Plant.ID),
  family = nbinom2,
  data = survey_model
)


#Model Diagnostics and checking of the residuals for nbinom2

simulationOutput_TimeSite_nbinom2 <- simulateResiduals(
  fittedModel = ModelTimeSite,
  n = 1000
)

plot(simulationOutput_TimeSite_nbinom2)

Uniformity_TimeSite_nbinom2 <- testUniformity(
  simulationOutput_TimeSite_nbinom2
)

Dispersion_TimeSite_nbinom2 <- testDispersion(
  simulationOutput_TimeSite_nbinom2
)

ZeroInflation_TimeSite_nbinom2 <- testZeroInflation(
  simulationOutput_TimeSite_nbinom2
)

Uniformity_TimeSite_nbinom2
Dispersion_TimeSite_nbinom2
ZeroInflation_TimeSite_nbinom2


#Compare nbinom1 and nbinom2
#nbinom1: variance increases linearly with the mean.
#nbinom2: variance increases quadratically with the mean.

ModelTimeSite_nbinom1 <- update(
  ModelTimeSite,
  family = nbinom1
)

ModelTimeSite_nbinom1$sdr$pdHess

AIC(
  ModelTimeSite,
  ModelTimeSite_nbinom1
)


#Model Diagnostics and checking of the residuals for nbinom1

simulationOutput_TimeSite_nbinom1 <- simulateResiduals(
  fittedModel = ModelTimeSite_nbinom1,
  n = 1000
)

plot(simulationOutput_TimeSite_nbinom1)

Uniformity_TimeSite_nbinom1 <- testUniformity(
  simulationOutput_TimeSite_nbinom1
)

Dispersion_TimeSite_nbinom1 <- testDispersion(
  simulationOutput_TimeSite_nbinom1
)

ZeroInflation_TimeSite_nbinom1 <- testZeroInflation(
  simulationOutput_TimeSite_nbinom1
)

Uniformity_TimeSite_nbinom1
Dispersion_TimeSite_nbinom1
ZeroInflation_TimeSite_nbinom1


#nbinom1 was chosen as it converges cleanly;
#pdHess is TRUE; diagnostics are satisfactory; AIC is meaningfully lower.

#Final model

ModelTimeSite_final <- ModelTimeSite_nbinom1


#Test the overall fixed effects

TimeSite_Anova <- car::Anova(
  ModelTimeSite_final,
  type = 2
)

TimeSite_Anova

# Time.Period:
# Does visitation vary throughout the day?

# Site:
# Does visitation differ between the two sites?


#################################################
#Question 1: Estimated visitation by time period
#################################################

#emmeans tells you what visitation the model predicts for each group
#after accounting for the other variables in the model.

Time_emmeans <- emmeans(
  ModelTimeSite_final,
  pairwise ~ Time.Period,
  type = "response",
  adjust = "tukey"
)

Time_emmeans

# Predicted visits per survey for each period

Time_emmeans$emmeans

# Tukey-adjusted pairwise comparisons

Time_emmeans$contrasts


#################################################
#Question 2: Estimated visitation by site
#################################################

Site_emmeans <- emmeans(
  ModelTimeSite_final,
  pairwise ~ Site,
  type = "response"
)

Site_emmeans

# Predicted visits per survey for each site

Site_emmeans$emmeans

# Site comparison as a response-scale ratio

Site_emmeans$contrasts


#Descriptive Summaries

survey_data %>%
  group_by(Time.Period) %>%
  summarise(
    Surveys = n(),
    Total_Visits = sum(Visits),
    Mean_Visits = mean(Visits),
    Median_Visits = median(Visits),
    Zero_Visit_Surveys = sum(Visits == 0),
    .groups = "drop"
  )

survey_data %>%
  group_by(Site) %>%
  summarise(
    Surveys = n(),
    Total_Visits = sum(Visits),
    Mean_Visits = mean(Visits),
    Median_Visits = median(Visits),
    Zero_Visit_Surveys = sum(Visits == 0),
    .groups = "drop"
  )

survey_data %>%
  group_by(Plant.ID, Site) %>%
  summarise(
    Surveys = n(),
    Total_Visits = sum(Visits),
    Mean_Visits = mean(Visits),
    .groups = "drop"
  )


#Check the random effects and the convergence

VarCorr(ModelTimeSite_final)

ModelTimeSite_final$sdr$pdHess

# TRUE which indicates a positive-definite Hessian and supports
# successful convergence.


#Check the value of date as a random effect

ModelTimeSite_no_date <- glmmTMB(
  Visits ~
    Time.Period +
    Site +
    (1 | Plant.ID),
  family = nbinom1,
  data = survey_model
)

AIC(
  ModelTimeSite_no_date,
  ModelTimeSite_final
)

#AIC was better in the model which included date as a random effect


#Final outputs to be saved

summary(ModelTimeSite_final)

TimeSite_Anova

Time_emmeans$emmeans
Time_emmeans$contrasts

Site_emmeans$emmeans
Site_emmeans$contrasts

VarCorr(ModelTimeSite_final)

AIC(ModelTimeSite_final)

ModelTimeSite_final$sdr$pdHess

plot(simulationOutput_TimeSite_nbinom1)

Uniformity_TimeSite_nbinom1
Dispersion_TimeSite_nbinom1
ZeroInflation_TimeSite_nbinom1

############################
#3. Network analysis metrics
############################
# Overall, site-level and Site × Time Period networks

#Prepare genuine plant-pollinator interactions

network_data <- pollinators %>%
  filter(
    Group != "None",
    Visits > 0,
    !is.na(Plant.ID),
    !is.na(Site),
    !is.na(Time.Period)
  ) %>%
  mutate(
    Species = trimws(Species),
    Group = trimws(Group),
    Plant.ID = factor(Plant.ID),
    Site = factor(Site),
    Time.Period = factor(Time.Period)
  )

# Change to scientific names

plant_names <- data.frame(
  Plant.ID = factor(1:14),
  
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

#Prepare species-level interaction data

# Use species-level identity wherever a reliable species name is
# available. Otherwise, retain the broader pollinator group.

species_network_data <- network_data %>%
  left_join(
    plant_names,
    by = "Plant.ID"
  ) %>%
  filter(
    !is.na(Plant.Species)
  ) %>%
  mutate(
    
    Pollinator.Label = case_when(
      !is.na(Species) &
        Species != "" &
        Species != "NA" &
        Species != "Unknown" ~
        
        sub(
          "^([A-Za-z])[A-Za-z]+ ",
          "\\1. ",
          Species
        ),
      
      TRUE ~ as.character(Group)
    ),
    
    # Abbreviate plant scientific names
    Plant.Label = sub(
      "^([A-Za-z])[A-Za-z]+ ",
      "\\1. ",
      Plant.Species
    ),
    
    Site = factor(
      Site,
      levels = c("1", "2"),
      labels = c("Site 1", "Site 2")
    )
  ) %>%
  filter(
    !is.na(Pollinator.Label),
    Pollinator.Label != ""
  )

#Function to create a plant × pollinator matrix

make_species_network <- function(data) {
  
  # Return an empty matrix if no interactions are present
  
  if (nrow(data) == 0) {
    
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
    data = as.numeric(web_table),
    nrow = nrow(web_table),
    ncol = ncol(web_table),
    dimnames = dimnames(web_table)
  )
  
  # Remove plants and pollinators with no interactions
  
  web <- web[
    rowSums(web) > 0,
    colSums(web) > 0,
    drop = FALSE
  ]
  
  web
}

#Create overall and site-level species networks

overall_species_web <- species_network_data %>%
  make_species_network()


site1_species_web <- species_network_data %>%
  filter(
    Site == "Site 1"
  ) %>%
  make_species_network()


site2_species_web <- species_network_data %>%
  filter(
    Site == "Site 2"
  ) %>%
  make_species_network()


#Inspect the matrices

site1_species_web
site2_species_web

#Site × time period species-level networks

site1_tp1_web <- species_network_data %>%
  filter(
    Site == "Site 1",
    Time.Period == 1
  ) %>%
  make_species_network()


site1_tp2_web <- species_network_data %>%
  filter(
    Site == "Site 1",
    Time.Period == 2
  ) %>%
  make_species_network()


site1_tp3_web <- species_network_data %>%
  filter(
    Site == "Site 1",
    Time.Period == 3
  ) %>%
  make_species_network()


site2_tp1_web <- species_network_data %>%
  filter(
    Site == "Site 2",
    Time.Period == 1
  ) %>%
  make_species_network()


site2_tp2_web <- species_network_data %>%
  filter(
    Site == "Site 2",
    Time.Period == 2
  ) %>%
  make_species_network()


site2_tp3_web <- species_network_data %>%
  filter(
    Site == "Site 2",
    Time.Period == 3
  ) %>%
  make_species_network()

#Network metric functions

#Clean the network matrix

clean_network_matrix <- function(web) {
  
  web <- as.matrix(web)
  
  storage.mode(web) <- "numeric"
  
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
#Measure of weighted nestedness

calculate_weighted_nodf <- function(web) {
  
  web <- clean_network_matrix(web)
  
  as.numeric(
    bipartite::networklevel(
      web,
      index = "weighted NODF"
    )
  )
}

#Connectance

calculate_connectance <- function(web) {
  
  web <- clean_network_matrix(web)
  
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
  
  web <- clean_network_matrix(web)
  
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

site_metrics <- tibble(
  Network = c(
    "Site 1",
    "Site 2"
  ),
  
  Plant_Nodes = c(
    nrow(site1_species_web),
    nrow(site2_species_web)
  ),
  
  Pollinator_Nodes = c(
    ncol(site1_species_web),
    ncol(site2_species_web)
  ),
  
  Total_Visits = c(
    sum(site1_species_web),
    sum(site2_species_web)
  ),
  
  Number_of_Links = c(
    sum(site1_species_web > 0),
    sum(site2_species_web > 0)
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
    calculate_modularity_meta(
      site1_species_web,
      module_reps = 100
    ),
    calculate_modularity_meta(
      site2_species_web,
      module_reps = 100
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


site_metrics

#Site × time period network metrics

site_time_metrics <- tibble(
  Network = c(
    "Site 1 - Time Period 1",
    "Site 1 - Time Period 2",
    "Site 1 - Time Period 3",
    "Site 2 - Time Period 1",
    "Site 2 - Time Period 2",
    "Site 2 - Time Period 3"
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

#Species-level module plots by site
#Scientific plant and pollinator names

#Detect modules

set.seed(123)

site1_modules <- bipartite::DIRT_LPA_wb_plus(
  site1_species_web,
  reps = 100
)


set.seed(123)

site2_modules <- bipartite::DIRT_LPA_wb_plus(
  site2_species_web,
  reps = 100
)

#Convert matrices into module webs

site1_module_web <- bipartite::convert2moduleWeb(
  site1_species_web,
  site1_modules
)


site2_module_web <- bipartite::convert2moduleWeb(
  site2_species_web,
  site2_modules
)

#Site 1 module plot

par(
  mar = c(13, 12, 8, 2),
  xpd = NA
)

bipartite::plotModuleWeb(
  site1_module_web,
  plotModules = TRUE,
  weighted = FALSE,
  displayAlabels = TRUE,
  displayBlabels = TRUE,
  labsize = 0.8,
  xlabel = "",
  ylabel = ""
)

mtext(
  "Site 1",
  side = 3,
  line = 2,
  font = 2,
  cex = 1.4
)

mtext(
  "Pollinator taxon",
  side = 1,
  line = 11,
  font = 2,
  cex = 1.1
)

mtext(
  "Plant species",
  side = 2,
  line = 10,
  font = 2,
  cex = 1.1
)

#Site 2 module plot

par(
  mar = c(13, 12, 8, 2),
  xpd = NA
)

bipartite::plotModuleWeb(
  site2_module_web,
  plotModules = TRUE,
  weighted = FALSE,
  displayAlabels = TRUE,
  displayBlabels = TRUE,
  labsize = 0.8,
  xlabel = "",
  ylabel = ""
)

mtext(
  "Site 2",
  side = 3,
  line = 2,
  font = 2,
  cex = 1.4
)

mtext(
  "Pollinator taxon",
  side = 1,
  line = 11,
  font = 2,
  cex = 1.1
)

mtext(
  "Plant species",
  side = 2,
  line = 10,
  font = 2,
  cex = 1.1
)

# Weighted nestedness by site

site1_weighted_nodf <- bipartite::networklevel(
  site1_species_web,
  index = "weighted NODF"
)

site2_weighted_nodf <- bipartite::networklevel(
  site2_species_web,
  index = "weighted NODF"
)


site_nestedness <- data.frame(
  Site = c(
    "Site 1",
    "Site 2"
  ),
  
  Weighted_NODF = round(
    c(
      as.numeric(site1_weighted_nodf),
      as.numeric(site2_weighted_nodf)
    ),
    2
  )
)


site_nestedness

#Prepare nestedness heatmap data

prepare_nestedness_data <- function(web) {
  
  #Sort plants and pollinators by decreasing interaction totals
  
  sorted_web <- bipartite::sortweb(
    web,
    sort.order = "dec"
  )
  
  #Remove empty rows and columns
  
  sorted_web <- sorted_web[
    rowSums(sorted_web) > 0,
    colSums(sorted_web) > 0,
    drop = FALSE
  ]
  
  #Convert matrix to long format
  
  nested_df <- as.data.frame(
    as.table(sorted_web)
  )
  
  names(nested_df) <- c(
    "Plant",
    "Pollinator",
    "Visits"
  )
  
  #Preserve the order returned by sortweb()
  
  nested_df <- nested_df %>%
    mutate(
      Plant = factor(
        Plant,
        levels = rev(
          rownames(sorted_web)
        )
      ),
      
      Pollinator = factor(
        Pollinator,
        levels = colnames(sorted_web)
      )
    )
  
  nested_df
}


#Actually create the plotting data

site1_nested_df <- prepare_nestedness_data(
  site1_species_web
)

site2_nested_df <- prepare_nestedness_data(
  site2_species_web
)


# Highlight restored plants

restored_labels <- c(
  "D. deltoides",
  "S. viscaria"
)


plant_label_function <- function(x) {
  
  ifelse(
    x %in% restored_labels,
    
    paste0(
      "<span style='color:#C2185B;'><b><i>",
      x,
      "</i></b></span>"
    ),
    
    paste0(
      "<i>",
      x,
      "</i>"
    )
  )
}


# Site 1 nestedness plot

site1_nestedness_plot <- ggplot(
  site1_nested_df,
  aes(
    x = Pollinator,
    y = Plant,
    fill = Visits
  )
) +
  
  geom_tile(
    colour = "white",
    linewidth = 0.4
  ) +
  
  scale_fill_gradient(
    low = "mintcream",
    high = "seagreen4",
    trans = "sqrt",
    name = "Number of visits"
  ) +
  
  scale_y_discrete(
    labels = plant_label_function
  ) +
  
  labs(
    title = paste0(
      "Site 1: weighted NODF = ",
      round(
        as.numeric(site1_weighted_nodf),
        2
      )
    ),
    x = "Pollinator taxon",
    y = "Plant species"
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    
    axis.text.y = ggtext::element_markdown(),
    
    axis.ticks = element_blank(),
    
    plot.margin = margin(
      t = 15,
      r = 15,
      b = 15,
      l = 15
    )
  )


site1_nestedness_plot

#Site 2 nestedness plot

site2_nestedness_plot <- ggplot(
  site2_nested_df,
  aes(
    x = Pollinator,
    y = Plant,
    fill = Visits
  )
) +
  
  geom_tile(
    colour = "white",
    linewidth = 0.4
  ) +
  
  scale_fill_gradient(
    low = "lavender",
    high = "mediumorchid3",
    trans = "sqrt",
    name = "Number of visits"
  ) +
  
  scale_y_discrete(
    labels = plant_label_function
  ) +
  
  labs(
    title = paste0(
      "Site 2: weighted NODF = ",
      round(
        as.numeric(site2_weighted_nodf),
        2
      )
    ),
    x = "Pollinator taxon",
    y = "Plant species"
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    
    axis.text.y = ggtext::element_markdown(),
    
    axis.ticks = element_blank(),
    
    plot.margin = margin(
      t = 15,
      r = 15,
      b = 15,
      l = 15
    )
  )


site2_nestedness_plot


# ============================================================
# Null-model standardisation of network metrics
# Species-level weighted networks by site
#
# Metrics:
# 1. Weighted NODF
# 2. Modularity Q
# 3. Connectance
#
# Null models:
# R2d = r2dtable
# VazNull = vaznull
#
# Replicates:
# 500 null networks
# ============================================================

#Helper functions

calculate_ses <- function(
    observed,
    null_values
) {
  
  null_values <- null_values[
    is.finite(null_values)
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
    observed - mean(null_values)
  ) / null_sd
}


calculate_two_tailed_p <- function(
    observed,
    null_values
) {
  
  null_values <- null_values[
    is.finite(null_values)
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
  ) / (
    length(null_values) + 1
  )
  
  upper_p <- (
    sum(
      null_values >= observed
    ) + 1
  ) / (
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

# Weighted NODF null-model standardisation

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
  
  observed <- calculate_modularity_meta(
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
      
      calculate_modularity_meta(
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
#VazNull preserves connectance, so only R2d is used for SES.

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

#Final null-model analysis
#500 null networks were generated for each null model.
#Modularity was calculated using repeated module detection.

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

#Rounded reporting table

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

#Inspect null distributions

summary(
  site1_modularity_results$r2d_values
)

summary(
  site1_modularity_results$vaz_values
)

summary(
  site2_modularity_results$r2d_values
)

summary(
  site2_modularity_results$vaz_values
)


hist(
  site1_modularity_results$r2d_values,
  main = "Site 1: R2d null modularity",
  xlab = "Modularity Q"
)

abline(
  v = site1_modularity_results$summary$Observed,
  lwd = 2
)


hist(
  site1_modularity_results$vaz_values,
  main = "Site 1: VazNull modularity",
  xlab = "Modularity Q"
)

abline(
  v = site1_modularity_results$summary$Observed,
  lwd = 2
)


hist(
  site2_modularity_results$r2d_values,
  main = "Site 2: R2d null modularity",
  xlab = "Modularity Q"
)

abline(
  v = site2_modularity_results$summary$Observed,
  lwd = 2
)


hist(
  site2_modularity_results$vaz_values,
  main = "Site 2: VazNull modularity",
  xlab = "Modularity Q"
)

abline(
  v = site2_modularity_results$summary$Observed,
  lwd = 2
)

#####################################
#4. Plant-pollinator Network Figures
#####################################

library(dplyr)
library(tidyr)
library(tibble)


#Plant species and categories

plant_key <- c(
  "1"  = "Hypochaeris radicata",
  "2"  = "Dianthus deltoides",
  "3"  = "Thymus drucei",
  "4"  = "Lotus corniculatus",
  "5"  = "Silene viscaria",
  "6"  = "Galium saxatile",
  "7"  = "Hypochaeris radicata",
  "8"  = "Thymus drucei",
  "9"  = "Lotus corniculatus",
  "10" = "Silene viscaria",
  "11" = "Teucrium scorodonia",
  "12" = "Galium verum",
  "13" = "Galium verum",
  "14" = "Cirsium arvense"
)

plant_order <- c(
  "Hypochaeris radicata",
  "Dianthus deltoides",
  "Thymus drucei",
  "Lotus corniculatus",
  "Silene viscaria",
  "Galium saxatile",
  "Teucrium scorodonia",
  "Galium verum",
  "Cirsium arvense"
)

restored_plants <- c(
  "Dianthus deltoides",
  "Silene viscaria"
)

restored_colour <- "#E7298A"
other_plant_colour <- "black"

time_period_order <- c(
  "1",
  "2",
  "3"
)


#Prepare network data

clean_network_text <- function(x) {
  
  x <- as.character(x)
  
  x <- gsub(
    "\u00A0",
    " ",
    x,
    fixed = TRUE
  )
  
  x <- gsub(
    "[[:space:]]+",
    " ",
    x
  )
  
  trimws(x)
}


network_figure_data <- pollinators %>%
  mutate(
    Species = clean_network_text(Species),
    Genus = clean_network_text(Genus),
    Family = clean_network_text(Family),
    Order = clean_network_text(Order),
    Group = clean_network_text(Group),
    
    Plant.ID = as.character(
      Plant.ID
    ),
    
    Site = as.character(
      Site
    ),
    
    Time.Period = as.character(
      Time.Period
    ),
    
    Plant.Species = unname(
      plant_key[
        Plant.ID
      ]
    ),
    
    Family.Network = case_when(
      !is.na(Family) &
        Family != "" ~
        Family,
      
      !is.na(Order) &
        Order != "" ~
        paste0(
          Order,
          " unidentified"
        ),
      
      TRUE ~
        "Unidentified"
    ),
    
    Pollinator = case_when(
      !is.na(Species) &
        Species != "" ~
        Species,
      
      !is.na(Genus) &
        Genus != "" ~
        paste0(
          Genus,
          " sp."
        ),
      
      !is.na(Family) &
        Family != "" ~
        paste0(
          Family,
          " family"
        ),
      
      !is.na(Order) &
        Order != "" ~
        paste0(
          Order,
          " order"
        ),
      
      TRUE ~
        NA_character_
    ),
    
    Visit.Rate = case_when(
      !is.na(Visits) &
        !is.na(Number.of.Flowers) &
        Number.of.Flowers > 0 ~
        Visits / Number.of.Flowers,
      
      TRUE ~
        NA_real_
    )
  )


network_interactions <- network_figure_data %>%
  filter(
    Group != "None",
    Visits > 0,
    !is.na(Visit.Rate),
    Visit.Rate > 0,
    !is.na(Plant.Species),
    Plant.Species != "",
    !is.na(Pollinator),
    Pollinator != "",
    !is.na(Site),
    Site != "",
    !is.na(Time.Period),
    Time.Period != ""
  )


network_flowers <- network_figure_data %>%
  filter(
    !is.na(Plant.Species),
    Plant.Species != "",
    !is.na(Abundance),
    Abundance > 0,
    !is.na(Site),
    Site != "",
    !is.na(Time.Period),
    Time.Period != ""
  ) %>%
  distinct(
    Date,
    Site,
    Block,
    Replication,
    Time.Period,
    Plant.ID,
    Start.Time,
    .keep_all = TRUE
  )


#Pollinator colours

pollinator_families <- sort(
  unique(
    network_interactions$Family.Network
  )
)

family_colours <- grDevices::hcl.colors(
  n = length(
    pollinator_families
  ),
  palette = "Dark 3"
)

names(
  family_colours
) <- pollinator_families


#Labels

abbreviate_pollinator <- function(x) {
  
  vapply(
    as.character(x),
    
    function(label) {
      
      if (
        grepl(
          " family$| order$| unidentified$| sp\\.$",
          label,
          ignore.case = TRUE
        )
      ) {
        
        return(
          label
        )
      }
      
      words <- strsplit(
        label,
        "[[:space:]]+"
      )[[1]]
      
      if (
        length(words) >= 2
      ) {
        
        return(
          paste0(
            substr(
              words[1],
              1,
              1
            ),
            ". ",
            paste(
              words[-1],
              collapse = " "
            )
          )
        )
      }
      
      label
    },
    
    character(1)
  )
}


abbreviate_plant <- function(x) {
  
  vapply(
    as.character(x),
    
    function(label) {
      
      words <- strsplit(
        label,
        "[[:space:]]+"
      )[[1]]
      
      if (
        length(words) >= 2
      ) {
        
        return(
          paste0(
            substr(
              words[1],
              1,
              1
            ),
            ". ",
            paste(
              words[-1],
              collapse = " "
            )
          )
        )
      }
      
      label
    },
    
    character(1)
  )
}


#Bar positions

make_bar_positions <- function(
    labels,
    weights,
    left = 0.04,
    right = 0.91,
    gap = 0.012
) {
  
  if (
    length(labels) == 0
  ) {
    
    return(
      tibble(
        Label = character(),
        Start = numeric(),
        End = numeric(),
        Centre = numeric()
      )
    )
  }
  
  weights[
    !is.finite(weights) |
      weights <= 0
  ] <- 0
  
  transformed <- sqrt(
    weights
  )
  
  if (
    max(
      transformed,
      na.rm = TRUE
    ) == 0
  ) {
    
    transformed <- rep(
      1,
      length(labels)
    )
  }
  
  minimum_width <- max(
    transformed
  ) * 0.05
  
  transformed <- pmax(
    transformed,
    minimum_width
  )
  
  available_width <-
    right -
    left -
    gap * (
      length(labels) - 1
    )
  
  widths <-
    transformed /
    sum(transformed) *
    available_width
  
  starts <- cumsum(
    c(
      left,
      head(
        widths,
        -1
      ) +
        gap
    )
  )
  
  ends <-
    starts +
    widths
  
  tibble(
    Label = labels,
    Start = starts,
    End = ends,
    Centre = (
      starts +
        ends
    ) / 2
  )
}


#Draw network panel

draw_network_panel <- function(
    site,
    time_period = NULL,
    panel_letter,
    panel_title,
    show_pollinator_labels = FALSE
) {
  
  site <- as.character(
    site
  )
  
  panel_interactions <- network_interactions %>%
    filter(
      Site == site
    )
  
  panel_flowers <- network_flowers %>%
    filter(
      Site == site
    )
  
  if (
    !is.null(
      time_period
    )
  ) {
    
    time_period <- as.character(
      time_period
    )
    
    panel_interactions <- panel_interactions %>%
      filter(
        Time.Period ==
          time_period
      )
    
    panel_flowers <- panel_flowers %>%
      filter(
        Time.Period ==
          time_period
      )
  }
  
  plot(
    NA,
    xlim = c(
      -0.08,
      1.30
    ),
    ylim = c(
      -0.28,
      1.30
    ),
    type = "n",
    axes = FALSE,
    xlab = "",
    ylab = "",
    xaxs = "i",
    yaxs = "i"
  )
  
  text(
    x = -0.04,
    y = 1.20,
    labels = paste0(
      panel_letter,
      ") ",
      panel_title
    ),
    adj = c(
      0,
      0.5
    ),
    cex = 1.15,
    font = 2
  )
  
  floral_summary <- panel_flowers %>%
    group_by(
      Plant.Species
    ) %>%
    summarise(
      Floral.Abundance = mean(
        Abundance,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    filter(
      is.finite(
        Floral.Abundance
      ),
      Floral.Abundance > 0
    )
  
  plants_with_visits <- unique(
    panel_interactions$Plant.Species
  )
  
  plants_present <- union(
    floral_summary$Plant.Species,
    plants_with_visits
  )
  
  plants_present <- plant_order[
    plant_order %in%
      plants_present
  ]
  
  plant_summary <- tibble(
    Plant.Species =
      plants_present
  ) %>%
    left_join(
      floral_summary,
      by = "Plant.Species"
    ) %>%
    mutate(
      Floral.Abundance = replace_na(
        Floral.Abundance,
        1
      )
    )
  
  plant_positions <- make_bar_positions(
    labels =
      plant_summary$Plant.Species,
    
    weights =
      plant_summary$Floral.Abundance,
    
    left =
      0.04,
    
    right =
      0.91,
    
    gap =
      0.015
  ) %>%
    rename(
      Plant.Species =
        Label
    )
  
  pollinator_summary <- panel_interactions %>%
    group_by(
      Pollinator
    ) %>%
    summarise(
      Total.Rate = sum(
        Visit.Rate,
        na.rm = TRUE
      ),
      
      Family.Network = first(
        Family.Network
      ),
      
      .groups = "drop"
    ) %>%
    arrange(
      desc(
        Total.Rate
      ),
      Pollinator
    )
  
  pollinator_positions <- make_bar_positions(
    labels =
      pollinator_summary$Pollinator,
    
    weights =
      pollinator_summary$Total.Rate,
    
    left =
      0.04,
    
    right =
      0.91,
    
    gap =
      0.008
  ) %>%
    rename(
      Pollinator =
        Label
    ) %>%
    left_join(
      pollinator_summary %>%
        select(
          Pollinator,
          Family.Network
        ),
      by = "Pollinator"
    )
  
  links <- panel_interactions %>%
    group_by(
      Plant.Species,
      Pollinator,
      Family.Network
    ) %>%
    summarise(
      Weight = sum(
        Visit.Rate,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    left_join(
      plant_positions %>%
        select(
          Plant.Species,
          Plant.X =
            Centre
        ),
      by = "Plant.Species"
    ) %>%
    left_join(
      pollinator_positions %>%
        select(
          Pollinator,
          Pollinator.X =
            Centre
        ),
      by = "Pollinator"
    ) %>%
    filter(
      !is.na(
        Plant.X
      ),
      !is.na(
        Pollinator.X
      )
    )
  
  if (
    nrow(
      links
    ) > 0
  ) {
    
    transformed_links <- sqrt(
      links$Weight
    )
    
    if (
      max(
        transformed_links
      ) ==
      min(
        transformed_links
      )
    ) {
      
      link_widths <- rep(
        2,
        nrow(
          links
        )
      )
      
    } else {
      
      link_widths <- 0.8 +
        3.2 *
        (
          transformed_links -
            min(
              transformed_links
            )
        ) /
        (
          max(
            transformed_links
          ) -
            min(
              transformed_links
            )
        )
    }
    
    for (
      i in seq_len(
        nrow(
          links
        )
      )
    ) {
      
      link_colour <- family_colours[
        links$Family.Network[i]
      ]
      
      if (
        is.na(
          link_colour
        )
      ) {
        
        link_colour <- "grey50"
      }
      
      segments(
        x0 =
          links$Plant.X[i],
        
        y0 =
          0.22,
        
        x1 =
          links$Pollinator.X[i],
        
        y1 =
          0.82,
        
        col =
          grDevices::adjustcolor(
            link_colour,
            alpha.f = 0.42
          ),
        
        lwd =
          link_widths[i],
        
        lend =
          "butt"
      )
    }
  }
  
  if (
    nrow(
      plant_positions
    ) > 0
  ) {
    
    plant_colours <- ifelse(
      plant_positions$Plant.Species %in%
        restored_plants,
      
      restored_colour,
      
      other_plant_colour
    )
    
    rect(
      xleft =
        plant_positions$Start,
      
      ybottom =
        0.14,
      
      xright =
        plant_positions$End,
      
      ytop =
        0.22,
      
      col =
        plant_colours,
      
      border =
        plant_colours
    )
    
    text(
      x =
        plant_positions$Centre,
      
      y =
        0.12,
      
      labels =
        abbreviate_plant(
          plant_positions$Plant.Species
        ),
      
      srt =
        45,
      
      adj = c(
        1,
        0.5
      ),
      
      cex =
        0.72,
      
      font =
        3
    )
  }
  
  if (
    nrow(
      pollinator_positions
    ) > 0
  ) {
    
    pollinator_colours <-
      family_colours[
        pollinator_positions$Family.Network
      ]
    
    pollinator_colours[
      is.na(
        pollinator_colours
      )
    ] <- "grey50"
    
    rect(
      xleft =
        pollinator_positions$Start,
      
      ybottom =
        0.82,
      
      xright =
        pollinator_positions$End,
      
      ytop =
        0.91,
      
      col =
        pollinator_colours,
      
      border =
        pollinator_colours
    )
    
    if (
      show_pollinator_labels
    ) {
      
      text(
        x =
          pollinator_positions$Centre,
        
        y =
          0.93,
        
        labels =
          abbreviate_pollinator(
            pollinator_positions$Pollinator
          ),
        
        srt =
          90,
        
        adj = c(
          0,
          0.5
        ),
        
        cex =
          0.63,
        
        col =
          pollinator_colours
      )
    }
  }
  
  panel_total_visits <- sum(
    panel_interactions$Visits,
    na.rm = TRUE
  )
  
  panel_plant_richness <- length(
    plants_present
  )
  
  panel_pollinator_richness <- n_distinct(
    panel_interactions$Pollinator,
    na.rm = TRUE
  )
  
  panel_links <- nrow(
    links
  )
  
  text(
    x = 0.98,
    y = 0.68,
    
    labels = paste0(
      "Visits: ",
      panel_total_visits,
      "\nPlants: ",
      panel_plant_richness,
      "\nPollinators: ",
      panel_pollinator_richness,
      "\nLinks: ",
      panel_links
    ),
    
    adj = c(
      0,
      1
    ),
    
    cex =
      0.72
  )
  
  if (
    nrow(
      panel_interactions
    ) == 0
  ) {
    
    text(
      x = 0.47,
      y = 0.52,
      labels =
        "No pollinator visits recorded",
      cex = 0.85
    )
  }
  
  invisible(
    NULL
  )
}


#Legend

draw_network_legend <- function(
    site
) {
  
  site_families <- network_interactions %>%
    filter(
      Site ==
        as.character(
          site
        )
    ) %>%
    distinct(
      Family.Network
    ) %>%
    arrange(
      Family.Network
    ) %>%
    pull(
      Family.Network
    )
  
  plot(
    NA,
    xlim = c(
      0,
      1
    ),
    ylim = c(
      0,
      1
    ),
    type = "n",
    axes = FALSE,
    xlab = "",
    ylab = ""
  )
  
  text(
    x = 0.5,
    y = 0.88,
    labels =
      "Pollinator family",
    font = 2,
    cex = 1
  )
  
  if (
    length(
      site_families
    ) > 0
  ) {
    
    columns <- min(
      4,
      length(
        site_families
      )
    )
    
    rows <- ceiling(
      length(
        site_families
      ) /
        columns
    )
    
    x_positions <- rep(
      seq(
        0.05,
        0.78,
        length.out =
          columns
      ),
      times =
        rows
    )[
      seq_along(
        site_families
      )
    ]
    
    y_positions <- rep(
      seq(
        0.66,
        0.46,
        length.out =
          rows
      ),
      each =
        columns
    )[
      seq_along(
        site_families
      )
    ]
    
    rect(
      xleft =
        x_positions,
      
      ybottom =
        y_positions - 0.03,
      
      xright =
        x_positions + 0.025,
      
      ytop =
        y_positions + 0.03,
      
      col =
        family_colours[
          site_families
        ],
      
      border =
        family_colours[
          site_families
        ]
    )
    
    text(
      x =
        x_positions + 0.035,
      
      y =
        y_positions,
      
      labels =
        site_families,
      
      adj = c(
        0,
        0.5
      ),
      
      cex =
        0.76
    )
  }
  
  text(
    x = 0.5,
    y = 0.27,
    labels =
      "Plant category",
    font = 2,
    cex = 1
  )
  
  rect(
    xleft = c(
      0.28,
      0.60
    ),
    
    ybottom =
      0.10,
    
    xright = c(
      0.305,
      0.625
    ),
    
    ytop =
      0.17,
    
    col = c(
      other_plant_colour,
      restored_colour
    ),
    
    border = c(
      other_plant_colour,
      restored_colour
    )
  )
  
  text(
    x = c(
      0.315,
      0.635
    ),
    
    y =
      0.135,
    
    labels = c(
      "Co-flowering plants",
      "Restored plants"
    ),
    
    adj = c(
      0,
      0.5
    ),
    
    cex =
      0.80
  )
}


#Complete page for each site

plot_site_network <- function(
    site
) {
  
  site <- as.character(
    site
  )
  
  layout(
    matrix(
      1:5,
      nrow = 5,
      ncol = 1
    ),
    heights = c(
      1.30,
      1,
      1,
      1,
      0.55
    )
  )
  
  par(
    oma = c(
      0,
      0,
      2.5,
      0
    )
  )
  
  par(
    mar = c(
      0.4,
      0.4,
      0.4,
      0.4
    )
  )
  
  draw_network_panel(
    site =
      site,
    
    panel_letter =
      "a",
    
    panel_title =
      "Overall",
    
    show_pollinator_labels =
      TRUE
  )
  
  draw_network_panel(
    site =
      site,
    
    time_period =
      "1",
    
    panel_letter =
      "b",
    
    panel_title =
      "Time period 1"
  )
  
  draw_network_panel(
    site =
      site,
    
    time_period =
      "2",
    
    panel_letter =
      "c",
    
    panel_title =
      "Time period 2"
  )
  
  draw_network_panel(
    site =
      site,
    
    time_period =
      "3",
    
    panel_letter =
      "d",
    
    panel_title =
      "Time period 3"
  )
  
  draw_network_legend(
    site
  )
  
  mtext(
    paste(
      "Site",
      site
    ),
    side = 3,
    outer = TRUE,
    line = 0.5,
    cex = 1.7,
    font = 2
  )
  
  invisible(
    NULL
  )
}


#Check network data

table(
  network_interactions$Site,
  network_interactions$Time.Period
)


#Preview Site 1

plot_site_network(
  "1"
)


#Save Site 1 and Site 2 on separate pages

pdf(
  "Plant_pollinator_networks_by_site.pdf",
  width = 13,
  height = 16,
  family = "sans",
  onefile = TRUE
)

plot_site_network(
  "1"
)

plot_site_network(
  "2"
)

dev.list()
graphics.off()

#Preview Site 1

plot_site_network(
  "1"
)


#Preview Site 2

plot_site_network(
  "2"
)


################################################################################
#5. Does the effect of temperature on pollinator visitation vary among
#pollinator groups?
#Negative-binomial GLMM with survey Date as a random intercept
################################################################################


#Prep data

analysis_data_temp <- analysis_data %>%
  filter(
    !(Group %in% c("Mining Bee", "Wasp")),
    !is.na(Visits),
    !is.na(Group),
    !is.na(Temperature..C.degrees.),
    !is.na(Site),
    !is.na(Time.Period),
    !is.na(Plant.ID),
    !is.na(Date)
  ) %>%
  mutate(
    Group = droplevels(factor(Group)),
    Site = factor(Site),
    Time.Period = factor(Time.Period),
    Plant.ID = droplevels(factor(Plant.ID)),
    Date = factor(Date)
  )


#Check group and date sample sizes

table(analysis_data_temp$Group)
table(analysis_data_temp$Date)


#Mean-centre temperature

mean_temperature <- mean(
  analysis_data_temp$Temperature..C.degrees.,
  na.rm = TRUE
)

analysis_data_temp <- analysis_data_temp %>%
  mutate(
    Temperature_c =
      Temperature..C.degrees. - mean_temperature
  )


#Linear temperature x pollinator-group model

ModelTempGroup_linear <- glmmTMB(
  Visits ~
    Temperature_c * Group +
    Site +
    Time.Period +
    Plant.ID +
    (1 | Date),
  family = nbinom2,
  data = analysis_data_temp
)

summary(ModelTempGroup_linear)

TempGroup_linear_Anova <- car::Anova(
  ModelTempGroup_linear,
  type = 2
)

TempGroup_linear_Anova


#Prep data for quadratic analysis

#Sweat bees are excluded from the quadratic analysis because
#their observations span too few distinct temperatures to
#estimate a group-specific quadratic curve reliably.

analysis_data_quad <- analysis_data_temp %>%
  filter(
    Group != "Sweat Bee"
  ) %>%
  droplevels()

table(analysis_data_quad$Group)
table(analysis_data_quad$Date)


#Re-centre temperature using this reduced dataset

mean_temperature_quad <- mean(
  analysis_data_quad$Temperature..C.degrees.,
  na.rm = TRUE
)

analysis_data_quad <- analysis_data_quad %>%
  mutate(
    Temperature_c =
      Temperature..C.degrees. - mean_temperature_quad
  )


#Linear model on the same quad-analysis dataset

ModelTempGroup_linear_quad_data <- glmmTMB(
  Visits ~
    Temperature_c * Group +
    Site +
    Time.Period +
    Plant.ID +
    (1 | Date),
  family = nbinom2,
  data = analysis_data_quad
)


#Full group-specific quadratic model

ModelTempGroup_quad <- glmmTMB(
  Visits ~
    (
      Temperature_c +
        I(Temperature_c^2)
    ) * Group +
    Site +
    Time.Period +
    Plant.ID +
    (1 | Date),
  family = nbinom2,
  data = analysis_data_quad
)

summary(ModelTempGroup_quad)

TempGroup_quad_Anova <- car::Anova(
  ModelTempGroup_quad,
  type = 2
)

TempGroup_quad_Anova

#Temperature x Group was significant
#Temperature^2 x Group was not significant


#Does a group-specific quadratic response improve fit?

anova(
  ModelTempGroup_linear_quad_data,
  ModelTempGroup_quad
)

AIC(
  ModelTempGroup_linear_quad_data,
  ModelTempGroup_quad
)

#quadratic model did not significantly improve fit
#AIC differed by <2 so the models had similar support


#Simplified Model
#Shared quadratic temperature response, but group-specific
#linear temperature slopes.

ModelTempGroup_final <- glmmTMB(
  Visits ~
    Temperature_c +
    I(Temperature_c^2) +
    Group +
    Temperature_c:Group +
    Site +
    Time.Period +
    Plant.ID +
    (1 | Date),
  family = nbinom2,
  data = analysis_data_quad
)

summary(ModelTempGroup_final)

TempGroup_final_Anova <- car::Anova(
  ModelTempGroup_final,
  type = 2
)

TempGroup_final_Anova

#Group, Time Period, Plant ID and Temperature x Group were significant
#Temperature^2 was not significant at p < 0.05


#Compare the simplified model with the full quadratic interaction

anova(
  ModelTempGroup_final,
  ModelTempGroup_quad
)

AIC(
  ModelTempGroup_final,
  ModelTempGroup_quad
)

#Final model had a slightly lower AIC
#Full Temperature^2 x Group interaction did not significantly improve fit
#Simplified model was retained


#Test the overall linear temperature x group interaction

#This reduced model removes group-specific linear temperature slopes
#while retaining the shared quadratic temperature response.

ModelTempGroup_no_linear_interaction <- glmmTMB(
  Visits ~
    Temperature_c +
    I(Temperature_c^2) +
    Group +
    Site +
    Time.Period +
    Plant.ID +
    (1 | Date),
  family = nbinom2,
  data = analysis_data_quad
)

anova(
  ModelTempGroup_no_linear_interaction,
  ModelTempGroup_final
)

AIC(
  ModelTempGroup_no_linear_interaction,
  ModelTempGroup_final
)

#A significant likelihood-ratio test means that the linear
#temperature effect differs among pollinator groups.

#Temperature x Group significantly improved the model
#so the effect of temperature differed among pollinator groups


#Check random date effect

VarCorr(ModelTempGroup_final)

ModelTempGroup_final$sdr$pdHess

#pdHess was TRUE


ModelTempGroup_no_date <- glmmTMB(
  Visits ~
    Temperature_c +
    I(Temperature_c^2) +
    Group +
    Temperature_c:Group +
    Site +
    Time.Period +
    Plant.ID,
  family = nbinom2,
  data = analysis_data_quad
)

AIC(
  ModelTempGroup_no_date,
  ModelTempGroup_final
)

#Slightly better with date


#Model diagnostics

TempGroup_residuals <- simulateResiduals(
  fittedModel = ModelTempGroup_final,
  n = 1000
)

plot(TempGroup_residuals)

TempGroup_uniformity <- testUniformity(
  TempGroup_residuals
)

TempGroup_dispersion <- testDispersion(
  TempGroup_residuals
)

TempGroup_zero_inflation <- testZeroInflation(
  TempGroup_residuals
)

TempGroup_uniformity
TempGroup_dispersion
TempGroup_zero_inflation

#Uniformity test was not significant


#Check nbinom1 vs nbinom2

ModelTempGroup_final_nbinom1 <- update(
  ModelTempGroup_final,
  family = nbinom1
)

ModelTempGroup_final_nbinom1$sdr$pdHess

AIC(
  ModelTempGroup_final,
  ModelTempGroup_final_nbinom1
)

#nbinom2 had a much lower AIC than nbinom1
#both converged but nbinom2 was retained


#Generate predictions within each group's observed range for plotting

group_temperature_ranges <- analysis_data_quad %>%
  group_by(Group) %>%
  summarise(
    min_temperature = min(
      Temperature..C.degrees.,
      na.rm = TRUE
    ),
    max_temperature = max(
      Temperature..C.degrees.,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


Temp_predictions_final_df <- group_temperature_ranges %>%
  group_by(Group) %>%
  reframe(
    Temperature = seq(
      min_temperature,
      max_temperature,
      length.out = 100
    )
  ) %>%
  ungroup() %>%
  mutate(
    Temperature_c =
      Temperature - mean_temperature_quad,
    
    Group = factor(
      Group,
      levels = levels(analysis_data_quad$Group)
    ),
    
    Site = factor(
      levels(analysis_data_quad$Site)[1],
      levels = levels(analysis_data_quad$Site)
    ),
    
    Time.Period = factor(
      levels(analysis_data_quad$Time.Period)[1],
      levels = levels(analysis_data_quad$Time.Period)
    ),
    
    Plant.ID = factor(
      levels(analysis_data_quad$Plant.ID)[1],
      levels = levels(analysis_data_quad$Plant.ID)
    ),
    
    #Any valid date is required in newdata. Setting re.form = NA
    #below excludes the date-specific random effect.
    
    Date = factor(
      levels(analysis_data_quad$Date)[1],
      levels = levels(analysis_data_quad$Date)
    )
  )


prediction_results <- predict(
  ModelTempGroup_final,
  newdata = Temp_predictions_final_df,
  type = "link",
  se.fit = TRUE,
  re.form = NA,
  allow.new.levels = TRUE
)


Temp_predictions_final_df <- Temp_predictions_final_df %>%
  mutate(
    predicted_visits = exp(
      prediction_results$fit
    ),
    
    lower_CL = exp(
      prediction_results$fit -
        1.96 * prediction_results$se.fit
    ),
    
    upper_CL = exp(
      prediction_results$fit +
        1.96 * prediction_results$se.fit
    )
  )

head(Temp_predictions_final_df)


#Plot temperature responses

pollinator_colours <- c(
  "Bumblebee" = "firebrick2",
  "Butterfly" = "lightskyblue",
  "Honeybee" = "gold1",
  "Hoverfly" = "darkorange",
  "Non-syrphid Fly" = "mediumpurple"
)


Temperature_plot <- ggplot(
  Temp_predictions_final_df,
  aes(
    x = Temperature,
    y = predicted_visits,
    colour = Group,
    fill = Group
  )
) +
  
  geom_point(
    data = analysis_data_quad,
    aes(
      x = Temperature..C.degrees.,
      y = Visits,
      colour = Group
    ),
    alpha = 0.3,
    size = 1.5,
    inherit.aes = FALSE,
    position = position_jitter(
      width = 0.2,
      height = 0
    )
  ) +
  
  geom_ribbon(
    aes(
      ymin = lower_CL,
      ymax = upper_CL
    ),
    alpha = 0.15,
    colour = NA
  ) +
  
  geom_line(
    linewidth = 1.3
  ) +
  
  facet_wrap(
    ~ Group,
    scales = "free_y"
  ) +
  
  scale_colour_manual(
    values = pollinator_colours
  ) +
  
  scale_fill_manual(
    values = pollinator_colours
  ) +
  
  labs(
    x = "Temperature (°C)",
    y = "Predicted number of visits per survey",
    title = "Temperature responses of pollinator groups"
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    
    strip.text = element_text(
      face = "bold"
    ),
    
    legend.position = "none"
  )


Temperature_plot

#######################################################################################
#6. Does visitation differ among plant species after accounting for floral abundance?
#######################################################################################

#Create species variable
#Convert Plant.ID to character first so the mapping works whether Plant.ID is currently numeric or a factor

survey_model <- survey_model |>
  mutate(
    Plant.ID = as.character(Plant.ID),
    
    Plant.Species = case_when(
      Plant.ID %in% c("1", "7")  ~ "Hypochaeris radicata",
      Plant.ID == "2"            ~ "Dianthus deltoides",
      Plant.ID %in% c("3", "8")  ~ "Thymus polytrichus",
      Plant.ID %in% c("4", "9")  ~ "Lotus corniculatus",
      Plant.ID %in% c("5", "10") ~ "Silene viscaria",
      Plant.ID == "6"            ~ "Galium saxatile",
      Plant.ID == "11"           ~ "Teucrium scorodonia",
      Plant.ID %in% c("12", "13") ~ "Galium verum",
      Plant.ID == "14"           ~ "Cirsium arvense",
      TRUE                       ~ NA_character_
    ),
    
    Plant.Common.Name = case_when(
      Plant.Species == "Hypochaeris radicata"  ~ "Yellow flatweed",
      Plant.Species == "Dianthus deltoides"    ~ "Maiden pink",
      Plant.Species == "Thymus polytrichus"    ~ "Wild thyme",
      Plant.Species == "Lotus corniculatus"    ~ "Bird's-foot trefoil",
      Plant.Species == "Silene viscaria"       ~ "Sticky catchfly",
      Plant.Species == "Galium saxatile"       ~ "Heath bedstraw",
      Plant.Species == "Teucrium scorodonia"   ~ "Woodland germander",
      Plant.Species == "Galium verum"          ~ "Lady's bedstraw",
      Plant.Species == "Cirsium arvense"     ~ "Creeping thistle",
      TRUE                                     ~ NA_character_
    )
  )


#Prepare analysis data

analysis_data_plant <- survey_model |>
  filter(
    Number.of.Flowers > 0,
    !is.na(Plant.Species)
  ) |>
  mutate(
    Plant.Species = factor(Plant.Species),
    Site = factor(Site),
    Time.Period = factor(Time.Period)
  )


# Check that all plant species are represented
table(
  analysis_data_plant$Plant.Species,
  useNA = "ifany"
)

# Check total visits for each species
analysis_data_plant |>
  group_by(Plant.Species) |>
  summarise(
    Surveys = n(),
    Total.Visits = sum(Visits),
    .groups = "drop"
  )


#Fit species model

# The offset accounts for the number of flowers available
ModelSpecies_nb <- glm.nb(
  Visits ~
    Plant.Species +
    Site +
    Time.Period +
    offset(log(Number.of.Flowers)),
  data = analysis_data_plant
)


#Fit null model without plant species

ModelAbundance_nb <- glm.nb(
  Visits ~
    Site +
    Time.Period +
    offset(log(Number.of.Flowers)),
  data = analysis_data_plant
)


#Compare both models

# Does including plant species significantly improve model fit?
anova(
  ModelAbundance_nb,
  ModelSpecies_nb,
  test = "Chisq"
)

# Compare model support
AIC(
  ModelAbundance_nb,
  ModelSpecies_nb
)


#######################################################################
#7. Which plant species are most attractive to pollinators?
#######################################################################

#Estimated visitation rates

#offset = 0 corresponds to log(1 flower)
#Therefore, estimates are predicted visits per flower per survey
Species_emmeans <- emmeans(
  ModelSpecies_nb,
  ~ Plant.Species,
  type = "response",
  offset = 0
)

Species_emmeans

#Create results table

Species_table <- as.data.frame(Species_emmeans) |>
  left_join(
    survey_model |>
      distinct(
        Plant.Species,
        Plant.Common.Name
      ),
    by = "Plant.Species"
  ) |>
  select(
    Plant.Common.Name,
    Plant.Species,
    response,
    SE,
    asymp.LCL,
    asymp.UCL
  ) |>
  rename(
    `Common name` = Plant.Common.Name,
    `Scientific name` = Plant.Species,
    `Predicted visits per flower per 20-minute survey` = response,
    `Standard error` = SE,
    `Lower 95% CI` = asymp.LCL,
    `Upper 95% CI` = asymp.UCL
  ) |>
  arrange(
    desc(`Predicted visits per flower per 20-minute survey`)
  ) |>
  mutate(
    across(
      where(is.numeric),
      ~ round(.x, 3)
    )
  )

Species_table

Species_plot <- as.data.frame(Species_emmeans) |>
  left_join(
    survey_model |>
      distinct(
        Plant.Species,
        Plant.Common.Name
      ),
    by = "Plant.Species"
  ) |>
  arrange(response) |>
  mutate(
    Plant.Label = paste0(
      Plant.Common.Name,
      "\n(",
      Plant.Species,
      ")"
    ),
    Plant.Label = factor(
      Plant.Label,
      levels = Plant.Label
    )
  )

Species_plot <- as.data.frame(Species_emmeans) |>
  dplyr::left_join(
    survey_model |>
      dplyr::distinct(
        Plant.Species,
        Plant.Common.Name
      ),
    by = "Plant.Species"
  ) |>
  dplyr::mutate(
    Type = dplyr::case_when(
      Plant.Species %in% c(
        "Dianthus deltoides",
        "Silene viscaria"
      ) ~ "Restored",
      
      TRUE ~ "Native"
    )
  )

# Create plotting data first
Species_plot <- Species_plot |>
  mutate(
    Plant.Label = factor(
      Plant.Species,
      levels = Plant.Species[order(response)]
    )
  )

Species_plot <- Species_plot |>
  mutate(
    Type = case_when(
      Plant.Species %in% c(
        "Dianthus deltoides",
        "Silene viscaria"
      ) ~ "Restored",
      TRUE ~ "Native"
    ),
    Plant.Label = factor(
      Plant.Species,
      levels = Plant.Species[order(response)]
    )
  )

# Plot
ggplot(
  Species_plot,
  aes(
    x = Plant.Label,
    y = response,
    colour = Type
  )
) +
  geom_errorbar(
    aes(
      ymin = asymp.LCL,
      ymax = asymp.UCL
    ),
    width = 0.15,
    linewidth = 0.8
  ) +
  geom_point(size = 3.5) +
  coord_flip() +
  scale_y_log10() +
  scale_colour_manual(
    values = c(
      "Native" = "forestgreen",
      "Restored" = "#FF3399"
    )
  ) +
  labs(
    x = NULL,
    y = "Adjusted predicted visitation rate",
    colour = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.y = element_text(face = "italic"),
    legend.position = "top"
  )

#######################################################################################################################
#8. Which environmental variables (temperature, solar radiation, wind speed, humidity)
#best predict pollinator visitation?
#######################################################################################################################

#Fix dates

survey_environment$Date[
  survey_environment$Date == as.Date("2026-06-06")
] <- as.Date("2026-06-05")

survey_environment$Date[
  survey_environment$Date == as.Date("2026-07-03")
] <- as.Date("2026-07-02")

#Prep data

survey_environment2 <- survey_environment %>%
  filter(
    !(Plant.ID %in% c("6", "13")),
    !is.na(Visits),
    !is.na(Number.of.Flowers),
    Number.of.Flowers > 0,
    !is.na(Temperature..C.degrees.),
    !is.na(Solar.Power..W.m.2.),
    !is.na(Humidity....RH.),
    !is.na(Wind.Speed..m.s.),
    !is.na(Date),
    !is.na(Plant.ID),
    !is.na(Time.Period)
  ) %>%
  mutate(
    Plant.ID = droplevels(factor(Plant.ID)),
    Date = factor(Date),
    Time.Period = factor(Time.Period),
    
    Temp_c = Temperature..C.degrees. -
      mean(
        Temperature..C.degrees.,
        na.rm = TRUE
      ),
    
    Solar_c = Solar.Power..W.m.2. -
      mean(
        Solar.Power..W.m.2.,
        na.rm = TRUE
      ),
    
    Humidity_c = Humidity....RH. -
      mean(
        Humidity....RH.,
        na.rm = TRUE
      ),
    
    Wind_c = Wind.Speed..m.s. -
      mean(
        Wind.Speed..m.s.,
        na.rm = TRUE
      )
  )


#Check final data

nrow(survey_environment2)

table(survey_environment2$Date)
table(survey_environment2$Plant.ID)

summary(
  survey_environment2$Number.of.Flowers
)

#Linear abiotic model

#Plant.ID controls for differences in attractiveness among plants.
#Date accounts for observations collected on the same survey day.
#Number of flowers is included as an offset so the model estimates
#visitation relative to the number of flowers available.

Environment_linear <- glmmTMB(
  Visits ~
    Temp_c +
    Solar_c +
    Humidity_c +
    Wind_c +
    Plant.ID +
    offset(
      log(Number.of.Flowers)
    ) +
    (1 | Date),
  family = nbinom2,
  data = survey_environment2
)

summary(Environment_linear)

car::Anova(
  Environment_linear,
  type = 2
)

#Test quadratic temperature response

Environment_temp_quad <- glmmTMB(
  Visits ~
    Temp_c +
    I(Temp_c^2) +
    Solar_c +
    Humidity_c +
    Wind_c +
    Plant.ID +
    offset(
      log(Number.of.Flowers)
    ) +
    (1 | Date),
  family = nbinom2,
  data = survey_environment2
)

anova(
  Environment_linear,
  Environment_temp_quad
)

AIC(
  Environment_linear,
  Environment_temp_quad
)

#If the likelihood-ratio test is significant and/or AIC is lower,
#this supports retaining a quadratic temperature response.

#Test quadratic solar-radiation response

Environment_solar_quad <- glmmTMB(
  Visits ~
    Temp_c +
    Solar_c +
    I(Solar_c^2) +
    Humidity_c +
    Wind_c +
    Plant.ID +
    offset(
      log(Number.of.Flowers)
    ) +
    (1 | Date),
  family = nbinom2,
  data = survey_environment2
)

anova(
  Environment_linear,
  Environment_solar_quad
)

AIC(
  Environment_linear,
  Environment_solar_quad
)

#If the likelihood-ratio test is not significant and AIC does not
#meaningfully improve, the quadratic solar term is not retained.

#Test quadratic humidity response

Environment_humidity_quad <- glmmTMB(
  Visits ~
    Temp_c +
    Solar_c +
    Humidity_c +
    I(Humidity_c^2) +
    Wind_c +
    Plant.ID +
    offset(
      log(Number.of.Flowers)
    ) +
    (1 | Date),
  family = nbinom2,
  data = survey_environment2
)

anova(
  Environment_linear,
  Environment_humidity_quad
)

AIC(
  Environment_linear,
  Environment_humidity_quad
)

#If the likelihood-ratio test is not significant and AIC does not
#meaningfully improve, the quadratic humidity term is not retained.

#Final environmental model

#Temperature is retained as a quadratic response.
#Solar radiation, humidity and wind speed are retained as linear effects.

Environment_final2 <- glmmTMB(
  Visits ~
    Temp_c +
    I(Temp_c^2) +
    Solar_c +
    Humidity_c +
    Wind_c +
    Plant.ID +
    offset(
      log(Number.of.Flowers)
    ) +
    (1 | Date),
  family = nbinom2,
  data = survey_environment2
)

summary(Environment_final2)

Environment_final_Anova <- car::Anova(
  Environment_final2,
  type = 2
)

Environment_final_Anova


#Does time period explain anything beyond climate?

#This model is identical to the final environmental model
#except that Time.Period is added.

Environment_final_time <- glmmTMB(
  Visits ~
    Temp_c +
    I(Temp_c^2) +
    Solar_c +
    Humidity_c +
    Wind_c +
    Plant.ID +
    Time.Period +
    offset(
      log(Number.of.Flowers)
    ) +
    (1 | Date),
  family = nbinom2,
  data = survey_environment2
)

anova(
  Environment_final2,
  Environment_final_time
)

AIC(
  Environment_final2,
  Environment_final_time
)

#If adding Time.Period does not significantly improve model fit,
#then time of day does not explain additional variation once
#the measured environmental variables are accounted for.

#Check whether Date explains meaningful variation

VarCorr(Environment_final2)

#Fit the identical model without the Date random effect.

Environment_final_no_date <- glmmTMB(
  Visits ~
    Temp_c +
    I(Temp_c^2) +
    Solar_c +
    Humidity_c +
    Wind_c +
    Plant.ID +
    offset(
      log(Number.of.Flowers)
    ),
  family = nbinom2,
  data = survey_environment2
)

AIC(
  Environment_final_no_date,
  Environment_final2
)

#With only six dates, use the variance estimate, model
#stability and AIC as the main evidence for retaining Date.

#If Environment_final2 has the lower AIC, this supports
#retaining Date as a random intercept.

#Model diagnostics

Environment_diagnostics <- simulateResiduals(
  fittedModel = Environment_final2,
  n = 1000
)

plot(Environment_diagnostics)

Environment_dispersion <- testDispersion(
  Environment_diagnostics
)

Environment_zero_inflation <- testZeroInflation(
  Environment_diagnostics
)

Environment_uniformity <- testUniformity(
  Environment_diagnostics
)

Environment_dispersion
Environment_zero_inflation
Environment_uniformity

#Check convergence

Environment_final2$sdr$pdHess

#TRUE indicates a positive-definite Hessian and supports
#successful model convergence.


#Compare nbinom1 and nbinom2

Environment_final_nbinom1 <- update(
  Environment_final2,
  family = nbinom1
)

Environment_final_nbinom1$sdr$pdHess

AIC(
  Environment_final2,
  Environment_final_nbinom1
)

#Retain the negative-binomial family with the lower AIC,
#provided that it converges successfully and diagnostics are satisfactory.

#Final outputs

summary(Environment_final2)

Environment_final_Anova

AIC(Environment_final2)

VarCorr(Environment_final2)

Environment_final2$sdr$pdHess

Environment_dispersion
Environment_zero_inflation
Environment_uniformity

############################################################################
#9. Does pollinator community composition differ between Site 1 and Site 2?
############################################################################


#Prepare pollinator community data

community_data <- pollinators %>%
  filter(
    Group != "None",
    Visits > 0
  ) %>%
  mutate(
    Group = trimws(Group),
    Species = trimws(Species),
    
    #Use species-level identity for bumblebees and butterflies
    #and broader pollinator groups for all other taxa
    
    Pollinator.Taxon = case_when(
      Group %in% c("Bumblebee", "Butterfly") &
        !is.na(Species) &
        Species != "" ~ Species,
      
      TRUE ~ Group
    ),
    
    #Create a unique sampling-occasion ID
    
    Sample.ID = interaction(
      Date,
      Site,
      Block,
      Replication,
      Time.Period,
      Plant.ID,
      drop = TRUE
    )
  )


#Check pollinator taxa used in the community analysis

sort(
  unique(
    community_data$Pollinator.Taxon
  )
)


#Build the community matrix

community_matrix_data <- community_data %>%
  group_by(
    Sample.ID,
    Site,
    Pollinator.Taxon
  ) %>%
  summarise(
    Visits = sum(
      Visits,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Pollinator.Taxon,
    values_from = Visits,
    values_fill = 0
  )


#Create metadata

community_metadata <- community_matrix_data %>%
  select(
    Sample.ID,
    Site
  ) %>%
  mutate(
    Site = factor(
      Site,
      levels = c(1, 2),
      labels = c("Site 1", "Site 2")
    )
  )


#Create numeric community matrix

community_matrix <- community_matrix_data %>%
  select(
    -Sample.ID,
    -Site
  ) %>%
  as.data.frame()

rownames(community_matrix) <-
  community_matrix_data$Sample.ID


#Check community matrix

dim(community_matrix)

table(community_metadata$Site)

head(community_matrix)

rowSums(community_matrix)

sum(
  rowSums(community_matrix) == 0
)

#No rows should have a total abundance of zero

if (
  any(
    rowSums(community_matrix) == 0
  )
) {
  stop(
    "At least one sampling occasion contains no pollinator visits."
  )
}


#Run NMDS

set.seed(123)

site_nmds <- metaMDS(
  community_matrix,
  distance = "bray",
  k = 2,
  trymax = 200,
  autotransform = FALSE,
  trace = FALSE
)

site_nmds


#Check stress

site_nmds$stress


#Extract NMDS scores

nmds_scores <- as.data.frame(
  scores(
    site_nmds,
    display = "sites"
  )
)

nmds_scores$Site <-
  community_metadata$Site


#Calculate the centroid of each site

site_centroids <- nmds_scores %>%
  group_by(Site) %>%
  summarise(
    NMDS1 = mean(NMDS1),
    NMDS2 = mean(NMDS2),
    .groups = "drop"
  )


#Create the stress label

stress_label <- paste0(
  "Stress = ",
  round(
    site_nmds$stress,
    3
  )
)


#Run PERMANOVA

site_permanova <- adonis2(
  community_matrix ~ Site,
  data = community_metadata,
  method = "bray",
  permutations = 999
)

site_permanova

#PERMANOVA tests whether pollinator community composition
#differs significantly between Site 1 and Site 2.


#Check homogeneity of multivariate dispersion

bray_dist <- vegdist(
  community_matrix,
  method = "bray"
)

site_dispersion <- betadisper(
  bray_dist,
  community_metadata$Site
)

anova(
  site_dispersion
)

permutest(
  site_dispersion,
  permutations = 999
)

#A non-significant dispersion test supports interpreting
#the PERMANOVA as a difference in community composition
#rather than a difference in within-site dispersion.


#Plot pollinator community NMDS

nmds_plot <- ggplot(
  nmds_scores,
  aes(
    x = NMDS1,
    y = NMDS2,
    colour = Site,
    fill = Site
  )
) +
  
  stat_ellipse(
    geom = "polygon",
    level = 0.95,
    alpha = 0.15,
    linewidth = 1
  ) +
  
  geom_point(
    size = 3,
    alpha = 0.8
  ) +
  
  geom_point(
    data = site_centroids,
    aes(
      x = NMDS1,
      y = NMDS2,
      colour = Site
    ),
    shape = 4,
    size = 5,
    stroke = 1.5,
    inherit.aes = FALSE
  ) +
  
  scale_colour_manual(
    values = c(
      "Site 1" = "#7B3294",
      "Site 2" = "#2E8B57"
    )
  ) +
  
  scale_fill_manual(
    values = c(
      "Site 1" = "#7B3294",
      "Site 2" = "#2E8B57"
    )
  ) +
  
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    label = stress_label,
    hjust = 1.1,
    vjust = -0.8,
    size = 4
  ) +
  
  coord_equal() +
  
  theme_classic(
    base_size = 14
  ) +
  
  theme(
    legend.position = "right",
    legend.title = element_text(
      face = "bold"
    ),
    plot.title = element_text(
      face = "bold",
      size = 16,
      hjust = 0.5
    ),
    axis.title = element_text(
      face = "bold"
    ),
    panel.border = element_rect(
      fill = NA,
      colour = "black",
      linewidth = 0.7
    )
  ) +
  
  labs(
    title = "Pollinator community composition by site",
    x = "NMDS1",
    y = "NMDS2",
    colour = "Site",
    fill = "Site"
  )


nmds_plot


############################################################################
#Plant-species pollinator community NMDS
#Do plant species cluster according to their pollinator communities,
#and do these patterns differ between Site 1 and Site 2?
############################################################################


#Prepare plant community data

plant_community_data <- pollinators %>%
  filter(
    Group != "None",
    Visits > 0
  ) %>%
  mutate(
    Group = trimws(Group),
    Species = trimws(Species),
    
    Site = factor(
      Site,
      levels = c(1, 2),
      labels = c("Site 1", "Site 2")
    ),
    
    Plant.ID = factor(
      Plant.ID
    ),
    
    #Use species-level identity for bumblebees and butterflies,
    #and broader pollinator groups for all other taxa
    
    Pollinator.Taxon = case_when(
      Group %in% c("Bumblebee", "Butterfly") &
        !is.na(Species) &
        Species != "" ~ Species,
      
      TRUE ~ Group
    )
  )


#Check taxa and plant representation

sort(
  unique(
    plant_community_data$Pollinator.Taxon
  )
)

table(
  plant_community_data$Site,
  plant_community_data$Plant.ID
)


#Create plant x pollinator community matrix

#Each row represents one Plant.ID at one site.
#Each column represents one pollinator taxon.
#Cell values are total numbers of visits.

plant_matrix_data <- plant_community_data %>%
  group_by(
    Site,
    Plant.ID,
    Pollinator.Taxon
  ) %>%
  summarise(
    Visits = sum(
      Visits,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Pollinator.Taxon,
    values_from = Visits,
    values_fill = 0
  )


#Create metadata

plant_metadata <- plant_matrix_data %>%
  select(
    Site,
    Plant.ID
  )


#Create numeric community matrix

plant_matrix <- plant_matrix_data %>%
  select(
    -Site,
    -Plant.ID
  ) %>%
  as.data.frame()

rownames(plant_matrix) <- paste(
  plant_matrix_data$Site,
  plant_matrix_data$Plant.ID,
  sep = "_Plant_"
)


#Check community matrix

dim(plant_matrix)

head(plant_matrix)

rowSums(plant_matrix)

sum(
  rowSums(plant_matrix) == 0
)

#No rows should have a total abundance of zero

if (
  any(
    rowSums(plant_matrix) == 0
  )
) {
  stop(
    "At least one plant-site row contains no pollinator visits."
  )
}


#Run NMDS

set.seed(123)

plant_nmds <- metaMDS(
  plant_matrix,
  distance = "bray",
  k = 2,
  trymax = 200,
  autotransform = FALSE,
  trace = FALSE
)

plant_nmds


#Check stress

plant_nmds$stress

plant_stress_label <- paste0(
  "Stress = ",
  round(
    plant_nmds$stress,
    3
  )
)


#Extract NMDS scores

plant_nmds_scores <- as.data.frame(
  scores(
    plant_nmds,
    display = "sites"
  )
)

plant_nmds_scores$Site <-
  plant_metadata$Site

plant_nmds_scores$Plant.ID <-
  plant_metadata$Plant.ID


#Add plant species names

plant_names_nmds <- c(
  "1" = "Hypochaeris radicata",
  "2" = "Dianthus deltoides",
  "3" = "Thymus drucei",
  "4" = "Lotus corniculatus",
  "5" = "Silene viscaria",
  "6" = "Galium saxatile",
  "7" = "Hypochaeris radicata",
  "8" = "Thymus drucei",
  "9" = "Lotus corniculatus",
  "10" = "Silene viscaria",
  "11" = "Teucrium scorodonia",
  "12" = "Galium verum",
  "13" = "Galium verum",
  "14" = "Cirsium arvense"
)

plant_nmds_scores <- plant_nmds_scores %>%
  mutate(
    Plant.Species = unname(
      plant_names_nmds[
        as.character(Plant.ID)
      ]
    )
  )


#Calculate site centroids

plant_site_centroids <- plant_nmds_scores %>%
  group_by(Site) %>%
  summarise(
    NMDS1 = mean(NMDS1),
    NMDS2 = mean(NMDS2),
    .groups = "drop"
  )


#Plot plant pollinator communities

plant_nmds_plot <- ggplot(
  plant_nmds_scores,
  aes(
    x = NMDS1,
    y = NMDS2,
    colour = Site,
    shape = Plant.ID
  )
) +
  
  geom_point(
    size = 4,
    alpha = 0.9
  ) +
  
  geom_point(
    data = plant_site_centroids,
    aes(
      x = NMDS1,
      y = NMDS2,
      colour = Site
    ),
    shape = 4,
    size = 6,
    stroke = 1.5,
    inherit.aes = FALSE
  ) +
  
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    label = plant_stress_label,
    hjust = 1.1,
    vjust = -0.8,
    size = 4
  ) +
  
  scale_colour_manual(
    values = c(
      "Site 1" = "#7B3294",
      "Site 2" = "#2E8B57"
    )
  ) +
  
  coord_equal() +
  
  theme_classic(
    base_size = 14
  ) +
  
  theme(
    legend.position = "right",
    legend.title = element_text(
      face = "bold"
    ),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    axis.title = element_text(
      face = "bold"
    ),
    panel.border = element_rect(
      fill = NA,
      colour = "black",
      linewidth = 0.7
    )
  ) +
  
  labs(
    title = "Plant species pollinator communities by site",
    x = "NMDS1",
    y = "NMDS2",
    colour = "Site",
    shape = "Plant ID"
  )


plant_nmds_plot


#Label points with plant species

plant_nmds_labelled_plot <- ggplot(
  plant_nmds_scores,
  aes(
    x = NMDS1,
    y = NMDS2,
    colour = Site
  )
) +
  
  geom_point(
    size = 4,
    alpha = 0.9
  ) +
  
  geom_text(
    aes(
      label = Plant.Species
    ),
    vjust = -0.8,
    size = 3.5,
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  
  geom_point(
    data = plant_site_centroids,
    aes(
      x = NMDS1,
      y = NMDS2,
      colour = Site
    ),
    shape = 4,
    size = 6,
    stroke = 1.5,
    inherit.aes = FALSE
  ) +
  
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    label = plant_stress_label,
    hjust = 1.1,
    vjust = -0.8,
    size = 4
  ) +
  
  scale_colour_manual(
    values = c(
      "Site 1" = "#7B3294",
      "Site 2" = "#2E8B57"
    )
  ) +
  
  coord_equal() +
  
  theme_classic(
    base_size = 14
  ) +
  
  theme(
    legend.position = "right",
    legend.title = element_text(
      face = "bold"
    ),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    axis.title = element_text(
      face = "bold"
    ),
    panel.border = element_rect(
      fill = NA,
      colour = "black",
      linewidth = 0.7
    )
  ) +
  
  labs(
    title = "Plant species pollinator communities by site",
    x = "NMDS1",
    y = "NMDS2",
    colour = "Site"
  )


plant_nmds_labelled_plot


#Test whether plant-associated pollinator communities differ between sites

plant_bray_dist <- vegdist(
  plant_matrix,
  method = "bray"
)

plant_site_permanova <- adonis2(
  plant_matrix ~ Site,
  data = plant_metadata,
  method = "bray",
  permutations = 999
)

plant_site_permanova


#Check homogeneity of dispersion between sites

plant_site_dispersion <- betadisper(
  plant_bray_dist,
  plant_metadata$Site
)

anova(
  plant_site_dispersion
)

permutest(
  plant_site_dispersion,
  permutations = 999
)
##############################################################################
#10. Which co-flowering plant species share pollinators with restored species?
# Analysed separately by site
##############################################################################
#Pollinator taxon:
#Use Species where a valid species name is available
#Otherwise use the broader Group

#Check required columns

required_columns <- c(
  "Species",
  "Group",
  "Visits",
  "Site",
  "Plant.ID"
)

missing_columns <- setdiff(
  required_columns,
  names(pollinators)
)

if (length(missing_columns) > 0) {
  stop(
    paste(
      "The pollinators dataframe is missing:",
      paste(missing_columns, collapse = ", ")
    )
  )
}

#Plant names

plant_names <- c(
  "1"  = "Hypochaeris radicata",
  "2"  = "Dianthus deltoides",
  "3"  = "Thymus drucei",
  "4"  = "Lotus corniculatus",
  "5"  = "Silene viscaria",
  "6"  = "Galium saxatile",
  "7"  = "Hypochaeris radicata",
  "8"  = "Thymus drucei",
  "9"  = "Lotus corniculatus",
  "10" = "Silene viscaria",
  "11" = "Teucrium scorodonia",
  "12" = "Galium verum",
  "13" = "Galium verum",
  "14" = "Cirsium arvense"
)

restored_species <- c(
  "Dianthus deltoides",
  "Silene viscaria"
)


#Prep data

sharing_data <- pollinators %>%
  mutate(
    Group = trimws(as.character(Group)),
    Species = trimws(as.character(Species)),
    Plant.ID = as.character(Plant.ID),
    
    Site = case_when(
      as.character(Site) %in% c("1", "Site 1") ~ "Site 1",
      as.character(Site) %in% c("2", "Site 2") ~ "Site 2",
      TRUE ~ NA_character_
    ),
    
    Plant.Species = unname(
      plant_names[Plant.ID]
    ),
    
    Pollinator.Taxon = case_when(
      !is.na(Species) &
        Species != "" &
        Species != "NA" &
        Species != "Unknown" ~ Species,
      
      TRUE ~ Group
    )
  ) %>%
  filter(
    !is.na(Site),
    !is.na(Plant.Species),
    !is.na(Pollinator.Taxon),
    Pollinator.Taxon != "",
    Group != "None",
    !is.na(Visits),
    Visits > 0
  )


# Check the final pollinator categories

print(
  sort(
    unique(
      sharing_data$Pollinator.Taxon
    )
  )
)


#Create plant x pollinator presence matrix

make_presence_matrix <- function(data, selected_site) {
  
  site_table <- data %>%
    filter(
      Site == selected_site
    ) %>%
    group_by(
      Plant.Species,
      Pollinator.Taxon
    ) %>%
    summarise(
      Total.Visits = sum(
        Visits,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    mutate(
      Present = as.integer(
        Total.Visits > 0
      )
    ) %>%
    select(
      Plant.Species,
      Pollinator.Taxon,
      Present
    ) %>%
    pivot_wider(
      names_from = Pollinator.Taxon,
      values_from = Present,
      values_fill = 0
    )
  
  if (nrow(site_table) == 0) {
    stop(
      paste(
        "No usable observations were found for",
        selected_site
      )
    )
  }
  
  site_matrix <- site_table %>%
    column_to_rownames(
      "Plant.Species"
    ) %>%
    as.matrix()
  
  storage.mode(site_matrix) <- "numeric"
  
  return(site_matrix)
}


site1_presence <- make_presence_matrix(
  sharing_data,
  "Site 1"
)

site2_presence <- make_presence_matrix(
  sharing_data,
  "Site 2"
)


#Match pollinator columns between sites

all_pollinator_taxa <- union(
  colnames(site1_presence),
  colnames(site2_presence)
)


add_missing_taxa <- function(
    matrix_object,
    required_taxa
) {
  
  missing_taxa <- setdiff(
    required_taxa,
    colnames(matrix_object)
  )
  
  if (length(missing_taxa) > 0) {
    
    added_columns <- matrix(
      0,
      nrow = nrow(matrix_object),
      ncol = length(missing_taxa),
      dimnames = list(
        rownames(matrix_object),
        missing_taxa
      )
    )
    
    matrix_object <- cbind(
      matrix_object,
      added_columns
    )
  }
  
  matrix_object[
    ,
    required_taxa,
    drop = FALSE
  ]
}


site1_presence <- add_missing_taxa(
  site1_presence,
  all_pollinator_taxa
)

site2_presence <- add_missing_taxa(
  site2_presence,
  all_pollinator_taxa
)


# Confirm columns match

stopifnot(
  identical(
    colnames(site1_presence),
    colnames(site2_presence)
  )
)

###################################################
#11. Calculate shared pollinator taxa across sites
###################################################

#Rows = Site 2 plants
#Columns = Site 1 plants
#Matrix multiplication counts how many pollinator taxa are present on both plants.
cross_site_shared_matrix <-
  site2_presence %*% t(site1_presence)


#Convert to long format

cross_site_shared <- as.data.frame(
  cross_site_shared_matrix
) %>%
  rownames_to_column(
    "Site2.Plant"
  ) %>%
  pivot_longer(
    cols = -Site2.Plant,
    names_to = "Site1.Plant",
    values_to = "Shared"
  )


#Calc pollinator richness for each plant

site1_richness <- rowSums(
  site1_presence > 0
)

site2_richness <- rowSums(
  site2_presence > 0
)


#Order plants of restored plants first followed by native plants alphabetically

site1_order <- c(
  intersect(
    restored_species,
    rownames(site1_presence)
  ),
  sort(
    setdiff(
      rownames(site1_presence),
      restored_species
    )
  )
)

site2_order <- c(
  intersect(
    restored_species,
    rownames(site2_presence)
  ),
  sort(
    setdiff(
      rownames(site2_presence),
      restored_species
    )
  )
)


cross_site_shared <- cross_site_shared %>%
  mutate(
    Site1.Plant = factor(
      Site1.Plant,
      levels = site1_order
    ),
    
    Site2.Plant = factor(
      Site2.Plant,
      levels = rev(site2_order)
    )
  )

#Axis labels
#Number in brackets = pollinator richness for that plant

site1_labels <- function(x) {
  
  x <- as.character(x)
  
  richness <- unname(
    site1_richness[x]
  )
  
  label_text <- paste0(
    x,
    " (",
    richness,
    ")"
  )
  
  ifelse(
    x %in% restored_species,
    
    paste0(
      "<span style='color:#7B3294'><b><i>",
      label_text,
      "</i></b></span>"
    ),
    
    paste0(
      "<i>",
      label_text,
      "</i>"
    )
  )
}


site2_labels <- function(x) {
  
  x <- as.character(x)
  
  richness <- unname(
    site2_richness[x]
  )
  
  label_text <- paste0(
    x,
    " (",
    richness,
    ")"
  )
  
  ifelse(
    x %in% restored_species,
    
    paste0(
      "<span style='color:#7B3294'><b><i>",
      label_text,
      "</i></b></span>"
    ),
    
    paste0(
      "<i>",
      label_text,
      "</i>"
    )
  )
}


#Create heatmap

cross_site_shared_heatmap <- ggplot(
  cross_site_shared,
  aes(
    x = Site1.Plant,
    y = Site2.Plant,
    fill = Shared
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.8
  ) +
  geom_text(
    aes(
      label = Shared
    ),
    size = 4
  ) +
  scale_fill_gradient(
    low = "lightgoldenrodyellow",
    high = "darkorange1",
    name = "Shared\npollinator taxa"
  ) +
  scale_x_discrete(
    labels = site1_labels,
    drop = FALSE
  ) +
  scale_y_discrete(
    labels = site2_labels,
    drop = FALSE
  ) +
  labs(
    title = "Pollinator overlap between plant species across sites",
    subtitle = paste0(
      "X-axis = Site 1 plants | ",
      "Y-axis = Site 2 plants"
    ),
    x = "Plant species at Site 1",
    y = "Plant species at Site 2",
    caption = paste0(
      "Cell values show the number of shared pollinator taxa. ",
      "Numbers after plant names show total pollinator richness."
    )
  ) +
  theme_classic(
    base_size = 13
  ) +
  theme(
    axis.text.x = ggtext::element_markdown(
      angle = 45,
      hjust = 1,
      size = 9
    ),
    axis.text.y = ggtext::element_markdown(
      size = 9
    ),
    axis.title = element_text(
      face = "bold"
    ),
    axis.ticks = element_blank(),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 16
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 11
    ),
    plot.caption = element_text(
      hjust = 0,
      size = 9
    ),
    legend.title = element_text(
      face = "bold"
    ),
    panel.border = element_rect(
      fill = NA,
      colour = "black",
      linewidth = 0.7
    )
  )

print(
  cross_site_shared_heatmap
)


##############################################################################
#12. Which co-flowering species have the greatest potential to compete with
# or facilitate pollination of restored species?
#Pollinator sharing is defined as:the same pollinator taxon visiting a restored and native plant species
# at the same site and during the same time period.
##############################################################################

#Restored plant species
  mutate(
    Plant.ID = as.character(Plant.ID),
    
    Plant.Species = case_when(
      Plant.ID %in% c("1", "7")   ~ "Hypochaeris radicata",
      Plant.ID == "2"             ~ "Dianthus deltoides",
      Plant.ID %in% c("3", "8")   ~ "Thymus drucei",
      Plant.ID %in% c("4", "9")   ~ "Lotus corniculatus",
      Plant.ID %in% c("5", "10")  ~ "Silene viscaria",
      Plant.ID == "6"             ~ "Galium saxatile",
      Plant.ID == "11"            ~ "Teucrium scorodonia",
      Plant.ID %in% c("12", "13") ~ "Galium verum",
      Plant.ID == "14"            ~ "Cirsium arvense",
      TRUE ~ NA_character_
    )
  ) %>%
  select(
    Site,
    Time.Period,
    Plant.Species
  ) %>%
  distinct()

if (!exists("plant_occasion")) {
  stop("The object 'plant_occasion' does not exist.")
}


restored_species <- c(
  "Dianthus deltoides",
  "Silene viscaria"
)


#Check required objects and columns
if (!exists("plant_occasion")) {
  stop("The object 'plant_occasion' does not exist.")
}

if (!exists("interaction_data")) {
  stop("The object 'interaction_data' does not exist.")
}


required_occasion_columns <- c(
  "Site",
  "Time.Period",
  "Plant.Species"
)

required_interaction_columns <- c(
  "Site",
  "Time.Period",
  "Plant.Species",
  "Pollinator.Taxon",
  "Visits"
)

``
missing_occasion_columns <- setdiff(
  required_occasion_columns,
  names(plant_occasion)
)

missing_interaction_columns <- setdiff(
  required_interaction_columns,
  names(interaction_data)
)


if (length(missing_occasion_columns) > 0) {
  stop(
    paste(
      "plant_occasion is missing:",
      paste(missing_occasion_columns, collapse = ", ")
    )
  )
}

if (length(missing_interaction_columns) > 0) {
  stop(
    paste(
      "interaction_data is missing:",
      paste(missing_interaction_columns, collapse = ", ")
    )
  )
}


#Clean and standardise data

plant_occasion_clean <- plant_occasion %>%
  mutate(
    Site = case_when(
      as.character(Site) %in% c("1", "Site 1") ~ "Site 1",
      as.character(Site) %in% c("2", "Site 2") ~ "Site 2",
      TRUE ~ as.character(Site)
    ),
    
    Time.Period = trimws(
      as.character(Time.Period)
    ),
    
    Plant.Species = trimws(
      as.character(Plant.Species)
    )
  ) %>%
  filter(
    !is.na(Site),
    !is.na(Time.Period),
    Time.Period != "",
    !is.na(Plant.Species),
    Plant.Species != ""
  ) %>%
  distinct(
    Site,
    Time.Period,
    Plant.Species
  )


interaction_data_clean <- interaction_data %>%
  mutate(
    Site = case_when(
      as.character(Site) %in% c("1", "Site 1") ~ "Site 1",
      as.character(Site) %in% c("2", "Site 2") ~ "Site 2",
      TRUE ~ as.character(Site)
    ),
    
    Time.Period = trimws(
      as.character(Time.Period)
    ),
    
    Plant.Species = trimws(
      as.character(Plant.Species)
    ),
    
    Pollinator.Taxon = trimws(
      as.character(Pollinator.Taxon)
    )
  ) %>%
  filter(
    !is.na(Site),
    !is.na(Time.Period),
    Time.Period != "",
    !is.na(Plant.Species),
    Plant.Species != "",
    !is.na(Pollinator.Taxon),
    Pollinator.Taxon != "",
    Pollinator.Taxon != "None",
    !is.na(Visits),
    Visits > 0
  )


#Identify all restored and native plants present within each site x time period

restored_periods <- plant_occasion_clean %>%
  filter(
    Plant.Species %in% restored_species
  ) %>%
  transmute(
    Site,
    Time.Period,
    Restored.Species = Plant.Species
  ) %>%
  distinct()


native_periods <- plant_occasion_clean %>%
  filter(
    !Plant.Species %in% restored_species
  ) %>%
  transmute(
    Site,
    Time.Period,
    Native.Species = Plant.Species
  ) %>%
  distinct()


#Calc co-flowering/temporal overlap
#A restored and native species overlap when both were present
#at the same site during the same time period

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


#Summarise restored-plant interactions w/i each site x time period

restored_period_interactions <- interaction_data_clean %>%
  filter(
    Plant.Species %in% restored_species
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
    Restored.Species = Plant.Species
  )


#Summarise native-plant interactions w/i each site x time period

native_period_interactions <- interaction_data_clean %>%
  filter(
    !Plant.Species %in% restored_species
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
    Native.Species = Plant.Species
  )


#Identify shared pollinator use
#same site, same time period, same pollinator taxon

period_sharing <- restored_period_interactions %>%
  inner_join(
    native_period_interactions,
    by = c(
      "Site",
      "Time.Period",
      "Pollinator.Taxon"
    )
  )

#Summarise sharing for each restored-native pair

shared_pollinator_summary <- period_sharing %>%
  group_by(
    Site,
    Restored.Species,
    Native.Species
  ) %>%
  summarise(
    # Number of distinct pollinator taxa shared in at least one time period
    Shared.Pollinator.Taxa = n_distinct(
      Pollinator.Taxon
    ),
    
    # Number of unique time-period × pollinator combinations
    Time.Period.Sharing.Events = n_distinct(
      interaction(
        Time.Period,
        Pollinator.Taxon,
        drop = TRUE
      )
    ),
    
    # Number of time periods during which any sharing occurred
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


#Complete potential_interactions table

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
      ~ replace_na(.x, 0)
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
    desc(Shared.Pollinator.Taxa),
    desc(Time.Period.Sharing.Events),
    Native.Species
  )


# View the complete results table

print(
  potential_interactions,
  n = Inf
)


#Screen all native species

all_native_species_screened <- potential_interactions %>%
  distinct(
    Native.Species
  ) %>%
  arrange(
    Native.Species
  )


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


cat("\nAll native species screened:\n")

print(
  all_native_species_screened,
  n = Inf
)


cat("\nNative species sharing at least one pollinator taxon:\n")

print(
  native_species_with_sharing,
  n = Inf
)


# Full screening table

native_screening_results <- potential_interactions %>%
  select(
    Site,
    Restored.Species,
    Native.Species,
    Overlapping.Time.Periods,
    Shared.Pollinator.Taxa,
    Time.Period.Sharing.Events,
    Time.Periods.With.Shared.Pollinators,
    Proportion.Time.Periods.With.Sharing,
    Shared.Any.Pollinators
  ) %>%
  arrange(
    Native.Species,
    Site,
    Restored.Species
  )


print(
  native_screening_results,
  n = Inf
)


#Plot Bubble Plot

native_order <- potential_interactions %>%
  filter(
    Shared.Any.Pollinators
  ) %>%
  distinct(
    Native.Species
  ) %>%
  arrange(
    Native.Species
  ) %>%
  pull(
    Native.Species
  ) %>%
  as.character()


if (length(native_order) == 0) {
  stop(
    paste(
      "No native species shared pollinator taxa",
      "with the restored species."
    )
  )
}


interaction_plot_data <- potential_interactions %>%
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
    
    Native.Species = factor(
      Native.Species,
      levels = rev(
        native_order
      )
    )
  )


#What rows will appear as bubbles

print(
  interaction_plot_data %>%
    select(
      Site,
      Restored.Species,
      Native.Species,
      Shared.Pollinator.Taxa,
      Time.Period.Sharing.Events,
      Time.Periods.With.Shared.Pollinators
    ),
  n = Inf
)


#X-axis range

x_max <- max(
  interaction_plot_data$Shared.Pollinator.Taxa,
  na.rm = TRUE
)

x_upper <- ceiling(
  x_max
) + 0.5


#Grey background for Maiden Pink at Site 2 as it isn't present

background_x <- seq(
  0,
  x_upper,
  length.out = 150
)

tile_width <- background_x[2] -
  background_x[1]


absence_background <- expand_grid(
  x = background_x,
  
  Native.Species = factor(
    rev(native_order),
    levels = rev(native_order)
  ),
  
  Site = factor(
    "Site 2",
    levels = c(
      "Site 1",
      "Site 2"
    )
  ),
  
  Restored.Species = factor(
    "Dianthus deltoides",
    levels = c(
      "Dianthus deltoides",
      "Silene viscaria"
    )
  )
)


middle_native_species <- native_order[
  ceiling(
    length(native_order) / 2
  )
]


absence_label <- data.frame(
  x = x_upper / 2,
  
  Native.Species = factor(
    middle_native_species,
    levels = rev(native_order)
  ),
  
  Site = factor(
    "Site 2",
    levels = c(
      "Site 1",
      "Site 2"
    )
  ),
  
  Restored.Species = factor(
    "Dianthus deltoides",
    levels = c(
      "Dianthus deltoides",
      "Silene viscaria"
    )
  ),
  
  label = "Dianthus deltoides
absent from Site 2"
)


#Create bubble plot

interaction_plot <- ggplot(
  interaction_plot_data,
  aes(
    x = Shared.Pollinator.Taxa,
    y = Native.Species,
    size = Time.Period.Sharing.Events
  )
) +
  
  geom_tile(
    data = absence_background,
    aes(
      x = x,
      y = Native.Species
    ),
    inherit.aes = FALSE,
    width = tile_width * 1.1,
    height = 1.05,
    fill = "grey95",
    colour = NA
  ) +
  
  geom_text(
    data = absence_label,
    aes(
      x = x,
      y = Native.Species,
      label = label
    ),
    inherit.aes = FALSE,
    colour = "grey30",
    fontface = "italic",
    size = 4
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
    space = "free_y",
    switch = "y",
    drop = FALSE
  ) +
  
  scale_y_discrete(
    position = "right",
    drop = TRUE
  ) +
  
  scale_x_continuous(
    breaks = seq(
      0,
      ceiling(x_max),
      by = 1
    ),
    
    limits = c(
      0,
      x_upper
    ),
    
    expand = expansion(
      mult = c(0, 0)
    )
  ) +
  
  scale_size_continuous(
    name = "Shared time periods",
    range = c(
      3,
      10
    )
  ) +
  
  labs(
    x = "Number of shared pollinator taxa",
    y = "Native plant species",
    
    title = paste0(
      "Shared pollinator use between restored ",
      "and native plant species"
    ),
    
    subtitle = paste0(
      "Larger points indicate more time-period × ",
      "pollinator sharing events"
    )
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    axis.text.y.right = element_text(
      face = "italic"
    ),
    
    axis.title.y.right = element_text(
      margin = margin(
        l = 10
      )
    ),
    
    strip.text.y.left = element_text(
      face = "italic",
      angle = 90
    ),
    
    strip.text.x = element_text(
      face = "bold"
    ),
    
    strip.placement = "outside",
    
    strip.background = element_rect(
      fill = "white",
      colour = "black"
    ),
    
    strip.switch.pad.grid = grid::unit(
      0.15,
      "cm"
    ),
    
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    
    legend.position = "right",
    
    plot.margin = margin(
      t = 10,
      r = 20,
      b = 10,
      l = 30
    )
  )

print(
  interaction_plot
)

#############################
#13. Maiden Pink Pollinators
#############################

# Pollinator taxa recorded visiting Maiden pink at Site 1
dianthus_taxa <- community_data %>%
  filter(
    Site == 1,
    Plant.ID == 2
  ) %>%
  pull(Pollinator.Taxon) %>%
  unique()

sort(dianthus_taxa)


#All pollinator taxa recorded anywhere at Site 2
site2_taxa <- community_data %>%
  filter(Site == 2) %>%
  pull(Pollinator.Taxon) %>%
  unique()

sort(site2_taxa)


#Maiden pink pollinators that were also present at Site 2
shared_taxa <- intersect(dianthus_taxa, site2_taxa)

shared_taxa
length(shared_taxa)


#Maiden pink pollinators that were absent from Site 2
missing_taxa <- setdiff(dianthus_taxa, site2_taxa)

missing_taxa
length(missing_taxa)


#Number of Maiden pink visits made by each missing taxon
community_data %>%
  filter(
    Site == 1,
    Plant.ID == 2,
    Pollinator.Taxon %in% missing_taxa
  ) %>%
  group_by(Pollinator.Taxon) %>%
  summarise(
    Total.Visits = sum(Visits, na.rm = TRUE),
    .groups = "drop"
  )


#Total Maiden pink visits and percentage made by missing taxa
community_data %>%
  filter(
    Site == 1,
    Plant.ID == 2
  ) %>%
  summarise(
    Total.Visits = sum(Visits, na.rm = TRUE),
    Missing.Taxa.Visits = sum(
      Visits[Pollinator.Taxon %in% missing_taxa],
      na.rm = TRUE
    ),
    Percentage.Missing = 
      100 * Missing.Taxa.Visits / Total.Visits
  )

community_data %>%
  filter(Pollinator.Taxon == "Thymelicus sylvestris") %>%
  group_by(Site, Plant.ID) %>%
  summarise(
    Visits = sum(Visits),
    .groups = "drop"
  ) %>%
  arrange(Site, desc(Visits))

community_data %>%
  filter(Site == 1) %>%
  group_by(Plant.ID) %>%
  summarise(
    Total.Visits = sum(Visits),
    Missing.Visits = sum(Visits[Pollinator.Taxon %in% missing_taxa]),
    Percent = 100 * Missing.Visits / Total.Visits,
    .groups = "drop"
  ) %>%
  arrange(desc(Percent))

community_data %>%
  filter(Plant.ID %in% c(1, 7)) %>%
  group_by(Site, Plant.ID) %>%
  summarise(
    Total.Visits = sum(Visits),
    Mean.Flowers = mean(Number.of.Flowers),
    Visits.Per.Flower = Total.Visits / Mean.Flowers,
    .groups = "drop"
  )

site1_connectance <- bipartite::networklevel(
  site1_species_web,
  index = "connectance"
)

site2_connectance <- bipartite::networklevel(
  site2_species_web,
  index = "connectance"
)

site1_connectance
site2_connectance

########################
#Extra Analysis/Appendix
########################

#############################################################################
#14. What percentage of each plant's visits came from each pollinator taxon?
#############################################################################

#Calculate percentages
plant_pollinator <- community_data %>%
  group_by(Plant.ID, Pollinator.Taxon) %>%
  summarise(
    Visits = sum(Visits),
    .groups = "drop"
  ) %>%
  group_by(Plant.ID) %>%
  mutate(
    Percent = 100 * Visits / sum(Visits)
  ) %>%
  arrange(Plant.ID, desc(Percent))
plant_pollinator
library(ggplot2)

ggplot(
  plant_pollinator,
  aes(
    x = Plant.ID,
    y = Percent,
    fill = Pollinator.Taxon
  )
) +
  geom_col() +
  labs(
    x = "Plant",
    y = "Percentage of visits"
  ) +
  theme_classic()

plant_pollinator %>%
  group_by(Plant.ID) %>%
  slice_max(
    Percent,
    n = 1
  )

######################################################
#15. How many pollinator taxa are needed to reach 80%?
######################################################

plant_pollinator <- community_data %>%
  group_by(Plant.ID, Pollinator.Taxon) %>%
  summarise(
    Visits = sum(Visits),
    .groups = "drop"
  )
dominance <- plant_pollinator %>%
  group_by(Plant.ID) %>%
  arrange(desc(Visits), .by_group = TRUE) %>%
  mutate(
    Percent = 100 * Visits / sum(Visits),
    CumPercent = cumsum(Percent),
    Rank = row_number()
  )
dominance80 <- dominance %>%
  filter(CumPercent >= 80) %>%
  slice(1) %>%
  ungroup() %>%
  select(
    Plant.ID,
    Pollinators.for.80 = Rank,
    CumPercent
  )

dominance80

################################
#15. Sticky Catchfly
#Pollinators and Visits per Site
################################

#Sticky catchfly is Plant ID 5 at Site 1
#and Plant ID 10 at Site 2.

catchfly <- community_data %>%
  filter(
    (Site == 1 & Plant.ID == 5) |
      (Site == 2 & Plant.ID == 10)
  )


#Summarise visitation and pollinator richness

catchfly_summary <- catchfly %>%
  group_by(
    Site,
    Plant.ID
  ) %>%
  summarise(
    Total.Visits = sum(
      Visits,
      na.rm = TRUE
    ),
    
    Mean.Flowers = mean(
      Number.of.Flowers,
      na.rm = TRUE
    ),
    
    Pollinator.Richness = n_distinct(
      Pollinator.Taxon
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    Visits.Per.Flower =
      Total.Visits / Mean.Flowers
  )

catchfly_summary


#Pollinator taxa visiting Sticky Catchfly at each site

catchfly_pollinators <- catchfly %>%
  group_by(
    Site,
    Pollinator.Taxon
  ) %>%
  summarise(
    Visits = sum(
      Visits,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  arrange(
    Site,
    desc(Visits)
  )

catchfly_pollinators


####################################
#16. Overall descriptive statistics
####################################

#Total visits

total_visits <- sum(
  pollinators$Visits,
  na.rm = TRUE
)


#Total number of pollinator records

n_records <- nrow(
  pollinators
)


#Number of plant IDs receiving pollinator visits

n_plant_ids <- pollinators %>%
  filter(
    Group != "None"
  ) %>%
  distinct(
    Plant.ID
  ) %>%
  nrow()


#Number of pollinator taxa at the taxonomic resolution
#used in the network analysis

n_pollinator_taxa <- species_network_data %>%
  distinct(
    Pollinator.Label
  ) %>%
  nrow()


#Number of unique survey occasions

n_surveys <- pollinators %>%
  distinct(
    Date,
    Site,
    Plant.ID,
    Time.Period,
    Replication,
    Start.Time
  ) %>%
  nrow()


total_visits
n_records
n_plant_ids
n_pollinator_taxa
n_surveys


#######################################
#17. Restored-plant visitor assemblages
#######################################

#Summarise pollinator taxa visiting the restored plant populations

focal_pollinators <- pollinators %>%
  filter(
    Plant.ID %in% c(
      2,
      5,
      10
    ),
    Group != "None"
  ) %>%
  mutate(
    
    Focal_Plant = case_when(
      Plant.ID == 2 ~
        "Dianthus deltoides, Site 1",
      
      Plant.ID == 5 ~
        "Silene viscaria, Site 1",
      
      Plant.ID == 10 ~
        "Silene viscaria, Site 2"
    ),
    
    Pollinator_Taxon = case_when(
      !is.na(Species) &
        Species != "" &
        Species != "None" ~
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


############################################
#18. Force-directed plant-pollinator network
############################################

#This section requires tidygraph and ggraph.
#Keep only if the force-directed network is used as a final figure.


#Create interaction edges

network_edges <- species_network_data %>%
  group_by(
    Plant.Label,
    Pollinator.Label
  ) %>%
  summarise(
    Visits = sum(
      Visits,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  filter(
    Visits > 0
  ) %>%
  transmute(
    from = paste0(
      "Plant:",
      Plant.Label
    ),
    to = paste0(
      "Pollinator:",
      Pollinator.Label
    ),
    Visits = Visits
  )


#Create plant nodes

plant_nodes <- species_network_data %>%
  distinct(
    Plant.Label
  ) %>%
  transmute(
    name = paste0(
      "Plant:",
      Plant.Label
    ),
    
    Label = Plant.Label,
    
    Node.Type = "Plant",
    
    Plant.Status = if_else(
      Plant.Label %in% c(
        "D. deltoides",
        "S. viscaria"
      ),
      "Restored plant",
      "Co-flowering plant"
    ),
    
    Pollinator.Group = NA_character_
  )


#Create pollinator nodes

pollinator_nodes <- species_network_data %>%
  group_by(
    Pollinator.Label
  ) %>%
  summarise(
    Pollinator.Group = first(
      na.omit(
        as.character(Group)
      )
    ),
    .groups = "drop"
  ) %>%
  transmute(
    name = paste0(
      "Pollinator:",
      Pollinator.Label
    ),
    
    Label = Pollinator.Label,
    
    Node.Type = "Pollinator",
    
    Plant.Status = NA_character_,
    
    Pollinator.Group =
      Pollinator.Group
  )


#Combine nodes

network_nodes <- bind_rows(
  plant_nodes,
  pollinator_nodes
)


#Check for duplicate node names

network_nodes %>%
  count(
    name
  ) %>%
  filter(
    n > 1
  )


#Create graph

pollinator_graph <- tbl_graph(
  nodes = network_nodes,
  edges = network_edges,
  directed = FALSE
) %>%
  activate(
    nodes
  ) %>%
  mutate(
    Interaction.Strength =
      centrality_degree(
        weights = Visits
      )
  )


#Plot network

set.seed(123)

pollinator_network_plot <- ggraph(
  pollinator_graph,
  layout = "fr"
) +
  
  geom_edge_link(
    aes(
      width = Visits,
      alpha = Visits
    ),
    colour = "grey65",
    show.legend = FALSE
  ) +
  
  scale_edge_width(
    range = c(
      0.3,
      3
    )
  ) +
  
  scale_edge_alpha(
    range = c(
      0.2,
      0.75
    )
  ) +
  
  geom_node_point(
    aes(
      filter =
        Node.Type == "Pollinator",
      
      size =
        Interaction.Strength,
      
      colour =
        Pollinator.Group
    ),
    shape = 16
  ) +
  
  geom_node_point(
    aes(
      filter =
        Node.Type == "Plant" &
        Plant.Status ==
        "Co-flowering plant",
      
      size =
        Interaction.Strength
    ),
    shape = 22,
    fill = "darkseagreen3",
    colour = "black",
    stroke = 1
  ) +
  
  geom_node_point(
    aes(
      filter =
        Node.Type == "Plant" &
        Plant.Status ==
        "Restored plant",
      
      size =
        Interaction.Strength
    ),
    shape = 22,
    fill = "deeppink2",
    colour = "black",
    stroke = 1.2
  ) +
  
  geom_node_text(
    aes(
      label = Label,
      filter =
        Node.Type == "Plant"
    ),
    repel = TRUE,
    fontface = "italic",
    size = 4
  ) +
  
  geom_node_text(
    aes(
      label = Label,
      filter =
        Node.Type == "Pollinator"
    ),
    repel = TRUE,
    size = 3
  ) +
  
  scale_size_continuous(
    range = c(
      3,
      10
    ),
    name = "Total visits"
  ) +
  
  labs(
    colour = "Pollinator group"
  ) +
  
  theme_void() +
  
  theme(
    legend.position = "right",
    plot.margin = margin(
      20,
      20,
      20,
      20
    )
  )

pollinator_network_plot


################
#Appendix tables
################


#Survey effort table

survey_table <- pollinators %>%
  group_by(
    Site,
    Plant.ID
  ) %>%
  summarise(
    `Number of Replicates` =
      n_distinct(
        Replication
      ),
    
    `Number of Surveys` =
      n_distinct(
        interaction(
          Date,
          Start.Time,
          Time.Period,
          drop = TRUE
        )
      ),
    
    .groups = "drop"
  ) %>%
  arrange(
    Site,
    Plant.ID
  )

survey_table


#Pollinator taxonomic list

pollinator_list <- pollinators %>%
  filter(
    Group != "None"
  ) %>%
  distinct(
    Order,
    Family,
    Species,
    Group
  ) %>%
  arrange(
    Order,
    Family,
    Species
  )

pollinator_list


#Check whether any pollinator species were assigned
#to more than one pollinator group

group_inconsistencies <- pollinators %>%
  filter(
    !is.na(Species),
    Group != "None"
  ) %>%
  group_by(
    Species
  ) %>%
  summarise(
    n_groups =
      n_distinct(
        Group
      ),
    
    groups =
      paste(
        sort(
          unique(Group)
        ),
        collapse = ", "
      ),
    
    .groups = "drop"
  ) %>%
  filter(
    n_groups > 1
  )

group_inconsistencies


###########################################
#Plant x pollinator-group interaction table
###########################################

all_groups <- c(
  "Bumblebee",
  "Butterfly",
  "Honeybee",
  "Hoverfly",
  "Mining Bee",
  "Non-syrphid Fly",
  "Sweat Bee",
  "Wasp"
)

interaction_table <- pollinators %>%
  filter(
    Group != "None"
  ) %>%
  mutate(
    Group = factor(
      Group,
      levels = all_groups
    )
  ) %>%
  group_by(
    Plant.ID,
    Group,
    .drop = FALSE
  ) %>%
  summarise(
    Visits = sum(
      Visits,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Group,
    values_from = Visits,
    values_fill = 0,
    names_expand = TRUE
  ) %>%
  arrange(
    Plant.ID
  )

interaction_table


#############################################
#Plant x pollinator-species interaction table
#############################################

species_interaction_table <- pollinators %>%
  filter(
    Group != "None",
    !is.na(Species),
    Species != ""
  ) %>%
  group_by(
    Plant.ID,
    Species
  ) %>%
  summarise(
    Visits = sum(
      Visits,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Species,
    values_from = Visits,
    values_fill = 0
  ) %>%
  arrange(
    Plant.ID
  )

print(
  species_interaction_table,
  width = Inf
)


#Long-format version containing interactions only

species_long <- species_interaction_table %>%
  pivot_longer(
    cols = -Plant.ID,
    names_to = "Pollinator species",
    values_to = "Visits"
  ) %>%
  filter(
    Visits > 0
  ) %>%
  arrange(
    Plant.ID,
    desc(Visits)
  )

species_long


################################
#Abiotic descriptive statistics
################################

abiotic_summary <- pollinators %>%
  summarise(
    
    Temperature_Mean =
      mean(
        Temperature..C.degrees.,
        na.rm = TRUE
      ),
    
    Temperature_SD =
      sd(
        Temperature..C.degrees.,
        na.rm = TRUE
      ),
    
    Temperature_Min =
      min(
        Temperature..C.degrees.,
        na.rm = TRUE
      ),
    
    Temperature_Max =
      max(
        Temperature..C.degrees.,
        na.rm = TRUE
      ),
    
    Solar_Mean =
      mean(
        Solar.Power..W.m.2.,
        na.rm = TRUE
      ),
    
    Solar_SD =
      sd(
        Solar.Power..W.m.2.,
        na.rm = TRUE
      ),
    
    Solar_Min =
      min(
        Solar.Power..W.m.2.,
        na.rm = TRUE
      ),
    
    Solar_Max =
      max(
        Solar.Power..W.m.2.,
        na.rm = TRUE
      ),
    
    Humidity_Mean =
      mean(
        Humidity....RH.,
        na.rm = TRUE
      ),
    
    Humidity_SD =
      sd(
        Humidity....RH.,
        na.rm = TRUE
      ),
    
    Humidity_Min =
      min(
        Humidity....RH.,
        na.rm = TRUE
      ),
    
    Humidity_Max =
      max(
        Humidity....RH.,
        na.rm = TRUE
      ),
    
    Wind_Mean =
      mean(
        Wind.Speed..m.s.,
        na.rm = TRUE
      ),
    
    Wind_SD =
      sd(
        Wind.Speed..m.s.,
        na.rm = TRUE
      ),
    
    Wind_Min =
      min(
        Wind.Speed..m.s.,
        na.rm = TRUE
      ),
    
    Wind_Max =
      max(
        Wind.Speed..m.s.,
        na.rm = TRUE
      )
  )

abiotic_summary


####################################
#Correlation among abiotic variables
####################################

abiotic <- pollinators %>%
  select(
    Temperature =
      Temperature..C.degrees.,
    
    Solar =
      Solar.Power..W.m.2.,
    
    Humidity =
      Humidity....RH.,
    
    Wind =
      Wind.Speed..m.s.
  )


#Pearson correlation coefficients

cor_matrix <- cor(
  abiotic,
  use = "complete.obs",
  method = "pearson"
)

round(
  cor_matrix,
  2
)

##################################
#Appendix: Overall pollinator groups
##################################

appendix_colours <- c(
  "#E78FB3",
  "#86BBD8",
  "#F2C14E",
  "#8FBC8F",
  "#A78BC1",
  "#E9A66A",
  "#72B7A8"
)


#Visitation by time period

Time_appendix <- pollinators %>%
  group_by(
    Time.Period
  ) %>%
  summarise(
    Total.Visits = sum(
      Visits,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  ggplot(
    aes(
      x = Time.Period,
      y = Total.Visits,
      fill = Time.Period
    )
  ) +
  geom_col(
    width = 0.7
  ) +
  scale_fill_manual(
    values = appendix_colours[1:3]
  ) +
  labs(
    x = "Time period",
    y = "Total pollinator visits"
  ) +
  theme_classic(
    base_size = 13
  ) +
  theme(
    legend.position = "none"
  )

Time_appendix


#Visitation by site

Site_appendix <- pollinators %>%
  group_by(
    Site
  ) %>%
  summarise(
    Total.Visits = sum(
      Visits,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  ggplot(
    aes(
      x = Site,
      y = Total.Visits,
      fill = Site
    )
  ) +
  geom_col(
    width = 0.65
  ) +
  scale_fill_manual(
    values = appendix_colours[4:5]
  ) +
  labs(
    x = "Site",
    y = "Total pollinator visits"
  ) +
  theme_classic(
    base_size = 13
  ) +
  theme(
    legend.position = "none"
  )

Site_appendix


#Visitation by pollinator group

Pollinator_appendix <- pollinators %>%
  filter(
    Group != "None"
  ) %>%
  group_by(
    Group
  ) %>%
  summarise(
    Total.Visits = sum(
      Visits,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  ggplot(
    aes(
      x = reorder(
        Group,
        Total.Visits
      ),
      y = Total.Visits,
      fill = Group
    )
  ) +
  geom_col(
    width = 0.7
  ) +
  coord_flip() +
  scale_fill_manual(
    values = rep(
      appendix_colours,
      length.out = 20
    )
  ) +
  labs(
    x = "Pollinator group",
    y = "Total pollinator visits"
  ) +
  theme_classic(
    base_size = 13
  ) +
  theme(
    legend.position = "none"
  )

Pollinator_appendix
