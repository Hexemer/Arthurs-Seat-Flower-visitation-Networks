##################################
#00 Data Setup
##################################


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
      paste(
        missing_columns,
        collapse = ", "
      )
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
    
    # Standardise group terminology from other fly to non-syrphid fly
    # (true fly not hoverfly)
    
    Group = case_when(
      Group == "Other Fly" ~ "Non-syrphid Fly",
      TRUE ~ Group
    )
  )


#Prep variables

pollinators$Group <- factor(
  trimws(
    pollinators$Group
  )
)

pollinators$Plant.ID <- factor(
  pollinators$Plant.ID
)

pollinators$Time.Period <- factor(
  pollinators$Time.Period
)

pollinators$Site <- factor(
  pollinators$Site
)

pollinators$Block <- factor(
  pollinators$Block
)

pollinators$Replication <- factor(
  pollinators$Replication
)


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


#Create pollinator list

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


#Create survey-level dataset

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
      
      TRUE ~
        Date
    ),
    
    Time.Period = factor(
      Time.Period
    ),
    
    Site = factor(
      Site
    ),
    
    Plant.ID = factor(
      Plant.ID
    ),
    
    Block = factor(
      Block
    ),
    
    Replication = factor(
      Replication
    )
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
    Date = droplevels(
      factor(
        Date
      )
    ),
    
    Time.Period = droplevels(
      factor(
        Time.Period
      )
    ),
    
    Site = droplevels(
      factor(
        Site
      )
    ),
    
    Plant.ID = droplevels(
      factor(
        Plant.ID
      )
    )
  )


#Check survey-level data

head(
  survey_data
)

# Number of individual surveys

nrow(
  survey_data
)

# Distribution of visits

summary(
  survey_data$Visits
)

# Number of zero-visit surveys

sum(
  survey_data$Visits == 0
)

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

table(
  survey_data$Date
)

# Check plant sample sizes

table(
  survey_data$Plant.ID
)

# Check plant IDs within sites

table(
  survey_data$Plant.ID,
  survey_data$Site
)


#Plant names

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

plant_names

