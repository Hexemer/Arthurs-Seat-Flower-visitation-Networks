#######################################################################################################################
# AIM 1
# Which insects visit maiden pink and sticky catchfly, and to what extent are their flower visitors
# shared with non-restored plant species?
#######################################################################################################################


###########
# Packages
###########

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggtext)
library(MASS)
library(emmeans)


############
# Load data
############

pollinators <- read.csv(
  "Raw Data/Pollinators.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)

#Although the dataframe and some internal object names use the term
#"pollinator", the observations represent flower visitation rather than
#confirmed pollination.


#####################
# Plant species names
#####################

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

restored_species <- c(
  "Dianthus deltoides",
  "Silene viscaria"
)


################################
# 1. Prepare flower-visitor data
################################

#Use species where a reliable species name is available.
#Otherwise use genus, family or order.

community_data <- pollinators %>%
  filter(
    Group != "None",
    Visits > 0
  ) %>%
  mutate(
    Plant.ID = as.character(Plant.ID),
    
    Site = case_when(
      as.character(Site) %in% c("1", "Site 1") ~ "Site 1",
      as.character(Site) %in% c("2", "Site 2") ~ "Site 2",
      TRUE ~ NA_character_
    ),
    
    Plant.Species = unname(
      plant_key[Plant.ID]
    ),
    
    Flower.Visitor.Taxon = case_when(
      !is.na(Species) &
        Species != "" &
        Species != "None" &
        Species != "Unknown" ~ trimws(Species),
      
      !is.na(Genus) &
        Genus != "" ~ paste0(
          trimws(Genus),
          " sp."
        ),
      
      !is.na(Family) &
        Family != "" ~ paste0(
          trimws(Family),
          " family"
        ),
      
      !is.na(Order) &
        Order != "" ~ paste0(
          trimws(Order),
          " order"
        ),
      
      TRUE ~ as.character(Group)
    )
  ) %>%
  filter(
    !is.na(Site),
    !is.na(Plant.Species),
    !is.na(Flower.Visitor.Taxon)
  )


#############################################
# 2. Raw visitation to restored plant species
#############################################

#Maiden pink = Plant ID 2 at Site 1
#Sticky catchfly = Plant ID 5 at Site 1 and Plant ID 10 at Site 2

restored_visitation <- community_data %>%
  filter(
    Plant.ID %in% c(
      "2",
      "5",
      "10"
    )
  ) %>%
  mutate(
    Restored.Population = case_when(
      Plant.ID == "2" ~
        "Dianthus deltoides, Site 1",
      
      Plant.ID == "5" ~
        "Silene viscaria, Site 1",
      
      Plant.ID == "10" ~
        "Silene viscaria, Site 2"
    )
  )


#Total raw visits to each restored population

restored_raw_totals <- restored_visitation %>%
  group_by(
    Restored.Population
  ) %>%
  summarise(
    Total.Visits = sum(
      Visits,
      na.rm = TRUE
    ),
    
    Flower.Visitor.Richness = n_distinct(
      Flower.Visitor.Taxon
    ),
    
    .groups = "drop"
  )

restored_raw_totals


###########################################################
# 3. Which insects visited maiden pink and sticky catchfly?
###########################################################

#Summarise flower-visitor taxa visiting the restored plant populations

restored_visitor_assemblages <- restored_visitation %>%
  group_by(
    Restored.Population,
    Flower.Visitor.Taxon,
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
    Restored.Population,
    desc(Visits)
  )

restored_visitor_assemblages


#Dominant flower visitor at each restored population

dominant_restored_visitors <- restored_visitor_assemblages %>%
  group_by(
    Restored.Population
  ) %>%
  slice_max(
    Visits,
    n = 1,
    with_ties = TRUE
  ) %>%
  ungroup()

dominant_restored_visitors


##################################
# Maiden pink flower visitors
##################################

#Prepare data

dianthus_plot_data <- maiden_pink_visitors %>%
  arrange(
    Visits
  ) %>%
  mutate(
    Flower.Visitor.Taxon = factor(
      Flower.Visitor.Taxon,
      levels = Flower.Visitor.Taxon
    )
  )


#Plot

MaidenPink_plot <- ggplot(
  dianthus_plot_data,
  aes(
    x = Visits,
    y = Flower.Visitor.Taxon,
    fill = Flower.Visitor.Taxon
  )
) +
  
  geom_col(
    width = 0.7,
    show.legend = FALSE
  ) +
  
  geom_text(
    aes(
      label = Visits
    ),
    hjust = -0.3,
    size = 4
  ) +
  
  scale_fill_manual(
    values = visitor_colours
  ) +
  
  scale_x_continuous(
    expand = expansion(
      mult = c(
        0,
        0.12
      )
    )
  ) +
  
  labs(
    x = "Number of visits",
    y = "Floral visitor taxon"
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    axis.text.y = element_text(
      size = 11
    ),
    
    axis.title = element_text(
      size = 13
    )
  )

MaidenPink_plot

visitor_colours <- c(
  "Bombus hortorum"        = "#F8766D",
  "Bombus hypnorum"        = "#ED8319",
  "Bombus lapidarius"      = "#D89B00",
  "Bombus lucorum"         = "#B6A800",
  "Bombus pascuorum"       = "#8DB000",
  "Bombus pratorum"        = "#45B500",
  "Bombus terrestris"      = "#00BA5A",
  "Honeybee"               = "#00BE8A",
  "Hoverfly"               = "#00B9AC",
  "Lycaena phlaeas"        = "#00AFCB",
  "Mining Bee"             = "#00A5E3",
  "Other Fly"              = "#3298ED",
  "Polyommatus icarus"     = "#8B80F9",
  "Sweat Bee"              = "#C66AF4",
  "Thymelicus sylvestris"  = "#E657E4",
  "Vanessa cardui"          = "#F24DB7",
  "Wasp"                    = "#F45A8D"
)

maiden_pink_visitors <- restored_visitor_assemblages %>%
  filter(
    Restored.Population ==
      "Dianthus deltoides, Site 1"
  ) %>%
  arrange(
    desc(Visits)
  )

maiden_pink_visitors


#Plot

) %>%
  group_by(
    Flower.Visitor.Taxon
  ) %>%
  summarise(
    Total.Visits = sum(
      Visits,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    Site.Status = if_else(
      Flower.Visitor.Taxon %in% site2_taxa,
      "Also recorded at Site 2",
      "Only recorded at Site 1"
    )
  ) %>%
  arrange(
    Total.Visits
  ) %>%
  mutate(
    Flower.Visitor.Taxon = factor(
      Flower.Visitor.Taxon,
      levels = Flower.Visitor.Taxon
    )
  )


#Plot

MaidenPink_plot <- ggplot(
  dianthus_plot_data,
  aes(
    x = Total.Visits,
    y = Flower.Visitor.Taxon,
    fill = Site.Status
  )
) +
  
  geom_col(
    width = 0.7
  ) +
  
  geom_text(
    aes(
      label = Total.Visits
    ),
    hjust = -0.3,
    size = 4
  ) +
  
  scale_fill_manual(
    values = c(
      "Only recorded at Site 1" = "#1B9E77",
      "Also recorded at Site 2" = "#7570B3"
    ),
    breaks = c(
      "Only recorded at Site 1",
      "Also recorded at Site 2"
    ),
    name = "Visitor distribution"
  ) +
  
  scale_x_continuous(
    expand = expansion(
      mult = c(
        0,
        0.12
      )
    )
  ) +
  
  labs(
    x = "Number of visits",
    y = "Floral visitor taxon"
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    axis.text.y = element_text(
      size = 11
    ),
    
    axis.title = element_text(
      size = 13
    ),
    
    legend.position = "bottom",
    
    legend.title = element_text(
      face = "bold"
    )
  )

####################################
# 5. Sticky catchfly flower visitors
####################################

catchfly_visitors <- restored_visitor_assemblages %>%
  filter(
    Restored.Population %in% c(
      "Silene viscaria, Site 1",
      "Silene viscaria, Site 2"
    )
  ) %>%
  mutate(
    Site = case_when(
      Restored.Population ==
        "Silene viscaria, Site 1" ~ "Site 1",
      
      Restored.Population ==
        "Silene viscaria, Site 2" ~ "Site 2"
    )
  )

catchfly_visitors


#Dominant sticky catchfly visitor at each site

catchfly_dominant_visitors <- catchfly_visitors %>%
  group_by(
    Site
  ) %>%
  slice_max(
    Visits,
    n = 1,
    with_ties = TRUE
  ) %>%
  ungroup()

catchfly_dominant_visitors


#Check whether any butterflies visited sticky catchfly

catchfly_butterflies <- catchfly_visitors %>%
  filter(
    Group == "Butterfly"
  )

catchfly_butterflies


#Plot

catchfly_plot_data <- catchfly_visitors %>%
  group_by(
    Flower.Visitor.Taxon
  ) %>%
  mutate(
    Total.Taxon.Visits = sum(
      Visits,
      na.rm = TRUE
    )
  ) %>%
  ungroup() %>%
  arrange(
    Total.Taxon.Visits
  ) %>%
  mutate(
    Flower.Visitor.Taxon = factor(
      Flower.Visitor.Taxon,
      levels = unique(
        Flower.Visitor.Taxon
      )
    )
  )


Catchfly_plot <- ggplot(
  catchfly_plot_data,
  aes(
    x = Visits,
    y = Flower.Visitor.Taxon,
    fill = Site
  )
) +
  
  geom_col(
    position = position_dodge(
      width = 0.75
    ),
    width = 0.68
  ) +
  
  geom_text(
    aes(
      label = Visits
    ),
    position = position_dodge(
      width = 0.75
    ),
    hjust = -0.25,
    size = 4
  ) +
  
  scale_x_continuous(
    expand = expansion(
      mult = c(
        0,
        0.12
      )
    )
  ) +
  
  labs(
    x = "Number of visits",
    y = "Flower-visitor taxon",
    fill = NULL
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    axis.text.y = element_text(
      size = 11
    ),
    
    axis.title = element_text(
      size = 13
    ),
    
    legend.position = "bottom"
  )

Catchfly_plot


#############################################################
# 6. Were butterflies present elsewhere in the study system?
#############################################################

#This checks whether the absence of butterfly visits to sticky catchfly
#simply reflects an absence of butterflies from the study system.

butterfly_visitation <- community_data %>%
  filter(
    Group == "Butterfly"
  ) %>%
  group_by(
    Site,
    Plant.Species,
    Flower.Visitor.Taxon
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

butterfly_visitation


#Total butterfly visits by site

butterflies_by_site <- community_data %>%
  filter(
    Group == "Butterfly"
  ) %>%
  group_by(
    Site
  ) %>%
  summarise(
    Butterfly.Visits = sum(
      Visits,
      na.rm = TRUE
    ),
    
    Butterfly.Taxa = n_distinct(
      Flower.Visitor.Taxon
    ),
    
    .groups = "drop"
  )

butterflies_by_site


##################################################################
# 7. Survey-level dataset for floral-abundance-adjusted visitation
##################################################################

#Raw visit totals do not account for differences in floral abundance
#or sampling structure, so visitation is also compared using a model
#with the number of flowers included as an offset.

survey_data <- pollinators %>%
  mutate(
    Date = as.Date(
      Date,
      format = "%d/%m/%Y"
    ),
    
    #Correct known date-entry errors
    
    Date = case_when(
      Date == as.Date("2026-06-06") ~
        as.Date("2026-06-05"),
      
      Date == as.Date("2026-07-03") ~
        as.Date("2026-07-02"),
      
      TRUE ~ Date
    ),
    
    Plant.ID = as.character(
      Plant.ID
    ),
    
    Plant.Species = unname(
      plant_key[Plant.ID]
    ),
    
    Site = factor(Site),
    
    Time.Period = factor(
      Time.Period
    )
  ) %>%
  
  group_by(
    Date,
    Block,
    Replication,
    Time.Period,
    Site,
    Plant.ID,
    Plant.Species,
    Start.Time
  ) %>%
  
  summarise(
    Visits = sum(
      Visits,
      na.rm = TRUE
    ),
    
    Number.of.Flowers = first(
      Number.of.Flowers
    ),
    
    Minutes = first(
      Minutes
    ),
    
    .groups = "drop"
  )


######################################################################################
# 8. Does visitation differ among plant species after accounting for floral abundance?
######################################################################################

#Plant IDs 6 and 13 received no flower-visitor visits.
#They are retained in descriptive summaries but excluded from
#the inferential count model.

analysis_data_plant <- survey_data %>%
  filter(
    !(Plant.ID %in% c(
      "6",
      "13"
    )),
    Number.of.Flowers > 0,
    !is.na(Plant.Species)
  ) %>%
  mutate(
    Plant.Species = factor(
      Plant.Species
    ),
    
    Site = factor(
      Site
    ),
    
    Time.Period = factor(
      Time.Period
    )
  )


#Fit negative-binomial model
#Number of flowers is included as an offset so the model estimates
#visitation relative to floral abundance.

ModelSpecies_nb <- MASS::glm.nb(
  Visits ~
    Plant.Species +
    Site +
    Time.Period +
    offset(
      log(Number.of.Flowers)
    ),
  data = analysis_data_plant
)

summary(
  ModelSpecies_nb
)


#Estimated marginal means for plant species

Species_emmeans <- emmeans::emmeans(
  ModelSpecies_nb,
  ~ Plant.Species,
  type = "response"
)

Species_emmeans


###############################################
# 9. Rank plant species by predicted visitation
###############################################

#This produces the ranking used to compare maiden pink and sticky catchfly
#after accounting for floral abundance.

plant_visitation_ranking <- as.data.frame(
  Species_emmeans
) %>%
  arrange(
    desc(response)
  ) %>%
  mutate(
    Rank = row_number()
  ) %>%
  dplyr::select(
    Rank,
    Plant.Species,
    response,
    SE,
    asymp.LCL,
    asymp.UCL
  )

plant_visitation_ranking


#Show only the restored species

restored_visitation_ranking <- plant_visitation_ranking %>%
  filter(
    Plant.Species %in%
      restored_species
  )

restored_visitation_ranking


###################################################
# 10. Plot predicted visitation among plant species
###################################################

Species_plot <- as.data.frame(
  Species_emmeans
) %>%
  mutate(
    Type = case_when(
      Plant.Species %in%
        restored_species ~
        "Restored",
      
      TRUE ~
        "Non-restored"
    ),
    
    Plant.Label = factor(
      Plant.Species,
      levels = Plant.Species[
        order(response)
      ]
    )
  )


ggplot(
  Species_plot,
  aes(
    x = response,
    y = Plant.Label
  )
) +
  
  geom_point(
    size = 3
  ) +
  
  geom_errorbarh(
    aes(
      xmin = asymp.LCL,
      xmax = asymp.UCL
    ),
    height = 0.2
  ) +
  
  labs(
    x = "Predicted visits per flower per 20-minute survey",
    y = "Plant species"
  ) +
  
  theme_classic(
    base_size = 13
  )


#################################################
# 11. Flower-visitor composition by plant species
#################################################

plant_visitor_composition <- community_data %>%
  group_by(
    Site,
    Plant.ID,
    Plant.Species,
    Flower.Visitor.Taxon
  ) %>%
  
  summarise(
    Visits = sum(
      Visits,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  
  group_by(
    Site,
    Plant.ID
  ) %>%
  
  mutate(
    Percent =
      100 * Visits /
      sum(Visits)
  ) %>%
  
  ungroup()


#Plot floral visitor composition by plant species
#Restored plant species are shown in dark magenta.

ggplot(
  plant_visitor_composition,
  aes(
    x = Plant.Species,
    y = Percent,
    fill = Flower.Visitor.Taxon
  )
) +
  
  geom_col() +
  
  facet_wrap(
    ~ Site,
    scales = "free_x"
  ) +
  
  scale_x_discrete(
    labels = function(x) {
      
      ifelse(
        x %in%
          restored_species,
        
        paste0(
          "<span style='color:#8B008B'><i>",
          x,
          "</i></span>"
        ),
        
        paste0(
          "<i>",
          x,
          "</i>"
        )
      )
    }
  ) +
  
  labs(
    x = "Plant species",
    y = "Percentage of flower-visitor visits",
    fill = "Flower-visitor taxon"
  ) +
  
  theme_classic() +
  
  theme(
    axis.text.x =
      ggtext::element_markdown(
        angle = 45,
        hjust = 1
      ),
    
    strip.text =
      element_text(
        face = "bold"
      )
  )


############################################################################
# 12. Flower visitors shared between restored and non-restored plant species
############################################################################

#For Aim 1, sharing is used to identify which flower-visitor taxa
#were recorded on both a restored and a non-restored plant species.
#The timing of sharing is examined separately under Aim 3.


#Create plant x flower-visitor presence data

visitor_presence <- community_data %>%
  group_by(
    Site,
    Plant.Species,
    Flower.Visitor.Taxon
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
  )


#Restored plant interactions

restored_interactions <- visitor_presence %>%
  filter(
    Plant.Species %in%
      restored_species
  ) %>%
  
  rename(
    Restored.Species =
      Plant.Species
  )


#Non-restored plant interactions

non_restored_interactions <- visitor_presence %>%
  filter(
    !(Plant.Species %in%
        restored_species)
  ) %>%
  
  rename(
    Non.Restored.Species =
      Plant.Species
  )


#Identify shared flower-visitor taxa

shared_visitors <- restored_interactions %>%
  inner_join(
    non_restored_interactions,
    by = c(
      "Site",
      "Flower.Visitor.Taxon"
    ),
    suffix = c(
      ".Restored",
      ".Non.Restored"
    )
  ) %>%
  
  dplyr::select(
    Site,
    Restored.Species,
    Non.Restored.Species,
    Flower.Visitor.Taxon
  ) %>%
  
  distinct()


shared_visitors


##########################################
# 13. Number of shared flower-visitor taxa
##########################################

shared_visitor_summary <- shared_visitors %>%
  group_by(
    Site,
    Restored.Species,
    Non.Restored.Species
  ) %>%
  
  summarise(
    Shared.Flower.Visitor.Taxa =
      n_distinct(
        Flower.Visitor.Taxon
      ),
    
    .groups = "drop"
  ) %>%
  
  arrange(
    Site,
    Restored.Species,
    desc(
      Shared.Flower.Visitor.Taxa
    )
  )

shared_visitor_summary


############################################
# 14. Which flower-visitor taxa were shared?
############################################

shared_visitor_taxa <- shared_visitors %>%
  arrange(
    Site,
    Restored.Species,
    Non.Restored.Species,
    Flower.Visitor.Taxon
  )

shared_visitor_taxa


#########################
# 15. Final Aim 1 outputs
#########################

#Raw visitation to restored plants

restored_raw_totals


#Flower-visitor assemblages of restored plants

restored_visitor_assemblages


#Dominant visitor to each restored population

dominant_restored_visitors


#Sticky catchfly dominant visitor at each site

catchfly_dominant_visitors


#Check butterfly visitation to sticky catchfly

catchfly_butterflies


#Butterfly visitation elsewhere in the study system

butterfly_visitation
butterflies_by_site


#Predicted visitation after accounting for floral abundance

plant_visitation_ranking
restored_visitation_ranking


#Shared flower visitors

shared_visitor_summary
shared_visitor_taxa


#Figures

MaidenPink_plot
Catchfly_plot
