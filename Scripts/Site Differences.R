# AIM 2
# To what extent do flower-visitation patterns and flower-visitor
# community composition vary between sites?


# Packages

library(dplyr)
library(tidyr)
library(ggplot2)
library(glmmTMB)
library(MASS)
library(car)
library(emmeans)
library(DHARMa)
library(vegan)


# Load data

pollinators <- read.csv(
  "Raw Data/Pollinators.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)

#Although the dataframe and some internal object names use the term
#"pollinator", the observations represent flower visitation rather than
#confirmed pollination.


# Plant species names

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


# 1. Prepare survey-level data

survey_data <- pollinators %>%
  mutate(
    Date = as.Date(
      Date,
      format = "%d/%m/%Y"
    ),
    
    #Correct confirmed date-entry errors
    
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
    #Sum visits from all flower-visitor groups recorded during
    #the same 20-minute plant survey
    
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


# Check survey-level data

nrow(survey_data)

table(
  survey_data$Site
)

table(
  survey_data$Plant.ID,
  survey_data$Site
)


# 2. Prepare site-visitation model data

#Plant IDs 6 and 13 received no flower-visitor visits in any survey.
#They are retained in descriptive summaries but excluded from
#inferential count models because they provide no within-plant
#variation in visitation.

survey_model <- survey_data %>%
  filter(
    !(Plant.ID %in% c(
      "6",
      "13"
    ))
  ) %>%
  mutate(
    Plant.ID = droplevels(Plant.ID),
    Date = droplevels(Date),
    Site = droplevels(Site),
    Time.Period = droplevels(Time.Period)
  )


# 3. Test block

#Block represents the first two weeks and second two weeks of sampling.

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

#Addition of block was not significant.


# 4. Flower visitation between sites

#Date and Plant ID are included as random effects.

ModelTimeSite <- glmmTMB(
  Visits ~
    Time.Period +
    Site +
    (1 | Date) +
    (1 | Plant.ID),
  family = nbinom2,
  data = survey_model
)


# Compare nbinom1 and nbinom2

ModelTimeSite_nbinom1 <- update(
  ModelTimeSite,
  family = nbinom1
)

AIC(
  ModelTimeSite,
  ModelTimeSite_nbinom1
)

ModelTimeSite_nbinom1$sdr$pdHess

#nbinom1 was retained because it converged cleanly,
#had satisfactory diagnostics and had the lower AIC.

ModelTimeSite_final <- ModelTimeSite_nbinom1


# Overall fixed-effect tests

TimeSite_Anova <- car::Anova(
  ModelTimeSite_final,
  type = 2
)

TimeSite_Anova


# Estimated visitation by site

Site_emmeans <- emmeans(
  ModelTimeSite_final,
  pairwise ~ Site,
  type = "response"
)

Site_emmeans

#Predicted visits per survey for each site

Site_emmeans$emmeans

#Site comparison as a response-scale ratio

Site_emmeans$contrasts


# 5. Site-model diagnostics

set.seed(123)

TimeSite_residuals <- simulateResiduals(
  fittedModel = ModelTimeSite_final,
  n = 1000
)

testUniformity(
  TimeSite_residuals
)

testDispersion(
  TimeSite_residuals
)

testZeroInflation(
  TimeSite_residuals
)


# Four-panel diagnostic figure

par(
  mfrow = c(2, 2),
  mar = c(6, 5, 3, 2),
  cex = 0.75
)

#QQ residual diagnostic

plotQQunif(
  TimeSite_residuals
)

#Residuals against fitted values

plotResiduals(
  TimeSite_residuals,
  form = fitted(ModelTimeSite_final)
)

#Dispersion diagnostic

testDispersion(
  TimeSite_residuals,
  plot = TRUE
)

#Zero-inflation diagnostic

testZeroInflation(
  TimeSite_residuals,
  plot = TRUE
)

#Reset plotting layout

par(
  mfrow = c(1, 1)
)


# 6. Descriptive visitation by site

site_visitation_summary <- survey_data %>%
  group_by(
    Site
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

site_visitation_summary


# Sampling effort

site_sampling_effort <- survey_data %>%
  group_by(
    Site
  ) %>%
  summarise(
    Surveys = n(),
    
    Floral.Units.Observed = sum(
      Number.of.Flowers,
      na.rm = TRUE
    ),
    
    Survey.Minutes = sum(
      Minutes,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

site_sampling_effort


# 7. Plant species shared between sites

plants_by_site <- pollinators %>%
  mutate(
    Plant.ID = as.character(
      Plant.ID
    ),
    
    Plant.Species = unname(
      plant_key[Plant.ID]
    ),
    
    Site = case_when(
      as.character(Site) %in% c(
        "1",
        "Site 1"
      ) ~ "Site 1",
      
      as.character(Site) %in% c(
        "2",
        "Site 2"
      ) ~ "Site 2",
      
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(Site),
    !is.na(Plant.Species)
  ) %>%
  distinct(
    Site,
    Plant.Species
  )


shared_plant_species <- plants_by_site %>%
  count(
    Plant.Species
  ) %>%
  filter(
    n == 2
  ) %>%
  pull(
    Plant.Species
  )

shared_plant_species

length(
  shared_plant_species
)


# 8. Flower-visitor community composition

#Use species-level identity for bumblebees and butterflies,
#and broader flower-visitor groups for all other taxa.

community_data <- pollinators %>%
  filter(
    Group != "None",
    Visits > 0
  ) %>%
  mutate(
    Group = trimws(
      Group
    ),
    
    Species = trimws(
      Species
    ),
    
    Flower.Visitor.Taxon = case_when(
      Group %in% c(
        "Bumblebee",
        "Butterfly"
      ) &
        !is.na(Species) &
        Species != "" ~
        Species,
      
      TRUE ~
        Group
    ),
    
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


# Check taxa used in community analysis

sort(
  unique(
    community_data$Flower.Visitor.Taxon
  )
)


# Build community matrix

community_matrix_data <- community_data %>%
  group_by(
    Sample.ID,
    Site,
    Flower.Visitor.Taxon
  ) %>%
  summarise(
    Visits = sum(
      Visits,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Flower.Visitor.Taxon,
    values_from = Visits,
    values_fill = 0
  )


# Metadata

community_metadata <- community_matrix_data %>%
  dplyr::select(
    Sample.ID,
    Site
  ) %>%
  mutate(
    Site = factor(
      Site,
      levels = c(
        1,
        2
      ),
      labels = c(
        "Site 1",
        "Site 2"
      )
    )
  )


# Numeric community matrix

community_matrix <- community_matrix_data %>%
  dplyr::select(
    -Sample.ID,
    -Site
  ) %>%
  as.data.frame()

rownames(
  community_matrix
) <- community_matrix_data$Sample.ID


# Check matrix

dim(
  community_matrix
)

table(
  community_metadata$Site
)

sum(
  rowSums(
    community_matrix
  ) == 0
)


# 9. NMDS

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


# Extract NMDS scores

nmds_scores <- as.data.frame(
  scores(
    site_nmds,
    display = "sites"
  )
)

nmds_scores$Site <-
  community_metadata$Site


# Site centroids

site_centroids <- nmds_scores %>%
  group_by(
    Site
  ) %>%
  summarise(
    NMDS1 = mean(NMDS1),
    NMDS2 = mean(NMDS2),
    .groups = "drop"
  )


# Stress label

stress_label <- paste0(
  "Stress = ",
  round(
    site_nmds$stress,
    3
  )
)


# 10. PERMANOVA

set.seed(123)

site_permanova <- adonis2(
  community_matrix ~ Site,
  data = community_metadata,
  method = "bray",
  permutations = 999
)

site_permanova

#PERMANOVA tests whether flower-visitor community composition
#differs between Site 1 and Site 2.


# 11. Multivariate dispersion

site_bray_dist <- vegdist(
  community_matrix,
  method = "bray"
)

site_dispersion <- betadisper(
  site_bray_dist,
  community_metadata$Site
)

anova(
  site_dispersion
)

set.seed(123)

site_dispersion_test <- permutest(
  site_dispersion,
  permutations = 999
)

site_dispersion_test

#A non-significant dispersion test supports interpreting
#the PERMANOVA as a difference in community composition
#rather than a difference in within-site dispersion.


# 12. NMDS figure

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
      "Site 1" = "#2E8B57",
      "Site 2" = "#7B3294"
    )
  ) +
  
  scale_fill_manual(
    values = c(
      "Site 1" = "#2E8B57",
      "Site 2" = "#7B3294"
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
  
  labs(
    title = "Flower visitor community composition by site",
    x = "NMDS1",
    y = "NMDS2",
    colour = "Site",
    fill = "Site"
  ) +
  
  theme_classic(
    base_size = 14
  ) +
  
  theme(
    legend.position = "right",
    
    legend.title = element_text(
      face = "bold"
    ),
    
    axis.title = element_text(
      face = "bold"
    ),
    
    panel.border = element_rect(
      fill = NA,
      colour = "black",
      linewidth = 0.7
    )
  )

nmds_plot


# 13. Flower-visitor groups by site

group_visitation_by_site <- pollinators %>%
  filter(
    Group != "None",
    Visits > 0
  ) %>%
  mutate(
    Site = factor(
      Site,
      levels = c(
        1,
        2
      ),
      labels = c(
        "Site 1",
        "Site 2"
      )
    )
  ) %>%
  group_by(
    Site,
    Group
  ) %>%
  summarise(
    Visits = sum(
      Visits,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  group_by(
    Site
  ) %>%
  mutate(
    Percent =
      100 * Visits /
      sum(Visits)
  ) %>%
  ungroup()

group_visitation_by_site


# Butterfly activity by site

butterfly_visits_by_site <- group_visitation_by_site %>%
  filter(
    Group == "Butterfly"
  )

butterfly_visits_by_site


# Hoverfly activity by site

hoverfly_visits_by_site <- group_visitation_by_site %>%
  filter(
    Group == "Hoverfly"
  )

hoverfly_visits_by_site


# 14. Floral abundance between sites

#Create one floral-abundance value for each plant population
#on each survey date to avoid counting repeated rows.

floral_abundance <- pollinators %>%
  mutate(
    Date = as.Date(
      Date,
      format = "%d/%m/%Y"
    ),
    
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
    
    Site = factor(
      Site,
      levels = c(
        1,
        2
      ),
      labels = c(
        "Site 1",
        "Site 2"
      )
    )
  ) %>%
  distinct(
    Date,
    Site,
    Plant.ID,
    Plant.Species,
    Abundance
  )


# Mean floral abundance by plant species and site

mean_floral_abundance <- floral_abundance %>%
  group_by(
    Site,
    Plant.Species
  ) %>%
  summarise(
    Mean.Abundance = mean(
      Abundance,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

mean_floral_abundance


# Thyme abundance between sites

thyme_abundance <- mean_floral_abundance %>%
  filter(
    Plant.Species ==
      "Thymus drucei"
  )

thyme_abundance


# Ratio of Site 1 to Site 2 thyme abundance

thyme_abundance_ratio <-
  thyme_abundance$Mean.Abundance[
    thyme_abundance$Site == "Site 1"
  ] /
  thyme_abundance$Mean.Abundance[
    thyme_abundance$Site == "Site 2"
  ]

thyme_abundance_ratio


# 15. Sticky catchfly between sites

#Use the finer taxonomic resolution used for descriptive visitor analyses.

visitor_data <- pollinators %>%
  filter(
    Group != "None",
    Visits > 0
  ) %>%
  mutate(
    Plant.ID = as.character(
      Plant.ID
    ),
    
    Plant.Species = unname(
      plant_key[Plant.ID]
    ),
    
    Site = case_when(
      as.character(Site) %in% c(
        "1",
        "Site 1"
      ) ~ "Site 1",
      
      as.character(Site) %in% c(
        "2",
        "Site 2"
      ) ~ "Site 2",
      
      TRUE ~ NA_character_
    ),
    
    Flower.Visitor.Taxon = case_when(
      !is.na(Species) &
        Species != "" &
        Species != "None" &
        Species != "Unknown" ~
        trimws(Species),
      
      !is.na(Genus) &
        Genus != "" ~
        paste0(
          trimws(Genus),
          " sp."
        ),
      
      !is.na(Family) &
        Family != "" ~
        paste0(
          trimws(Family),
          " family"
        ),
      
      !is.na(Order) &
        Order != "" ~
        paste0(
          trimws(Order),
          " order"
        ),
      
      TRUE ~
        as.character(Group)
    )
  )


# Sticky catchfly visitor assemblage by site

catchfly_visitors <- visitor_data %>%
  filter(
    Plant.Species ==
      "Silene viscaria"
  ) %>%
  group_by(
    Site,
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

catchfly_visitors


# Total sticky catchfly visits and visitor richness

catchfly_site_summary <- catchfly_visitors %>%
  group_by(
    Site
  ) %>%
  summarise(
    Total.Visits = sum(
      Visits
    ),
    
    Visitor.Richness = n_distinct(
      Flower.Visitor.Taxon
    ),
    
    .groups = "drop"
  )

catchfly_site_summary


# Dominant visitor at each site

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


# 16. Were sticky catchfly visitors present at the opposite site?

#All flower-visitor taxa recorded anywhere at each site

site1_taxa <- visitor_data %>%
  filter(
    Site == "Site 1"
  ) %>%
  pull(
    Flower.Visitor.Taxon
  ) %>%
  unique()

site2_taxa <- visitor_data %>%
  filter(
    Site == "Site 2"
  ) %>%
  pull(
    Flower.Visitor.Taxon
  ) %>%
  unique()


# Sticky catchfly taxa at each site

catchfly_site1_taxa <- catchfly_visitors %>%
  filter(
    Site == "Site 1"
  ) %>%
  pull(
    Flower.Visitor.Taxon
  )

catchfly_site2_taxa <- catchfly_visitors %>%
  filter(
    Site == "Site 2"
  ) %>%
  pull(
    Flower.Visitor.Taxon
  )


# Site 1 sticky catchfly visitors also occurring somewhere at Site 2

catchfly_site1_shared <- intersect(
  catchfly_site1_taxa,
  site2_taxa
)

catchfly_site1_shared

length(
  catchfly_site1_shared
)

100 *
  length(catchfly_site1_shared) /
  length(catchfly_site1_taxa)


# Site 2 sticky catchfly visitors also occurring somewhere at Site 1

catchfly_site2_shared <- intersect(
  catchfly_site2_taxa,
  site1_taxa
)

catchfly_site2_shared

length(
  catchfly_site2_shared
)

100 *
  length(catchfly_site2_shared) /
  length(catchfly_site2_taxa)


# 17. Final Aim 2 outputs

# 5. Site-model diagnostics

set.seed(123)

TimeSite_residuals <- simulateResiduals(
  fittedModel = ModelTimeSite_final,
  n = 1000
)

Uniformity_TimeSite <- testUniformity(
  TimeSite_residuals
)

Dispersion_TimeSite <- testDispersion(
  TimeSite_residuals
)

ZeroInflation_TimeSite <- testZeroInflation(
  TimeSite_residuals
)

Uniformity_TimeSite
Dispersion_TimeSite
ZeroInflation_TimeSite


# Four-panel diagnostic figure

par(
  mfrow = c(2, 2),
  mar = c(6, 5, 3, 2),
  cex = 0.75
)

# QQ residual diagnostic

plotQQunif(
  TimeSite_residuals
)

# Residuals against fitted values

plotResiduals(
  TimeSite_residuals,
  form = fitted(ModelTimeSite_final)
)

# Dispersion diagnostic

testDispersion(
  TimeSite_residuals,
  plot = TRUE
)

# Zero-inflation diagnostic

testZeroInflation(
  TimeSite_residuals,
  plot = TRUE
)

# Reset plotting layout

par(
  mfrow = c(1, 1)
)

TimeSite_Anova

Site_emmeans$emmeans
Site_emmeans$contrasts

Uniformity_TimeSite
Dispersion_TimeSite
ZeroInflation_TimeSite

site_visitation_summary
site_sampling_effort

shared_plant_species

site_nmds$stress
site_permanova
site_dispersion_test

group_visitation_by_site
butterfly_visits_by_site
hoverfly_visits_by_site

thyme_abundance
thyme_abundance_ratio

catchfly_site_summary
catchfly_dominant_visitors

catchfly_site1_shared
catchfly_site2_shared

nmds_plot

