##################################
#Site Comparison
##################################

# Aim 2:
# Determine the extent to which pollinator visitation
# and pollinator community composition vary between sites.


source("Overall Survey and Data.R")

source("Visitation Time and Site.R")


#Site visitation

# Site visitation is tested in the shared Time.Period + Site GLMM
# fitted in Visitation Time and Site.R.

# Display the overall Site effect.

TimeSite_Anova


# Display estimated visitation at each site.

Site_emmeans$emmeans


# Display the response-scale comparison between sites.

Site_emmeans$contrasts


#Descriptive site summaries

site_descriptive_summary <- survey_data %>%
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

site_descriptive_summary


#Prepare pollinator community data

# Species-level identity is retained for bumblebees and butterflies.
# Broader pollinator groups are used for all other taxa.

community_data <- pollinators %>%
  filter(
    Group != "None",
    Visits > 0
  ) %>%
  
  mutate(
    Group = trimws(
      as.character(Group)
    ),
    
    Species = trimws(
      as.character(Species)
    ),
    
    Pollinator.Taxon = case_when(
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
    
    # Create a unique sampling-occasion ID.
    
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


#Check pollinator taxa used in the analysis

sort(
  unique(
    community_data$Pollinator.Taxon
  )
)


#Build community matrix

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


#Create numeric community matrix

community_matrix <- community_matrix_data %>%
  
  dplyr::select(
    -Sample.ID,
    -Site
  ) %>%
  
  as.data.frame()


rownames(
  community_matrix
) <- community_matrix_data$Sample.ID


#Check community matrix

dim(
  community_matrix
)

table(
  community_metadata$Site
)

head(
  community_matrix
)

rowSums(
  community_matrix
)

sum(
  rowSums(
    community_matrix
  ) == 0
)


# No row should have zero total pollinator abundance.

if (
  any(
    rowSums(
      community_matrix
    ) == 0
  )
) {
  
  stop(
    "At least one sampling occasion contains no pollinator visits."
  )
}


#Run NMDS

set.seed(
  123
)

site_nmds <- metaMDS(
  community_matrix,
  distance = "bray",
  k = 2,
  trymax = 200,
  autotransform = FALSE,
  trace = FALSE
)

site_nmds


#Check NMDS stress

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


#Calculate site centroids

site_centroids <- nmds_scores %>%
  group_by(
    Site
  ) %>%
  
  summarise(
    NMDS1 = mean(
      NMDS1
    ),
    
    NMDS2 = mean(
      NMDS2
    ),
    
    .groups = "drop"
  )


#Create stress label

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

# PERMANOVA tests whether pollinator community composition
# differs significantly between Site 1 and Site 2.


#Check homogeneity of multivariate dispersion

bray_dist <- vegdist(
  community_matrix,
  method = "bray"
)

site_dispersion <- betadisper(
  bray_dist,
  community_metadata$Site
)


site_dispersion_anova <- anova(
  site_dispersion
)

site_dispersion_permutation <- permutest(
  site_dispersion,
  permutations = 999
)

site_dispersion_anova

site_dispersion_permutation

# A non-significant dispersion test supports interpreting
# the PERMANOVA as a test of differences in community composition
# rather than differences in within-site dispersion.


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


#Final Aim 2 outputs

Site_emmeans$emmeans

Site_emmeans$contrasts

site_descriptive_summary

site_nmds$stress

site_permanova

site_dispersion_anova

site_dispersion_permutation

nmds_plot
