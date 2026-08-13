# TEMPORAL DIFFERENCES
# Temporal variation in flower visitation and flower-visitor sharing


# Packages

library(dplyr)
library(tidyr)
library(ggplot2)
library(glmmTMB)
library(MASS)
library(car)
library(emmeans)
library(DHARMa)


# Load data

pollinators <- read.csv(
  "Raw Data/Pollinators.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)

#Although the dataframe and some internal object names use the term
#"pollinator", the observations represent flower visitation rather than
#confirmed pollination.


# Prepare survey-level data

survey_data <- pollinators %>%
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
  ) %>%
  
  mutate(
    Date = droplevels(
      factor(Date)
    ),
    
    Time.Period = droplevels(
      factor(Time.Period)
    ),
    
    Site = droplevels(
      factor(Site)
    ),
    
    Plant.ID = droplevels(
      factor(Plant.ID)
    )
  )


# Prepare model data

#Plant IDs 6 and 13 received no flower-visitor visits in any survey.
#They are retained in descriptive summaries but excluded from
#inferential count models.

survey_model <- survey_data %>%
  filter(
    !(Plant.ID %in% c(
      "6",
      "13"
    ))
  ) %>%
  mutate(
    Plant.ID = droplevels(
      Plant.ID
    ),
    
    Date = droplevels(
      Date
    ),
    
    Site = droplevels(
      Site
    ),
    
    Time.Period = droplevels(
      Time.Period
    )
  )


# Temporal visitation model

#The same model is used to examine both site and time-period effects.
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


# Estimated visitation by time period

Time_emmeans <- emmeans(
  ModelTimeSite_final,
  pairwise ~ Time.Period,
  type = "response",
  adjust = "tukey"
)

Time_emmeans


# Predicted visits per survey

Time_emmeans$emmeans


# Tukey-adjusted pairwise comparisons

Time_emmeans$contrasts


# Contrast comparing Time Period 2 with Time Periods 1 and 3

Time_emm <- emmeans(
  ModelTimeSite_final,
  ~ Time.Period
)

Time_period2_contrast <- contrast(
  Time_emm,
  method = list(
    "Period 2 vs Periods 1 and 3" =
      c(
        -0.5,
        1,
        -0.5
      )
  ),
  adjust = "none"
)

Time_period2_contrast


# Descriptive visitation by time period

time_visitation_summary <- survey_data %>%
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

time_visitation_summary


# Model diagnostics

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

#QQ residual diagnostic

plotQQunif(
  TimeSite_residuals
)

#Residuals against fitted values

plotResiduals(
  TimeSite_residuals,
  form = fitted(
    ModelTimeSite_final
  )
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


# Temporal visitation figure

Time_plot_data <- as.data.frame(
  emmeans(
    ModelTimeSite_final,
    ~ Time.Period,
    type = "response"
  )
)

Time_plot_data$Time.Period <- factor(
  Time_plot_data$Time.Period,
  levels = c(
    "1",
    "2",
    "3"
  ),
  labels = c(
    "Time period 1",
    "Time period 2",
    "Time period 3"
  )
)


Time_plot <- ggplot(
  Time_plot_data,
  aes(
    x = Time.Period,
    y = response,
    colour = Time.Period
  )
) +
  
  geom_errorbar(
    aes(
      ymin = asymp.LCL,
      ymax = asymp.UCL
    ),
    width = 0.15,
    linewidth = 1
  ) +
  
  geom_point(
    size = 4.5
  ) +
  
  scale_colour_manual(
    values = c(
      "Time period 1" = "#CC0033",
      "Time period 2" = "#FF9900",
      "Time period 3" = "#FFCC00"
    )
  ) +
  
  labs(
    x = "Time period",
    y = "Estimated floral visitor visits per 20-minute survey"
  ) +
  
  theme_classic(
    base_size = 14
  ) +
  
  theme(
    legend.position = "none",
    
    axis.title = element_text(
      face = "bold"
    ),
    
    axis.text = element_text(
      colour = "black"
    ),
    
    panel.border = element_rect(
      fill = NA,
      colour = "black",
      linewidth = 0.7
    )
  )

Time_plot


# Final outputs

TimeSite_Anova

Time_emmeans$emmeans
Time_emmeans$contrasts

Time_period2_contrast

time_visitation_summary

Uniformity_TimeSite
Dispersion_TimeSite
ZeroInflation_TimeSite

Time_plot
