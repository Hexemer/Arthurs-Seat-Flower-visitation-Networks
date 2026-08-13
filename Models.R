#Import Data

pollinators <- read.csv(
  "Pollinators.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)

# Rename pollinator groups

pollinators$Group[
  trimws(pollinators$Group) == "Other Fly"
] <- "Non-syrphid Fly"

community_data$Sample.ID <- as.character(community_data$Sample.ID)

community_data$Site[community_data$Pollinator.ID == 16] <- 1

community_data$Sample.ID[community_data$Pollinator.ID == 16] <-
  "04/06/2026.2.1.1.1.5"

community_data %>%
  filter(Pollinator.ID == 16)

#Check required columns

required_columns <- c(
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


#Prep variables

pollinators$Group <- factor(trimws(pollinators$Group))
pollinators$Plant.ID <- factor(pollinators$Plant.ID)
pollinators$Time.Period <- factor(pollinators$Time.Period)
pollinators$Site <- factor(pollinators$Site)
pollinators$Block <- factor(pollinators$Block)
pollinators$Replication <- factor(pollinators$Replication)

analysis_data <- droplevels(
  subset(pollinators, Group != "None")
)

analysis_data$Group <- relevel(
  analysis_data$Group,
  ref = "Bumblebee"
)


#Load negative binomial functions

library(MASS)

#Block or no block

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

###############################################################
# QUESTIONS 1 AND 2
#
# 1. Does pollinator visitation vary throughout the day?
# 2. Does pollinator visitation vary between sites?
###############################################################


# =============================================================
# 1. LOAD PACKAGES
# =============================================================

library(dplyr)
library(glmmTMB)
library(car)
library(emmeans)
library(DHARMa)


# Do not include install.packages() in the final analysis script.
# Packages only need to be installed once.


# =============================================================
# 2. CREATE ONE ROW PER 20-MINUTE PLANT SURVEY
# =============================================================

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


# =============================================================
# 3. CHECK THE SURVEY-LEVEL DATA
# =============================================================

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


# =============================================================
# 4. PREPARE THE INFERENTIAL DATASET
# =============================================================

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


# =============================================================
# 5. FIT THE TIME-PERIOD AND SITE MODEL
# =============================================================

# Date accounts for non-independence among surveys conducted
# on the same sampling day.
#
# Plant.ID accounts for repeated surveys of the same plant
# population and differences in baseline attractiveness.

ModelTimeSite <- glmmTMB(
  Visits ~
    Time.Period +
    Site +
    (1 | Date) +
    (1 | Plant.ID),
  family = nbinom2,
  data = survey_model
)


summary(ModelTimeSite)


# =============================================================
# 6. TEST THE OVERALL FIXED EFFECTS
# =============================================================

TimeSite_Anova <- car::Anova(
  ModelTimeSite,
  type = 2
)

TimeSite_Anova


# Read the relevant rows as:
#
# Time.Period:
# Does visitation vary throughout the day?
#
# Site:
# Does visitation differ between the two sites?


# =============================================================
# 7. QUESTION 1: ESTIMATED VISITATION BY TIME PERIOD
# =============================================================

Time_emmeans <- emmeans(
  ModelTimeSite,
  pairwise ~ Time.Period,
  type = "response",
  adjust = "tukey"
)

Time_emmeans

# Predicted visits per survey for each period
Time_emmeans$emmeans

# Tukey-adjusted pairwise comparisons
Time_emmeans$contrasts


# Optional compact letter display
# install.packages("multcomp") only if not already installed

# library(multcomp)
# cld(
#   emmeans(
#     ModelTimeSite,
#     ~ Time.Period,
#     type = "response"
#   ),
#   adjust = "tukey"
# )


# =============================================================
# 8. QUESTION 2: ESTIMATED VISITATION BY SITE
# =============================================================

Site_emmeans <- emmeans(
  ModelTimeSite,
  pairwise ~ Site,
  type = "response"
)

Site_emmeans

# Predicted visits per survey for each site
Site_emmeans$emmeans

# Site comparison as a response-scale ratio
Site_emmeans$contrasts


# =============================================================
# 9. DESCRIPTIVE SUMMARIES
# =============================================================

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


# =============================================================
# 10. CHECK THE RANDOM EFFECTS AND CONVERGENCE
# =============================================================

VarCorr(ModelTimeSite)

ModelTimeSite$sdr$pdHess

# TRUE indicates a positive-definite Hessian and supports
# successful convergence.


# =============================================================
# 11. MODEL DIAGNOSTICS
# =============================================================

simulationOutput_TimeSite <- simulateResiduals(
  fittedModel = ModelTimeSite,
  n = 1000
)

plot(simulationOutput_TimeSite)

Uniformity_TimeSite <- testUniformity(
  simulationOutput_TimeSite
)

Dispersion_TimeSite <- testDispersion(
  simulationOutput_TimeSite
)

ZeroInflation_TimeSite <- testZeroInflation(
  simulationOutput_TimeSite
)

Uniformity_TimeSite
Dispersion_TimeSite
ZeroInflation_TimeSite


# =============================================================
# 12. OPTIONAL: COMPARE NBINOM1 AND NBINOM2
# =============================================================

ModelTimeSite_nbinom1 <- update(
  ModelTimeSite,
  family = nbinom1
)

ModelTimeSite_nbinom1$sdr$pdHess

AIC(
  ModelTimeSite,
  ModelTimeSite_nbinom1
)

# Only select nbinom1 if:
# 1. it converges cleanly;
# 2. pdHess is TRUE;
# 3. diagnostics are satisfactory;
# 4. AIC is meaningfully lower.


# =============================================================
# 13. OPTIONAL: CHECK THE VALUE OF THE DATE RANDOM EFFECT
# =============================================================

ModelTimeSite_no_date <- glmmTMB(
  Visits ~
    Time.Period +
    Site +
    (1 | Plant.ID),
  family = nbinom2,
  data = survey_model
)

AIC(
  ModelTimeSite_no_date,
  ModelTimeSite
)

# Date should normally remain because surveys from the same
# sampling date are not independent. This comparison is only
# a supporting check.


# =============================================================
# 14. FINAL OUTPUTS TO SAVE
# =============================================================

summary(ModelTimeSite)

TimeSite_Anova

Time_emmeans$emmeans
Time_emmeans$contrasts

Site_emmeans$emmeans
Site_emmeans$contrasts

VarCorr(ModelTimeSite)

AIC(ModelTimeSite)

ModelTimeSite$sdr$pdHess

simulationOutput_TimeSite_nbinom1 <- simulateResiduals(
  fittedModel = ModelTimeSite_nbinom1,
  n = 1000
)

plot(simulationOutput_TimeSite_nbinom1)

testUniformity(simulationOutput_TimeSite_nbinom1)
testDispersion(simulationOutput_TimeSite_nbinom1)
testZeroInflation(simulationOutput_TimeSite_nbinom1)

summary(ModelTimeSite_nbinom1)

car::Anova(
  ModelTimeSite_nbinom1,
  type = 2
)
summary(ModelTimeSite_nbinom1)

car::Anova(ModelTimeSite_nbinom1, type = 2)

Time_emmeans <- emmeans(
  ModelTimeSite_nbinom1,
  pairwise ~ Time.Period,
  type = "response"
)

Time_emmeans$emmeans

Site_emmeans <- emmeans(
  ModelTimeSite_nbinom1,
  pairwise ~ Site,
  type = "response"
)

Site_emmeans$emmeans

###########################################################################################
# 3. Does the effect of temperature on pollinator visitation vary among pollinator groups?
###########################################################################################

# Load packages
library(MASS)
library(car)
library(emmeans)
library(ggplot2)
library(dplyr)

# ============================================================
# PREPARE DATA
# ============================================================

# Check the number of observations in each pollinator group
table(analysis_data$Group)

# Remove groups with too few observations for the main analysis
analysis_data_temp <- droplevels(
  subset(
    analysis_data,
    !(Group %in% c("Mining Bee", "Wasp"))
  )
)

table(analysis_data_temp$Group)

# Mean-centre temperature
mean_temperature <- mean(
  analysis_data_temp$Temperature..C.degrees.,
  na.rm = TRUE
)

analysis_data_temp <- analysis_data_temp |>
  mutate(
    Temperature_c =
      Temperature..C.degrees. - mean_temperature
  )

# ============================================================
# LINEAR TEMPERATURE × GROUP MODEL
# ============================================================

ModelTempGroup_linear <- glm.nb(
  Visits ~
    Temperature_c * Group +
    Site +
    Time.Period +
    Plant.ID,
  data = analysis_data_temp
)

summary(ModelTempGroup_linear)

# Test the overall effects and interaction
TempGroup_linear_Anova <- car::Anova(
  ModelTempGroup_linear,
  type = 2
)

TempGroup_linear_Anova

# Estimate the temperature slope for each pollinator group
Temp_slopes <- emtrends(
  ModelTempGroup_linear,
  ~ Group,
  var = "Temperature_c"
)

Temp_slopes

# Pairwise comparisons among group-specific temperature slopes
Temp_slope_pairs <- pairs(
  Temp_slopes,
  adjust = "tukey"
)

Temp_slope_pairs

# ============================================================
# QUADRATIC TEMPERATURE × GROUP MODEL
# ============================================================

# Sweat bees are excluded because they contain only seven observations
# spanning two unique temperature values, which is insufficient to
# estimate a group-specific quadratic response.

analysis_data_quad <- droplevels(
  subset(
    analysis_data_temp,
    Group != "Sweat Bee"
  )
)

table(analysis_data_quad$Group)

# Re-centre temperature using the quadratic-analysis dataset
mean_temperature_quad <- mean(
  analysis_data_quad$Temperature..C.degrees.,
  na.rm = TRUE
)

analysis_data_quad <- analysis_data_quad |>
  mutate(
    Temperature_c =
      Temperature..C.degrees. - mean_temperature_quad
  )

# Fit a linear model to the same reduced dataset
ModelTempGroup_linear_quad_data <- glm.nb(
  Visits ~
    Temperature_c * Group +
    Site +
    Time.Period +
    Plant.ID,
  data = analysis_data_quad
)

# Fit the quadratic interaction model
ModelTempGroup_quad <- glm.nb(
  Visits ~
    (
      Temperature_c +
        I(Temperature_c^2)
    ) * Group +
    Site +
    Time.Period +
    Plant.ID,
  data = analysis_data_quad
)

summary(ModelTempGroup_quad)

# Test the effects in the quadratic model
TempGroup_quad_Anova <- car::Anova(
  ModelTempGroup_quad,
  type = 2
)

TempGroup_quad_Anova

# ============================================================
# COMPARE LINEAR AND QUADRATIC MODELS
# ============================================================

# Likelihood ratio test
Temp_model_comparison <- anova(
  ModelTempGroup_linear_quad_data,
  ModelTempGroup_quad,
  test = "Chisq"
)

Temp_model_comparison

# Compare AIC values
AIC(
  ModelTempGroup_linear_quad_data,
  ModelTempGroup_quad
)

# Retain the quadratic model only if it significantly improves model fit
# and has a meaningfully lower AIC.

# ============================================================
# CHECK MODEL DIAGNOSTICS
# ============================================================

library(DHARMa)

# Linear model diagnostics
Linear_residuals <- simulateResiduals(
  fittedModel = ModelTempGroup_linear
)

plot(Linear_residuals)

testUniformity(Linear_residuals)
testDispersion(Linear_residuals)
testZeroInflation(Linear_residuals)

# Quadratic model diagnostics
Quadratic_residuals <- simulateResiduals(
  fittedModel = ModelTempGroup_quad
)

plot(Quadratic_residuals)

testUniformity(Quadratic_residuals)
testDispersion(Quadratic_residuals)
testZeroInflation(Quadratic_residuals)

# ============================================================
# CREATE PREDICTIONS FROM THE QUADRATIC MODEL
# ============================================================

temperature_sequence <- seq(
  min(
    analysis_data_quad$Temperature..C.degrees.,
    na.rm = TRUE
  ),
  max(
    analysis_data_quad$Temperature..C.degrees.,
    na.rm = TRUE
  ),
  length.out = 100
)

# Convert the original temperature sequence to centred temperature
temperature_sequence_c <-
  temperature_sequence - mean_temperature_quad

Temp_predictions_quad <- emmeans(
  ModelTempGroup_quad,
  ~ Group | Temperature_c,
  at = list(
    Temperature_c = temperature_sequence_c
  ),
  nuisance = c(
    "Site",
    "Time.Period",
    "Plant.ID"
  ),
  type = "response"
)

Temp_predictions_quad_df <- as.data.frame(
  Temp_predictions_quad
)

# Convert centred temperature back to degrees Celsius for plotting
Temp_predictions_quad_df <- Temp_predictions_quad_df |>
  mutate(
    Temperature = Temperature_c + mean_temperature_quad
  )

head(Temp_predictions_quad_df)

# ============================================================
# PREDICT FINAL TEMPERATURE MODEL WITHIN OBSERVED RANGES
# ============================================================

library(emmeans)
library(dplyr)
library(ggplot2)

# Obtain the observed temperature range for each pollinator group
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


# Generate predictions separately for each group
Temp_predictions_final_df <- lapply(
  seq_len(nrow(group_temperature_ranges)),
  function(i) {
    
    current_group <- as.character(
      group_temperature_ranges$Group[i]
    )
    
    temperature_values <- seq(
      from = group_temperature_ranges$min_temperature[i],
      to = group_temperature_ranges$max_temperature[i],
      length.out = 100
    )
    
    # Convert original temperatures to centred temperatures
    centred_values <- temperature_values -
      mean_temperature_quad
    
    predictions <- emmeans(
      ModelTempGroup_final,
      ~ Temperature_c,
      at = list(
        Temperature_c = centred_values,
        Group = current_group
      ),
      type = "response"
    )
    
    predictions_df <- as.data.frame(predictions)
    
    predictions_df %>%
      mutate(
        Group = current_group,
        
        # Convert centred temperature back to degrees Celsius
        Temperature =
          Temperature_c + mean_temperature_quad
      )
  }
) %>%
  bind_rows()


# Ensure the group order is retained
Temp_predictions_final_df$Group <- factor(
  Temp_predictions_final_df$Group,
  levels = levels(analysis_data_quad$Group)
)



# ============================================================
# FINAL TEMPERATURE MODEL
# ============================================================

# The quadratic interaction between temperature and pollinator group
# did not significantly improve model fit, so a simpler model with
# a shared quadratic response was retained.

ModelTempGroup_final <- glm.nb(
  Visits ~
    Temperature_c +
    I(Temperature_c^2) +
    Group +
    Temperature_c:Group +
    Site +
    Time.Period +
    Plant.ID,
  data = analysis_data_quad
)

summary(ModelTempGroup_final)

# Confirm that removing the quadratic interaction does not
# significantly reduce model fit

anova(
  ModelTempGroup_final,
  ModelTempGroup_quad,
  test = "Chisq"
)

AIC(
  ModelTempGroup_final,
  ModelTempGroup_quad
)

car::Anova(
  ModelTempGroup_final,
  type = 2
)

# ============================================================
# PREDICT FINAL TEMPERATURE MODEL WITHIN OBSERVED RANGES
# ============================================================

library(emmeans)
library(dplyr)
library(ggplot2)

# Obtain the observed temperature range for each pollinator group
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


# Generate predictions separately for each group
Temp_predictions_final_df <- lapply(
  seq_len(nrow(group_temperature_ranges)),
  function(i) {
    
    current_group <- as.character(
      group_temperature_ranges$Group[i]
    )
    
    temperature_values <- seq(
      from = group_temperature_ranges$min_temperature[i],
      to = group_temperature_ranges$max_temperature[i],
      length.out = 100
    )
    
    # Convert original temperatures to centred temperatures
    centred_values <- temperature_values -
      mean_temperature_quad
    
    predictions <- emmeans(
      ModelTempGroup_final,
      ~ Temperature_c,
      at = list(
        Temperature_c = centred_values,
        Group = current_group
      ),
      type = "response"
    )
    
    predictions_df <- as.data.frame(predictions)
    
    predictions_df %>%
      mutate(
        Group = current_group,
        
        # Convert centred temperature back to degrees Celsius
        Temperature =
          Temperature_c + mean_temperature_quad
      )
  }
) %>%
  bind_rows()


# Ensure the group order is retained
Temp_predictions_final_df$Group <- factor(
  Temp_predictions_final_df$Group,
  levels = levels(analysis_data_quad$Group)
)



# ============================================================
# PREDICTIONS FROM THE FINAL SIMPLIFIED QUADRATIC MODEL
# ============================================================

library(dplyr)
library(ggplot2)

# Check that the final model exists
if (!exists("ModelTempGroup_final")) {
  stop("ModelTempGroup_final does not exist. Run the final model first.")
}

# Calculate the observed temperature range for each group
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

# Create 100 temperature values within each group's observed range
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
    # Convert temperature to the centred scale used in the model
    Temperature_c = Temperature - mean_temperature_quad,
    
    # Set other predictors to reference levels
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
    
    Group = factor(
      Group,
      levels = levels(analysis_data_quad$Group)
    )
  )

# Generate predictions on the model's log scale
prediction_results <- predict(
  ModelTempGroup_final,
  newdata = Temp_predictions_final_df,
  type = "link",
  se.fit = TRUE
)

# Convert predictions and confidence intervals to the response scale
Temp_predictions_final_df <- Temp_predictions_final_df %>%
  mutate(
    response = exp(prediction_results$fit),
    lower_CL = exp(
      prediction_results$fit -
        1.96 * prediction_results$se.fit
    ),
    upper_CL = exp(
      prediction_results$fit +
        1.96 * prediction_results$se.fit
    )
  )

# Check that the prediction data were created
head(Temp_predictions_final_df)

# ============================================================
# PLOT FINAL QUADRATIC TEMPERATURE RESPONSES
# ============================================================

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
    y = response,
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
  labs(
    x = "Temperature (°C)",
    y = "Number of visits per survey",
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
  ) +
  scale_colour_manual(
    values = pollinator_colours
  ) +
  scale_fill_manual(
    values = pollinator_colours
  )

Temperature_plot

Temperature_plot

#######################################################################
#Temperature with date as a random effect
#######################################################################
################################################################################
# QUESTION 3
# Does the effect of temperature on pollinator visitation vary among
# pollinator groups?
#
# Negative-binomial GLMM with survey Date as a random intercept
################################################################################


# ============================================================
# 1. LOAD PACKAGES
# ============================================================

library(dplyr)
library(glmmTMB)
library(car)
library(DHARMa)
library(ggplot2)


# ============================================================
# 2. CORRECT DATES IN THE SOURCE DATA, IF NOT ALREADY DONE
# ============================================================

# Only retain this section if analysis_data still contains the incorrect dates.

analysis_data$Date[
  analysis_data$Date == as.Date("2026-06-06")
] <- as.Date("2026-06-05")

analysis_data$Date[
  analysis_data$Date == as.Date("2026-07-03")
] <- as.Date("2026-07-02")


# ============================================================
# 3. PREPARE DATA FOR THE LINEAR ANALYSIS
# ============================================================

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


# Check group and date sample sizes
table(analysis_data_temp$Group)
table(analysis_data_temp$Date)


# Mean-centre temperature
mean_temperature <- mean(
  analysis_data_temp$Temperature..C.degrees.,
  na.rm = TRUE
)

analysis_data_temp <- analysis_data_temp %>%
  mutate(
    Temperature_c =
      Temperature..C.degrees. - mean_temperature
  )


# ============================================================
# 4. LINEAR TEMPERATURE × POLLINATOR-GROUP MODEL
# ============================================================

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

car::Anova(
  ModelTempGroup_linear,
  type = 2
)


# ============================================================
# 5. PREPARE DATA FOR THE QUADRATIC ANALYSIS
# ============================================================

# Sweat bees are excluded from the quadratic analysis because
# their observations span too few distinct temperatures to
# estimate a group-specific quadratic curve reliably.

analysis_data_quad <- analysis_data_temp %>%
  filter(
    Group != "Sweat Bee"
  ) %>%
  droplevels()


table(analysis_data_quad$Group)
table(analysis_data_quad$Date)


# Re-centre temperature using this reduced dataset
mean_temperature_quad <- mean(
  analysis_data_quad$Temperature..C.degrees.,
  na.rm = TRUE
)

analysis_data_quad <- analysis_data_quad %>%
  mutate(
    Temperature_c =
      Temperature..C.degrees. - mean_temperature_quad
  )


# ============================================================
# 6. LINEAR MODEL ON THE SAME QUADRATIC-ANALYSIS DATASET
# ============================================================

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


# ============================================================
# 7. FULL GROUP-SPECIFIC QUADRATIC MODEL
# ============================================================

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

car::Anova(
  ModelTempGroup_quad,
  type = 2
)


# ============================================================
# 8. DOES A GROUP-SPECIFIC QUADRATIC RESPONSE IMPROVE FIT?
# ============================================================

anova(
  ModelTempGroup_linear_quad_data,
  ModelTempGroup_quad
)

AIC(
  ModelTempGroup_linear_quad_data,
  ModelTempGroup_quad
)


# ============================================================
# 9. SIMPLIFIED MODEL
#
# Shared quadratic temperature response, but group-specific
# linear temperature slopes.
# ============================================================

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

car::Anova(
  ModelTempGroup_final,
  type = 2
)


# Compare the simplified model with the full quadratic interaction
anova(
  ModelTempGroup_final,
  ModelTempGroup_quad
)

AIC(
  ModelTempGroup_final,
  ModelTempGroup_quad
)


# Retain ModelTempGroup_final if removing the Group × Temperature²
# interaction does not significantly worsen model fit.


# ============================================================
# 10. TEST THE OVERALL LINEAR TEMPERATURE × GROUP INTERACTION
# ============================================================

# This reduced model removes group-specific linear temperature slopes
# while retaining the shared quadratic temperature response.

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


# Interpretation:
# A significant likelihood-ratio test means that the linear
# temperature effect differs among pollinator groups.


# ============================================================
# 11. CHECK DATE RANDOM EFFECT
# ============================================================

VarCorr(ModelTempGroup_final)

ModelTempGroup_final$sdr$pdHess


# Do not automatically remove Date. This comparison is only a
# supplementary check.

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


# ============================================================
# 12. MODEL DIAGNOSTICS
# ============================================================

TempGroup_residuals <- simulateResiduals(
  fittedModel = ModelTempGroup_final,
  n = 1000
)

plot(TempGroup_residuals)

testUniformity(TempGroup_residuals)
testDispersion(TempGroup_residuals)
testZeroInflation(TempGroup_residuals)


# ============================================================
# 13. CHECK NBINOM1 VERSUS NBINOM2
# ============================================================

ModelTempGroup_final_nbinom1 <- update(
  ModelTempGroup_final,
  family = nbinom1
)

ModelTempGroup_final_nbinom1$sdr$pdHess

AIC(
  ModelTempGroup_final,
  ModelTempGroup_final_nbinom1
)

# Only prefer nbinom1 if it converges cleanly and its DHARMa
# diagnostics are satisfactory.


# ============================================================
# 14. GENERATE PREDICTIONS WITHIN EACH GROUP'S OBSERVED RANGE
# ============================================================

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
    
    # Any valid date is required in newdata. Setting re.form = NA
    # below excludes the date-specific random effect.
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


# ============================================================
# 15. PLOT TEMPERATURE RESPONSES
# ============================================================

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

#######################################################################
# 4. Does visitation differ among plant species after accounting
#    for floral abundance?
#######################################################################

library(MASS)
library(dplyr)
library(emmeans)
library(ggplot2)


# ============================================================
# CREATE PLANT SPECIES VARIABLE
# ============================================================

# Convert Plant.ID to character first so the mapping works
# whether Plant.ID is currently numeric or a factor
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


# ============================================================
# PREPARE ANALYSIS DATA
# ============================================================

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


# ============================================================
# FIT SPECIES MODEL
# ============================================================

# The offset accounts for the number of flowers available
ModelSpecies_nb <- glm.nb(
  Visits ~
    Plant.Species +
    Site +
    Time.Period +
    offset(log(Number.of.Flowers)),
  data = analysis_data_plant
)


# ============================================================
# FIT NULL MODEL WITHOUT PLANT SPECIES
# ============================================================

ModelAbundance_nb <- glm.nb(
  Visits ~
    Site +
    Time.Period +
    offset(log(Number.of.Flowers)),
  data = analysis_data_plant
)


# ============================================================
# COMPARE MODELS
# ============================================================

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
# 5. Which plant species are most attractive to pollinators?
#######################################################################

# ============================================================
# ESTIMATED VISITATION RATES
# ============================================================

# offset = 0 corresponds to log(1 flower)
# Therefore, estimates are predicted visits per flower per survey
Species_emmeans <- emmeans(
  ModelSpecies_nb,
  ~ Plant.Species,
  type = "response",
  offset = 0
)

Species_emmeans


# ============================================================
# CREATE RESULTS TABLE
# ============================================================

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
#6. Which environmental variables (temperature, solar radiation, wind speed, humidity) best predict pollinator visitation?
#######################################################################################################################
# ============================================================
# QUESTION 6
# Which environmental variables best predict pollinator
# visitation?
#
# Response: visits
# Offset: number of flowers
# Fixed controls: plant identity
# Random intercept: survey date
# Error distribution: negative binomial
# ============================================================

#Fix dates

survey_environment$Date[
  survey_environment$Date == as.Date("2026-06-06")
] <- as.Date("2026-06-05")

survey_environment$Date[
  survey_environment$Date == as.Date("2026-07-03")
] <- as.Date("2026-07-02")

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
      mean(Temperature..C.degrees., na.rm = TRUE),
    
    Solar_c = Solar.Power..W.m.2. -
      mean(Solar.Power..W.m.2., na.rm = TRUE),
    
    Humidity_c = Humidity....RH. -
      mean(Humidity....RH., na.rm = TRUE),
    
    Wind_c = Wind.Speed..m.s. -
      mean(Wind.Speed..m.s., na.rm = TRUE)
  )

table(survey_environment2$Date)

# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------
install.packages("glmmTMB")


library(dplyr)
library(glmmTMB)
library(car)
library(DHARMa)


# ------------------------------------------------------------
# 2. Prepare the data
# ------------------------------------------------------------

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
      mean(Temperature..C.degrees., na.rm = TRUE),
    
    Solar_c = Solar.Power..W.m.2. -
      mean(Solar.Power..W.m.2., na.rm = TRUE),
    
    Humidity_c = Humidity....RH. -
      mean(Humidity....RH., na.rm = TRUE),
    
    Wind_c = Wind.Speed..m.s. -
      mean(Wind.Speed..m.s., na.rm = TRUE)
  )


# Check the final data
nrow(survey_environment2)
table(survey_environment2$Date)
table(survey_environment2$Plant.ID)
summary(survey_environment2$Number.of.Flowers)


# ------------------------------------------------------------
# 3. Linear abiotic model
#
# Plant.ID controls for differences in attractiveness among
# plants. Date accounts for observations collected on the
# same survey day.
# ------------------------------------------------------------

Environment_linear <- glmmTMB(
  Visits ~
    Temp_c +
    Solar_c +
    Humidity_c +
    Wind_c +
    Plant.ID +
    offset(log(Number.of.Flowers)) +
    (1 | Date),
  family = nbinom2,
  data = survey_environment2
)

summary(Environment_linear)


# ------------------------------------------------------------
# 4. Test quadratic temperature response
# ------------------------------------------------------------

Environment_temp_quad <- glmmTMB(
  Visits ~
    Temp_c +
    I(Temp_c^2) +
    Solar_c +
    Humidity_c +
    Wind_c +
    Plant.ID +
    offset(log(Number.of.Flowers)) +
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


# ------------------------------------------------------------
# 5. Test quadratic solar-radiation response
# ------------------------------------------------------------

Environment_solar_quad <- glmmTMB(
  Visits ~
    Temp_c +
    Solar_c +
    I(Solar_c^2) +
    Humidity_c +
    Wind_c +
    Plant.ID +
    offset(log(Number.of.Flowers)) +
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


# ------------------------------------------------------------
# 6. Test quadratic humidity response
# ------------------------------------------------------------

Environment_humidity_quad <- glmmTMB(
  Visits ~
    Temp_c +
    Solar_c +
    Humidity_c +
    I(Humidity_c^2) +
    Wind_c +
    Plant.ID +
    offset(log(Number.of.Flowers)) +
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


# ------------------------------------------------------------
# 7. Candidate final climate model
#
# Keep temperature and solar quadratic terms if supported by
# the comparisons above.
# ------------------------------------------------------------

Environment_final2 <- glmmTMB(
  Visits ~
    Temp_c +
    I(Temp_c^2) +
    Solar_c +
    Humidity_c +
    Wind_c +
    Plant.ID +
    offset(log(Number.of.Flowers)) +
    (1 | Date),
  family = nbinom2,
  data = survey_environment2
)

summary(Environment_final2)

car::Anova(
  Environment_final2,
  type = 2
)

anova(Environment_final2,
      Environment_final)

AIC(Environment_final2,
    Environment_final)

# ------------------------------------------------------------
# 8. Does time period explain anything beyond climate?
# ------------------------------------------------------------

Environment_final_time <- glmmTMB(
  Visits ~
    Temp_c +
    I(Temp_c^2) +
    Solar_c +
    I(Solar_c^2) +
    Humidity_c +
    Wind_c +
    Plant.ID +
    Time.Period +
    offset(log(Number.of.Flowers)) +
    (1 | Date),
  family = nbinom2,
  data = survey_environment2
)

anova(
  Environment_final,
  Environment_final_time
)

AIC(
  Environment_final,
  Environment_final_time
)

survey_environment2 %>%
  filter(Date %in% c("2026-06-06", "2026-07-03")) %>%
  select(
    Date,
    Site,
    Plant.ID,
    Visits,
    Number.of.Flowers,
    Time.Period
  )

# ------------------------------------------------------------
# 9. Check whether Date explains meaningful variation
# ------------------------------------------------------------

VarCorr(Environment_final)

Environment_final_no_date <- glmmTMB(
  Visits ~
    Temp_c +
    I(Temp_c^2) +
    Solar_c +
    I(Solar_c^2) +
    Humidity_c +
    Wind_c +
    Plant.ID +
    offset(log(Number.of.Flowers)),
  family = nbinom2,
  data = survey_environment2
)

AIC(
  Environment_final_no_date,
  Environment_final
)

# With only six dates, use the variance estimate, model
# stability and AIC as the main evidence for retaining Date.


# ------------------------------------------------------------
# 10. Model diagnostics
# ------------------------------------------------------------

Environment_diagnostics <- simulateResiduals(
  Environment_final,
  n = 1000
)

plot(Environment_diagnostics)

testDispersion(Environment_diagnostics)
testZeroInflation(Environment_diagnostics)
testUniformity(Environment_diagnostics)


# ------------------------------------------------------------
# 11. Check convergence
# ------------------------------------------------------------

Environment_final$sdr$pdHess

# TRUE indicates a positive-definite Hessian and supports
# successful model convergence.


# ------------------------------------------------------------
# 12. Compare nbinom1 and nbinom2 if necessary
# ------------------------------------------------------------

Environment_final_nbinom1 <- update(
  Environment_final,
  family = nbinom1
)

AIC(
  Environment_final,
  Environment_final_nbinom1
)


# ------------------------------------------------------------
# 13. Final outputs
# ------------------------------------------------------------

summary(Environment_final)

car::Anova(
  Environment_final,
  type = 2
)

AIC(Environment_final)

VarCorr(Environment_final)

########################################################################
#Does pollinator community composition differ between Site 1 and Site 2?
########################################################################

library(dplyr)
library(tidyr)
library(vegan)

community_data <- pollinators %>%
  filter(
    Group != "None",
    Visits > 0
  ) %>%
  mutate(
    Group = trimws(Group),
    Species = trimws(Species),
    
    Pollinator.Taxon = case_when(
      Group %in% c("Bumblebee", "Butterfly") &
        !is.na(Species) &
        Species != "" ~ Species,
      
      TRUE ~ Group
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
sort(unique(community_data$Pollinator.Taxon))

#Build the community matrix

community_matrix_data <- community_data %>%
  group_by(
    Sample.ID,
    Site,
    Pollinator.Taxon
  ) %>%
  summarise(
    Visits = sum(Visits),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Pollinator.Taxon,
    values_from = Visits,
    values_fill = 0
  )
#Create the meta data
community_metadata <- community_matrix_data %>%
  select(
    Sample.ID,
    Site
  )
community_matrix <- community_matrix_data %>%
  select(
    -Sample.ID,
    -Site
  ) %>%
  as.data.frame()

rownames(community_matrix) <- community_matrix_data$Sample.ID

dim(community_matrix)

table(community_metadata$Site)

head(community_matrix)

rowSums(community_matrix) == 0

sum(rowSums(community_matrix) == 0)

#NMDS

install.packages("vegan")
library(vegan)

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

#Check the stress
site_nmds$stress

#Extract the scores
nmds_scores <- as.data.frame(scores(site_nmds, display = "sites"))

nmds_scores$Site <- community_metadata$Site

# Make Site a factor with nice labels
nmds_scores$Site <- factor(
  community_metadata$Site,
  levels = c(1, 2),
  labels = c("Site 1", "Site 2")
)

# Calculate the centroid of each site
site_centroids <- aggregate(
  cbind(NMDS1, NMDS2) ~ Site,
  data = nmds_scores,
  FUN = mean
)

# Create the stress label
stress_label <- paste0(
  "Stress = ",
  round(site_nmds$stress, 3)
)

#Plot the ordination
library(ggplot2)

ggplot(nmds_scores,
       aes(x = NMDS1,
           y = NMDS2,
           colour = Site)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(group = Site), level = 0.95) +
  theme_classic() +
  labs(
    title = "Pollinator community composition by site",
    x = "NMDS1",
    y = "NMDS2",
    colour = "Site"
  )

#Run the statistical test
site_permanova <- adonis2(
  community_matrix ~ Site,
  data = community_metadata,
  method = "bray",
  permutations = 999
)

site_permanova

bray_dist <- vegdist(community_matrix, method = "bray")

site_dispersion <- betadisper(bray_dist, community_metadata$Site)

anova(site_dispersion)

permutest(site_dispersion)

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
      "Site 1" = "#7B3294",   # Purple
      "Site 2" = "#2E8B57"    # Green
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
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold"),
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

##Bray curtis for plant species at the site level
#Do they cluster together in the nmds space

########################################################################
# NMDS OF PLANT SPECIES BY SITE
# Question:
# Do plant species cluster according to their pollinator communities,
# and do these patterns differ between Site 1 and Site 2?
########################################################################

# ============================================================
# LOAD PACKAGES
# ============================================================

library(dplyr)
library(tidyr)
library(tibble)
library(vegan)
library(ggplot2)


# ============================================================
# PREPARE POLLINATOR DATA
# ============================================================

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
    Plant.ID = factor(Plant.ID),
    
    # Use species-level identity for bumblebees and butterflies,
    # and broader pollinator groups for all other taxa
    Pollinator.Taxon = case_when(
      Group %in% c("Bumblebee", "Butterfly") &
        !is.na(Species) &
        Species != "" ~ Species,
      
      TRUE ~ Group
    )
  )


# ============================================================
# CHECK TAXA AND PLANT REPRESENTATION
# ============================================================

sort(unique(plant_community_data$Pollinator.Taxon))

table(
  plant_community_data$Site,
  plant_community_data$Plant.ID
)


# ============================================================
# CREATE PLANT × POLLINATOR COMMUNITY MATRIX
# ============================================================
# Each row represents one plant species/Plant.ID at one site.
# Each column represents one pollinator taxon.
# Cell values are total numbers of visits.

plant_matrix_data <- plant_community_data %>%
  group_by(
    Site,
    Plant.ID,
    Pollinator.Taxon
  ) %>%
  summarise(
    Visits = sum(Visits, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Pollinator.Taxon,
    values_from = Visits,
    values_fill = 0
  )


# ============================================================
# CREATE METADATA
# ============================================================

plant_metadata <- plant_matrix_data %>%
  select(
    Site,
    Plant.ID
  )


# ============================================================
# CREATE NUMERIC COMMUNITY MATRIX
# ============================================================

plant_matrix <- plant_matrix_data %>%
  select(
    -Site,
    -Plant.ID
  ) %>%
  as.data.frame()


# Add informative row names

rownames(plant_matrix) <- paste(
  plant_matrix_data$Site,
  plant_matrix_data$Plant.ID,
  sep = "_Plant_"
)


# ============================================================
# CHECK COMMUNITY MATRIX
# ============================================================

dim(plant_matrix)

head(plant_matrix)

rowSums(plant_matrix)

sum(rowSums(plant_matrix) == 0)

# No rows should have a total abundance of zero
if (any(rowSums(plant_matrix) == 0)) {
  stop(
    "At least one plant-site row contains no pollinator visits."
  )
}


# ============================================================
# RUN NMDS
# ============================================================

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


# ============================================================
# CHECK STRESS
# ============================================================

plant_nmds$stress

stress_label <- paste0(
  "Stress = ",
  round(plant_nmds$stress, 3)
)


# ============================================================
# EXTRACT NMDS SCORES
# ============================================================

plant_nmds_scores <- as.data.frame(
  scores(
    plant_nmds,
    display = "sites"
  )
)

plant_nmds_scores$Site <- plant_metadata$Site
plant_nmds_scores$Plant.ID <- plant_metadata$Plant.ID


# ============================================================
# OPTIONAL: ADD PLANT SPECIES NAMES
# ============================================================

plant_names <- c(
  "1"  = "Hypochaeris radicata",
  "2"  = "Dianthus deltoides",
  "3"  = "Thymus polytrichus",
  "4"  = "Lotus corniculatus",
  "5"  = "Silene viscaria",
  "6"  = "Galium saxatile",
  "7"  = "Hypochaeris radicata",
  "8"  = "Thymus polytrichus",
  "9"  = "Lotus corniculatus",
  "10" = "Silene viscaria",
  "11" = "Teucrium scorodonia",
  "12" = "Galium verum",
  "13" = "Galium verum",
  "14" = "Cirsium arvense"
)

plant_nmds_scores <- plant_nmds_scores %>%
  mutate(
    Plant.Species = unname(
      plant_names[as.character(Plant.ID)]
    )
  )


# ============================================================
# CALCULATE SITE CENTROIDS
# ============================================================

plant_site_centroids <- plant_nmds_scores %>%
  group_by(Site) %>%
  summarise(
    NMDS1 = mean(NMDS1),
    NMDS2 = mean(NMDS2),
    .groups = "drop"
  )


# ============================================================
# PLOT 1: COLOUR BY SITE, SHAPE BY PLANT ID
# ============================================================

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
    label = stress_label,
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


# ============================================================
# PLOT 2: LABEL POINTS WITH PLANT SPECIES
# ============================================================

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
    label = stress_label,
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


# ============================================================
# TEST WHETHER COMMUNITY COMPOSITION DIFFERS BETWEEN SITES
# ============================================================

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


# ============================================================
# CHECK HOMOGENEITY OF DISPERSION BETWEEN SITES
# ============================================================

plant_site_dispersion <- betadisper(
  plant_bray_dist,
  plant_metadata$Site
)

anova(plant_site_dispersion)

permutest(
  plant_site_dispersion,
  permutations = 999
)


# ============================================================
# OPTIONAL: TEST BOTH SITE AND PLANT ID
# ============================================================
# Use this cautiously because many Plant.ID levels occur at only one site.

plant_site_plant_permanova <- adonis2(
  plant_matrix ~ Site + Plant.ID,
  data = plant_metadata,
  method = "bray",
  permutations = 999
)

plant_site_plant_permanova


# ============================================================
# SAVE FIGURES
# ============================================================

ggsave(
  "Plant_pollinator_community_NMDS.png",
  plot = plant_nmds_plot,
  width = 9,
  height = 6,
  dpi = 300
)

ggsave(
  "Plant_pollinator_community_NMDS_labelled.png",
  plot = plant_nmds_labelled_plot,
  width = 10,
  height = 7,
  dpi = 300
)

##############################################################################
# Which co-flowering plant species share pollinators with restored species?
# Analysed separately by site
##############################################################################
##############################################################################
# CROSS-SITE HEATMAP OF SHARED POLLINATOR TAXA
#
# X-axis = plant species at Site 1
# Y-axis = plant species at Site 2
#
# Cell values = number of pollinator taxa shared between each pair of plants
#
# Pollinator taxon:
# - use Species where a valid species name is available
# - otherwise use the broader Group
##############################################################################

# Install once if needed:
# install.packages("ggtext")

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggtext)


# ============================================================
# 1. CHECK REQUIRED COLUMNS
# ============================================================

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


# ============================================================
# 2. PLANT NAMES
# ============================================================

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


# ============================================================
# 3. PREPARE DATA
# ============================================================

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


# ============================================================
# 4. CREATE PLANT × POLLINATOR PRESENCE MATRIX
# ============================================================

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


# ============================================================
# 5. MATCH POLLINATOR COLUMNS BETWEEN SITES
# ============================================================

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


# ============================================================
# 6. CALCULATE SHARED POLLINATOR TAXA ACROSS SITES
# ============================================================
#
# Rows = Site 2 plants
# Columns = Site 1 plants
#
# Matrix multiplication counts how many pollinator taxa
# are present on both plants.
# ============================================================

cross_site_shared_matrix <-
  site2_presence %*% t(site1_presence)


# ============================================================
# 7. CONVERT TO LONG FORMAT
# ============================================================

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


# ============================================================
# 8. CALCULATE POLLINATOR RICHNESS FOR EACH PLANT
# ============================================================

site1_richness <- rowSums(
  site1_presence > 0
)

site2_richness <- rowSums(
  site2_presence > 0
)


# ============================================================
# 9. ORDER PLANTS
# Restored plants first, followed by native plants alphabetically
# ============================================================

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


# ============================================================
# 10. AXIS LABELS
# Number in brackets = pollinator richness for that plant
# ============================================================

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


# ============================================================
# 11. CREATE HEATMAP
# ============================================================

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

# ============================================================
# 12. DISPLAY AND SAVE
# ============================================================

print(
  cross_site_shared_heatmap
)

ggsave(
  filename = "Cross_site_shared_pollinator_heatmap.png",
  plot = cross_site_shared_heatmap,
  width = 11,
  height = 8,
  dpi = 300
)

##############################################################################
# QUESTION 3
# Which co-flowering species have the greatest potential to compete with
# or facilitate pollination of restored species?
##############################################################################
##############################################################################
# QUESTION 3
#
# Which co-flowering species have the greatest potential to compete with
# or facilitate pollination of restored species?
#
# Pollinator sharing is defined as:
# the same pollinator taxon visiting a restored and native plant species
# at the same site and during the same time period.
##############################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)


# ============================================================
# 1. RESTORED PLANT SPECIES
# ============================================================

restored_species <- c(
  "Dianthus deltoides",
  "Silene viscaria"
)


# ============================================================
# 2. CHECK REQUIRED OBJECTS AND COLUMNS
# ============================================================

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


# ============================================================
# 3. CLEAN AND STANDARDISE DATA
# ============================================================

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


# ============================================================
# 4. IDENTIFY ALL RESTORED AND NATIVE PLANTS
#    PRESENT WITHIN EACH SITE × TIME PERIOD
# ============================================================

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


# ============================================================
# 5. CALCULATE CO-FLOWERING / TEMPORAL OVERLAP
#
# A restored and native species overlap when both were present
# at the same site during the same time period.
# ============================================================

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


# ============================================================
# 6. SUMMARISE RESTORED-PLANT INTERACTIONS
#    WITHIN EACH SITE × TIME PERIOD
# ============================================================

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


# ============================================================
# 7. SUMMARISE NATIVE-PLANT INTERACTIONS
#    WITHIN EACH SITE × TIME PERIOD
# ============================================================

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


# ============================================================
# 8. IDENTIFY SHARED POLLINATOR USE
#
# The join requires:
# - the same site;
# - the same time period;
# - the same pollinator taxon.
# ============================================================

period_sharing <- restored_period_interactions %>%
  inner_join(
    native_period_interactions,
    by = c(
      "Site",
      "Time.Period",
      "Pollinator.Taxon"
    )
  )


# ============================================================
# 9. SUMMARISE SHARING FOR EACH RESTORED–NATIVE PAIR
# ============================================================

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


# ============================================================
# 10. CREATE COMPLETE potential_interactions TABLE
#
# This retains all co-flowering restored–native pairs,
# including those with zero detected pollinator sharing.
# ============================================================

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


# ============================================================
# 11. SCREEN ALL NATIVE SPECIES
# ============================================================

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


# ============================================================
# 12. PREPARE BUBBLE-PLOT DATA
#
# The analysis screened all native species.
# The figure displays only pairs with detected sharing.
# ============================================================

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


# Check exactly which rows will appear as bubbles

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


# ============================================================
# 13. SET X-AXIS RANGE
# ============================================================

x_max <- max(
  interaction_plot_data$Shared.Pollinator.Taxa,
  na.rm = TRUE
)

x_upper <- ceiling(
  x_max
) + 0.5


# ============================================================
# 14. CREATE GREY BACKGROUND FOR DIANTHUS × SITE 2
# ============================================================

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


# ============================================================
# 15. CREATE BUBBLE PLOT
# ============================================================

interaction_plot <- ggplot(
  interaction_plot_data,
  aes(
    x = Shared.Pollinator.Taxa,
    y = Native.Species,
    size = Time.Period.Sharing.Events
  )
) +
  
  # Shade the missing Dianthus × Site 2 panel
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
  
  # Label the missing panel
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
  
  # Observed restored–native sharing pairs
  geom_point(
    shape = 21,
    fill = "plum1",
    colour = "grey20",
    stroke = 0.7,
    alpha = 0.9
  ) +
  
  # Restored species strips on the left
  facet_grid(
    Restored.Species ~ Site,
    scales = "free_y",
    space = "free_y",
    switch = "y",
    drop = FALSE
  ) +
  
  # Native species labels on the right
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


# ============================================================
# 16. DISPLAY AND SAVE
# ============================================================

print(
  interaction_plot
)


ggsave(
  filename = "Restored_native_shared_pollinator_bubble_plot.png",
  plot = interaction_plot,
  width = 11,
  height = 8,
  dpi = 300
)

# NETWORK-LEVEL METRICS
# Overall, site-level and Site × Time Period networks
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(bipartite)

# ------------------------------------------------------------
# 1. Prepare genuine plant-pollinator interactions
# ------------------------------------------------------------

network_data <- pollinators %>%
  filter(
    Group != "None",
    Visits > 0,
    !is.na(Family),
    trimws(Family) != "",
    !is.na(Plant.ID),
    !is.na(Site),
    !is.na(Time.Period)
  ) %>%
  mutate(
    Family = trimws(Family),
    Plant.ID = factor(Plant.ID),
    Site = factor(Site),
    Time.Period = factor(Time.Period)
  )


# ------------------------------------------------------------
# 2. Function to construct a quantitative interaction matrix
# ------------------------------------------------------------

make_network <- function(data) {
  
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
    Visits ~ Plant.ID + Family,
    data = data
  )
  
  # Convert the table into a plain numeric matrix
  web <- matrix(
    data = as.numeric(web_table),
    nrow = nrow(web_table),
    ncol = ncol(web_table),
    dimnames = dimnames(web_table)
  )
  
  # Remove plants and pollinator families with no interactions
  web <- web[
    rowSums(web) > 0,
    colSums(web) > 0,
    drop = FALSE
  ]
  
  web
}


# ------------------------------------------------------------
# 3. Function to extract modularity safely
# ------------------------------------------------------------

calculate_modularity <- function(web, reps = 100) {
  
  set.seed(123)
  
  module_result <- tryCatch(
    bipartite::DIRT_LPA_wb_plus(
      web,
      reps = reps
    ),
    error = function(e) NULL
  )
  
  if (is.null(module_result)) {
    return(NA_real_)
  }
  
  if ("modularity" %in% names(module_result)) {
    return(
      as.numeric(module_result[["modularity"]])
    )
  }
  
  if ("Q" %in% names(module_result)) {
    return(
      as.numeric(module_result[["Q"]])
    )
  }
  
  # Fallback in case the returned object has a different structure
  numeric_values <- suppressWarnings(
    as.numeric(
      unlist(module_result)
    )
  )
  
  numeric_values <- numeric_values[
    is.finite(numeric_values)
  ]
  
  if (length(numeric_values) == 0) {
    return(NA_real_)
  }
  
  tail(
    numeric_values,
    1
  )
}


# ------------------------------------------------------------
# 4. Function to calculate network-level metrics
# ------------------------------------------------------------

calculate_metrics <- function(web) {
  
  # Deal with completely empty matrices
  if (length(web) == 0) {
    
    return(
      tibble(
        Plant_Nodes = 0,
        Pollinator_Nodes = 0,
        Total_Visits = 0,
        Number_of_Links = 0,
        Connectance = NA_real_,
        Weighted_NODF = NA_real_,
        Modularity_Q = NA_real_
      )
    )
  }
  
  # Convert to a plain numeric matrix
  clean_web <- matrix(
    data = as.numeric(web),
    nrow = nrow(web),
    ncol = ncol(web),
    dimnames = dimnames(web)
  )
  
  # Remove empty plant and pollinator nodes
  clean_web <- clean_web[
    rowSums(clean_web) > 0,
    colSums(clean_web) > 0,
    drop = FALSE
  ]
  
  plant_nodes <- nrow(clean_web)
  pollinator_nodes <- ncol(clean_web)
  
  # Small networks cannot provide meaningful estimates
  # of nestedness or modularity
  if (
    plant_nodes < 2 ||
    pollinator_nodes < 2
  ) {
    
    return(
      tibble(
        Plant_Nodes = plant_nodes,
        Pollinator_Nodes = pollinator_nodes,
        Total_Visits = sum(clean_web),
        Number_of_Links = sum(clean_web > 0),
        Connectance = NA_real_,
        Weighted_NODF = NA_real_,
        Modularity_Q = NA_real_
      )
    )
  }
  
  # Calculate connectance
  connectance_value <- bipartite::networklevel(
    clean_web,
    index = "connectance"
  )
  
  # Calculate weighted nestedness
  nestedness_value <- bipartite::networklevel(
    clean_web,
    index = "weighted NODF"
  )
  
  # Calculate modularity
  modularity_value <- calculate_modularity(
    clean_web,
    reps = 100
  )
  
  # Return all metrics
  tibble(
    Plant_Nodes = plant_nodes,
    Pollinator_Nodes = pollinator_nodes,
    Total_Visits = sum(clean_web),
    Number_of_Links = sum(clean_web > 0),
    Connectance = as.numeric(connectance_value),
    Weighted_NODF = as.numeric(nestedness_value),
    Modularity_Q = as.numeric(modularity_value)
  )
}

# ------------------------------------------------------------
# 5. Create overall, site-level and temporal networks
# ------------------------------------------------------------

networks <- list(
  
  "Overall" = make_network(
    network_data
  ),
  
  "Site 1" = make_network(
    network_data %>%
      filter(Site == 1)
  ),
  
  "Site 2" = make_network(
    network_data %>%
      filter(Site == 2)
  ),
  
  "Site 1 - Time Period 1" = make_network(
    network_data %>%
      filter(
        Site == 1,
        Time.Period == 1
      )
  ),
  
  "Site 1 - Time Period 2" = make_network(
    network_data %>%
      filter(
        Site == 1,
        Time.Period == 2
      )
  ),
  
  "Site 1 - Time Period 3" = make_network(
    network_data %>%
      filter(
        Site == 1,
        Time.Period == 3
      )
  ),
  
  "Site 2 - Time Period 1" = make_network(
    network_data %>%
      filter(
        Site == 2,
        Time.Period == 1
      )
  ),
  
  "Site 2 - Time Period 2" = make_network(
    network_data %>%
      filter(
        Site == 2,
        Time.Period == 2
      )
  ),
  
  "Site 2 - Time Period 3" = make_network(
    network_data %>%
      filter(
        Site == 2,
        Time.Period == 3
      )
  )
)
# ------------------------------------------------------------
# 6. Calculate metrics for every network
# ------------------------------------------------------------

network_metrics <- bind_rows(
  lapply(
    names(networks),
    function(network_name) {
      
      calculate_metrics(
        networks[[network_name]]
      ) %>%
        mutate(
          Network = network_name,
          .before = 1
        )
    }
  )
) %>%
  mutate(
    Connectance = round(Connectance, 3),
    Weighted_NODF = round(Weighted_NODF, 2),
    Modularity_Q = round(Modularity_Q, 3)
  )

site_metrics <- network_metrics %>%
  dplyr::filter(
    Network %in% c(
      "Site 1",
      "Site 2"
    )
  ) %>%
  dplyr::select(
    Network,
    Plant_Nodes,
    Pollinator_Nodes,
    Total_Visits,
    Number_of_Links,
    Weighted_NODF,
    Modularity_Q
  )

site_metrics

# ============================================================
# SPECIES-LEVEL MODULE PLOTS BY SITE
# Scientific plant and pollinator names
# ============================================================

library(dplyr)
library(bipartite)

# ------------------------------------------------------------
# 1. Scientific plant-name lookup
# ------------------------------------------------------------

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
    "Cirsium arvense"       # 14
  ),
  
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 2. Prepare interaction data using the same taxonomic
#    resolution as the NMDS
# ------------------------------------------------------------

species_network_data <- network_data %>%
  left_join(
    plant_names,
    by = "Plant.ID"
  ) %>%
  filter(
    Group != "None",
    !is.na(Plant.Species)
  ) %>%
  mutate(
    Species = trimws(Species),
    
    # Use species-level identity wherever available.
    # Otherwise, retain the broader pollinator group.
    Pollinator.Label = case_when(
      !is.na(Species) &
        trimws(Species) != "" &
        trimws(Species) != "NA" &
        trimws(Species) != "Unknown" ~
        
        sub(
          "^([A-Za-z])[A-Za-z]+ ",
          "\\1. ",
          trimws(Species)
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


# ------------------------------------------------------------
# 3. Function to create a plant × pollinator matrix
# ------------------------------------------------------------

make_species_network <- function(data) {
  
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

# ------------------------------------------------------------
# 4. Create separate Site 1 and Site 2 matrices
# ------------------------------------------------------------

site1_species_web <- species_network_data %>%
  filter(Site == "Site 1") %>%
  make_species_network()

site2_species_web <- species_network_data %>%
  filter(Site == "Site 2") %>%
  make_species_network()


# Inspect the matrices
site1_species_web
site2_species_web

# ------------------------------------------------------------
# 5. Detect modules
# ------------------------------------------------------------

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

# ------------------------------------------------------------
# 6. Convert matrices into module webs
# ------------------------------------------------------------

site1_module_web <- bipartite::convert2moduleWeb(
  site1_species_web,
  site1_modules
)

site2_module_web <- bipartite::convert2moduleWeb(
  site2_species_web,
  site2_modules
)
# ------------------------------------------------------------
# Site 1
# ------------------------------------------------------------

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
  outer = TRUE,
  line = 1,
  font = 2,
  cex = 1.5
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
# ------------------------------------------------------------
# 8. Site 2 module plot
# ------------------------------------------------------------

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

# Calculate modules
site1_modules <- bipartite::computeModules(site1_species_web)

site2_modules <- bipartite::computeModules(site2_species_web)

site1_modules@likelihood

site2_modules@likelihood

# ============================================================
# WEIGHTED NESTEDNESS BY SITE
# ============================================================

library(dplyr)
library(ggplot2)
library(ggtext)
library(bipartite)

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

site1_species_web <- make_species_network(
  species_network_data %>%
    filter(Site == "Site 1")
)

site2_species_web <- make_species_network(
  species_network_data %>%
    filter(Site == "Site 2")
)

# ============================================================
# PREPARE NESTEDNESS HEATMAP DATA
# ============================================================

prepare_nestedness_data <- function(web) {
  
  # Sort plants and pollinators by decreasing interaction totals
  sorted_web <- bipartite::sortweb(
    web,
    sort.order = "dec"
  )
  
  # Remove empty rows and columns
  sorted_web <- sorted_web[
    rowSums(sorted_web) > 0,
    colSums(sorted_web) > 0,
    drop = FALSE
  ]
  
  # Convert matrix to long format
  nested_df <- as.data.frame(
    as.table(sorted_web)
  )
  
  names(nested_df) <- c(
    "Plant",
    "Pollinator",
    "Visits"
  )
  
  # Preserve the order returned by sortweb()
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


# Actually create the plotting data
site1_nested_df <- prepare_nestedness_data(
  site1_species_web
)

site2_nested_df <- prepare_nestedness_data(
  site2_species_web
)


# ============================================================
# HIGHLIGHT RESTORED PLANTS
# ============================================================

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


# ============================================================
# SITE 1 PLOT
# ============================================================

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


# ============================================================
# SITE 2 PLOT
# ============================================================

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

species_network_data %>%
  filter(Group == "Hoverfly") %>%
  distinct(
    Genus,
    Species,
    Pollinator.Label
  )

# ============================================================
# NULL-MODEL STANDARDISATION OF NETWORK METRICS
# Species-level weighted networks by site
#
# Metrics:
#   1. Weighted NODF
#   2. Modularity Q
#   3. Connectance
#
# Null models:
#   - R2d: r2dtable
#   - VazNull: vaznull
#
# Replicates:
#   - 500 null networks
#   - metaComputeModules N = 500
# ============================================================


# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

library(dplyr)
library(tibble)
library(bipartite)


# ------------------------------------------------------------
# 2. Check that the site matrices already exist
# ------------------------------------------------------------

if (!exists("site1_species_web")) {
  stop("site1_species_web does not exist. Run the species-network preparation code first.")
}

if (!exists("site2_species_web")) {
  stop("site2_species_web does not exist. Run the species-network preparation code first.")
}


# ------------------------------------------------------------
# 3. Clean a network matrix
# ------------------------------------------------------------

clean_network_matrix <- function(web) {
  
  web <- as.matrix(web)
  storage.mode(web) <- "numeric"
  
  web <- web[
    rowSums(web) > 0,
    colSums(web) > 0,
    drop = FALSE
  ]
  
  if (nrow(web) < 2 || ncol(web) < 2) {
    stop("The network must contain at least two non-empty rows and columns.")
  }
  
  web
}


# ------------------------------------------------------------
# 4. Weighted NODF
# ------------------------------------------------------------

calculate_weighted_nodf <- function(web) {
  
  web <- clean_network_matrix(web)
  
  as.numeric(
    bipartite::networklevel(
      web,
      index = "weighted NODF"
    )
  )
}


# ------------------------------------------------------------
# 5. Connectance
# ------------------------------------------------------------

calculate_connectance <- function(web) {
  
  web <- clean_network_matrix(web)
  
  as.numeric(
    bipartite::networklevel(
      web,
      index = "connectance"
    )
  )
}


# ------------------------------------------------------------
# 6. Modularity using metaComputeModules
# ------------------------------------------------------------

calculate_modularity_meta <- function(
    web,
    module_reps = 500
) {
  
  web <- clean_network_matrix(web)
  
  result <- bipartite::metaComputeModules(
    web,
    N = module_reps,
    method = "Beckett"
  )
  
  as.numeric(result@likelihood)
}


# ------------------------------------------------------------
# 7. Helper functions for null-model statistics
# ------------------------------------------------------------

calculate_ses <- function(
    observed,
    null_values
) {
  
  null_values <- null_values[
    is.finite(null_values)
  ]
  
  null_sd <- sd(null_values)
  
  if (
    length(null_values) < 2 ||
    !is.finite(null_sd) ||
    null_sd == 0
  ) {
    return(NA_real_)
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
  
  if (length(null_values) == 0) {
    return(NA_real_)
  }
  
  lower_p <- (
    sum(null_values <= observed) + 1
  ) / (
    length(null_values) + 1
  )
  
  upper_p <- (
    sum(null_values >= observed) + 1
  ) / (
    length(null_values) + 1
  )
  
  min(
    1,
    2 * min(lower_p, upper_p)
  )
}


# ------------------------------------------------------------
# 8. Weighted NODF null-model standardisation
# ------------------------------------------------------------

standardise_weighted_nodf <- function(
    web,
    null_reps = 500,
    seed = 123
) {
  
  web <- clean_network_matrix(web)
  
  observed <- calculate_weighted_nodf(web)
  
  set.seed(seed)
  
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
  
  set.seed(seed)
  
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
    
    Null_Mean_R2d = mean(r2d_values, na.rm = TRUE),
    Null_SD_R2d = sd(r2d_values, na.rm = TRUE),
    z_R2d = calculate_ses(
      observed,
      r2d_values
    ),
    p_R2d = calculate_two_tailed_p(
      observed,
      r2d_values
    ),
    
    Null_Mean_VazNull = mean(vaz_values, na.rm = TRUE),
    Null_SD_VazNull = sd(vaz_values, na.rm = TRUE),
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


# ------------------------------------------------------------
# 9. Modularity null-model standardisation
# ------------------------------------------------------------

standardise_modularity <- function(
    web,
    null_reps = 500,
    module_reps = 500,
    seed = 123
) {
  
  web <- clean_network_matrix(web)
  
  observed <- calculate_modularity_meta(
    web,
    module_reps = module_reps
  )
  
  set.seed(seed)
  
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
  
  set.seed(seed)
  
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
    
    Null_Mean_R2d = mean(r2d_values, na.rm = TRUE),
    Null_SD_R2d = sd(r2d_values, na.rm = TRUE),
    z_R2d = calculate_ses(
      observed,
      r2d_values
    ),
    p_R2d = calculate_two_tailed_p(
      observed,
      r2d_values
    ),
    
    Null_Mean_VazNull = mean(vaz_values, na.rm = TRUE),
    Null_SD_VazNull = sd(vaz_values, na.rm = TRUE),
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


# ------------------------------------------------------------
# 10. Connectance null-model standardisation
#
# VazNull preserves connectance, so only R2d is used for SES.
# ------------------------------------------------------------

standardise_connectance <- function(
    web,
    null_reps = 500,
    seed = 123
) {
  
  web <- clean_network_matrix(web)
  
  observed <- calculate_connectance(web)
  
  set.seed(seed)
  
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
    
    Null_Mean_R2d = mean(r2d_values, na.rm = TRUE),
    Null_SD_R2d = sd(r2d_values, na.rm = TRUE),
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


# ============================================================
# 11. TEST RUN
#
# Run this section first.
# It uses only 5 null networks and 5 modularity repeats.
# ============================================================

# ------------------------------------------------------------
# Confirm that the two matrices exist
# ------------------------------------------------------------

if (!exists("site1_species_web")) {
  stop(
    "site1_species_web does not exist. ",
    "Run the species-network preparation code first."
  )
}

if (!exists("site2_species_web")) {
  stop(
    "site2_species_web does not exist. ",
    "Run the species-network preparation code first."
  )
}


# ------------------------------------------------------------
# Clean local copies of the matrices for inspection
# ------------------------------------------------------------

site1_test_web <- clean_network_matrix(site1_species_web)
site2_test_web <- clean_network_matrix(site2_species_web)


# ------------------------------------------------------------
# Check that the site matrices are genuinely different
# ------------------------------------------------------------

cat("\nAre the Site 1 and Site 2 matrices identical?\n")

print(
  identical(
    site1_test_web,
    site2_test_web
  )
)


cat("\nMatrix dimensions:\n")

print(
  data.frame(
    Site = c("Site 1", "Site 2"),
    Plant_Nodes = c(
      nrow(site1_test_web),
      nrow(site2_test_web)
    ),
    Pollinator_Nodes = c(
      ncol(site1_test_web),
      ncol(site2_test_web)
    ),
    Total_Visits = c(
      sum(site1_test_web),
      sum(site2_test_web)
    ),
    Number_of_Links = c(
      sum(site1_test_web > 0),
      sum(site2_test_web > 0)
    )
  )
)


# ------------------------------------------------------------
# Site 1 test analyses
# ------------------------------------------------------------

site1_nodf_test <- standardise_weighted_nodf(
  site1_species_web,
  null_reps = 5,
  seed = 123
)

site1_modularity_test <- standardise_modularity(
  site1_species_web,
  null_reps = 5,
  module_reps = 5,
  seed = 123
)

site1_connectance_test <- standardise_connectance(
  site1_species_web,
  null_reps = 5,
  seed = 123
)


# ------------------------------------------------------------
# Site 2 test analyses
# ------------------------------------------------------------

site2_nodf_test <- standardise_weighted_nodf(
  site2_species_web,
  null_reps = 5,
  seed = 456
)

site2_modularity_test <- standardise_modularity(
  site2_species_web,
  null_reps = 5,
  module_reps = 5,
  seed = 456
)

site2_connectance_test <- standardise_connectance(
  site2_species_web,
  null_reps = 5,
  seed = 456
)


# ------------------------------------------------------------
# Combine the test results
# ------------------------------------------------------------

site1_test_results <- dplyr::bind_rows(
  site1_nodf_test$summary,
  site1_modularity_test$summary,
  site1_connectance_test$summary
) %>%
  dplyr::mutate(
    Site = "Site 1",
    .before = 1
  )

site2_test_results <- dplyr::bind_rows(
  site2_nodf_test$summary,
  site2_modularity_test$summary,
  site2_connectance_test$summary
) %>%
  dplyr::mutate(
    Site = "Site 2",
    .before = 1
  )

test_results <- dplyr::bind_rows(
  site1_test_results,
  site2_test_results
)


# ------------------------------------------------------------
# Print all test results
# ------------------------------------------------------------

cat("\nTest-run null-model results:\n")

print(
  test_results,
  width = Inf
)


# ------------------------------------------------------------
# Directly calculate the three observed metrics
# ------------------------------------------------------------

observed_test_metrics <- tibble::tibble(
  Site = c("Site 1", "Site 2"),
  
  Weighted_NODF = c(
    calculate_weighted_nodf(site1_species_web),
    calculate_weighted_nodf(site2_species_web)
  ),
  
  Modularity_Q = c(
    calculate_modularity_meta(
      site1_species_web,
      module_reps = 5
    ),
    calculate_modularity_meta(
      site2_species_web,
      module_reps = 5
    )
  ),
  
  Connectance = c(
    calculate_connectance(site1_species_web),
    calculate_connectance(site2_species_web)
  )
)

cat("\nObserved site metrics:\n")

print(
  observed_test_metrics,
  width = Inf
)


# ------------------------------------------------------------
# Manual connectance check
#
# Connectance =
# observed links / possible links
# ------------------------------------------------------------

manual_connectance_check <- tibble::tibble(
  Site = c("Site 1", "Site 2"),
  
  Plant_Nodes = c(
    nrow(site1_test_web),
    nrow(site2_test_web)
  ),
  
  Pollinator_Nodes = c(
    ncol(site1_test_web),
    ncol(site2_test_web)
  ),
  
  Observed_Links = c(
    sum(site1_test_web > 0),
    sum(site2_test_web > 0)
  )
) %>%
  dplyr::mutate(
    Possible_Links = Plant_Nodes * Pollinator_Nodes,
    Manual_Connectance = Observed_Links / Possible_Links
  )

cat("\nManual connectance check:\n")

print(
  manual_connectance_check,
  width = Inf
)


# ------------------------------------------------------------
# Optional: inspect the matrices
# ------------------------------------------------------------

cat("\nSite 1 matrix:\n")
print(site1_test_web)

cat("\nSite 2 matrix:\n")
print(site2_test_web)

# ============================================================
# 12. FINAL ANALYSIS
#
# Warning:x
# Modularity with null_reps = 500 and module_reps = 500
# may take a very long time.
# ============================================================


# ------------------------------------------------------------
# Site 1
# ------------------------------------------------------------

site1_nodf_results <- standardise_weighted_nodf(
  site1_species_web,
  null_reps = 500
)

site1_modularity_results <- standardise_modularity(
  site1_species_web,
  null_reps = 500,
  module_reps = 500
)

site1_connectance_results <- standardise_connectance(
  site1_species_web,
  null_reps = 500
)


# ------------------------------------------------------------
# Site 2
# ------------------------------------------------------------

site2_nodf_results <- standardise_weighted_nodf(
  site2_species_web,
  null_reps = 500
)

site2_modularity_results <- standardise_modularity(
  site2_species_web,
  null_reps = 500,
  module_reps = 500
)

site2_connectance_results <- standardise_connectance(
  site2_species_web,
  null_reps = 500
)


# ============================================================
# 13. Combine the results
# ============================================================

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


# ============================================================
# 14. Rounded reporting table
# ============================================================

network_null_results_report <- network_null_results %>%
  mutate(
    Observed = round(Observed, 3),
    
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


# ============================================================
# 15. Inspect null distributions
# ============================================================

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

site1_modularity_results <- standardise_modularity(
  site1_species_web,
  null_reps = 500,
  module_reps = 5
)

site2_modularity_results <- standardise_modularity(
  site2_species_web,
  null_reps = 500,
  module_reps = 5
)

site1_modularity_results$summary

site2_modularity_results$summary

### MAIDEN PINK

# Pollinator taxa recorded visiting Maiden pink at Site 1
dianthus_taxa <- community_data %>%
  filter(
    Site == 1,
    Plant.ID == 2
  ) %>%
  pull(Pollinator.Taxon) %>%
  unique()

sort(dianthus_taxa)


# All pollinator taxa recorded anywhere at Site 2
site2_taxa <- community_data %>%
  filter(Site == 2) %>%
  pull(Pollinator.Taxon) %>%
  unique()

sort(site2_taxa)


# Maiden pink pollinators that were also present at Site 2
shared_taxa <- intersect(dianthus_taxa, site2_taxa)

shared_taxa
length(shared_taxa)


# Maiden pink pollinators that were absent from Site 2
missing_taxa <- setdiff(dianthus_taxa, site2_taxa)

missing_taxa
length(missing_taxa)


# Number of Maiden pink visits made by each missing taxon
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


# Total Maiden pink visits and percentage made by missing taxa
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

#"What percentage of each plant's visits came from each pollinator taxon?"

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

#How many pollinator taxa are needed to reach 80%?
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

#Shannon Diversity
library(vegan)

community_data %>%
  group_by(Plant.ID, Pollinator.Taxon) %>%
  summarise(Visits = sum(Visits), .groups="drop") %>%
  tidyr::pivot_wider(names_from = Pollinator.Taxon,
                     values_from = Visits,
                     values_fill = 0)

#Sticky Catchfly
catchfly <- community_data %>%
  filter(Plant.ID %in% c(5, 10))
catchfly_summary <- catchfly %>%
  group_by(Site, Plant.ID) %>%
  summarise(
    Total.Visits = sum(Visits),
    Mean.Flowers = mean(Number.of.Flowers),
    Pollinator.Richness = n_distinct(Pollinator.Taxon),
    .groups = "drop"
  ) %>%
  mutate(
    Visits.Per.Flower = Total.Visits / Mean.Flowers
  )

catchfly_summary
community_data %>%
  filter(Site == 2, Plant.ID == 5)
catchfly_summary <- catchfly %>%
  group_by(Site, Plant.ID) %>%
  summarise(
    Total.Visits = sum(Visits),
    Mean.Flowers = mean(Number.of.Flowers),
    Pollinator.Richness = n_distinct(Pollinator.Taxon),
    .groups = "drop"
  ) %>%
  mutate(
    Visits.Per.Flower = Total.Visits / Mean.Flowers
  )

catchfly_summary
community_data %>%
  filter(
    (Site == 1 & Plant.ID == 5) |
      (Site == 2 & Plant.ID == 10)
  ) %>%
  group_by(Site, Pollinator.Taxon) %>%
  summarise(Visits = sum(Visits), .groups = "drop") %>%
  arrange(Site, desc(Visits))


unique(pollinators$Group)

library(MASS)

#With and without time period

library(MASS)

pollinators$Temperature_centred <- scale(
  pollinators$Temperature..C.degrees.,
  center = TRUE,
  scale = FALSE
)[, 1]

pollinators$Solar_centred <- scale(
  pollinators$Solar.Power..W.m.2.,
  center = TRUE,
  scale = FALSE
)[, 1]

pollinators$Group <- factor(pollinators$Group)
pollinators$Plant.ID <- factor(pollinators$Plant.ID)
pollinators$Time.Period <- factor(pollinators$Time.Period)

abiotic_data <- pollinators[
  complete.cases(
    pollinators[, c(
      "Visits",
      "Temperature_centred",
      "Solar_centred",
      "Humidity....RH.",
      "Wind.Speed..m.s.",
      "Group",
      "Plant.ID",
      "Time.Period"
    )]
  ),
]

#Model without time period
library(MASS)

AbioticModel <- glm.nb(
  Visits ~
    Temperature_centred +
    I(Temperature_centred^2) +
    Solar_centred +
    I(Solar_centred^2) +
    Humidity....RH. +
    Wind.Speed..m.s. +
    Group +
    Plant.ID,
  data = abiotic_data
)

#Model with time period
AbioticTimeModel <- glm.nb(
  Visits ~
    Temperature_centred +
    I(Temperature_centred^2) +
    Solar_centred +
    I(Solar_centred^2) +
    Humidity....RH. +
    Wind.Speed..m.s. +
    Group +
    Plant.ID +
    Time.Period,
  data = abiotic_data
)

#Compare
anova(
  AbioticModel,
  AbioticTimeModel,
  test = "Chisq"
)

AIC(
  AbioticModel,
  AbioticTimeModel
)

library(emmeans)

emmeans(AbioticTimeModel, pairwise ~ Time.Period, type = "response")


library(ggplot2)

Time_plot <- as.data.frame(Time_emmeans$emmeans)

ggplot(
  Time_plot,
  aes(x = Time.Period, y = response)
) +
  geom_point(size = 4) +
  geom_errorbar(
    aes(
      ymin = asymp.LCL,
      ymax = asymp.UCL
    ),
    width = 0.15,
    linewidth = 0.8
  ) +
  labs(
    x = "Time period",
    y = "Predicted pollinator visits\nper 20-minute survey"
  ) +
  theme_classic(base_size = 16)

ggplot(Time_plot, aes(Time.Period, response)) +
  
  geom_errorbar(
    aes(ymin = asymp.LCL,
        ymax = asymp.UCL),
    width = 0.12,
    linewidth = 1.2,
    colour = "#6BAED6"
  ) +
  
geom_point(
  size = 5,
  colour = "#2171B5"
) +
  
  theme_classic(base_size = 18)
library(ggplot2)

Time_plot <- as.data.frame(Time_emmeans$emmeans)

library(ggplot2)

Time_plot <- as.data.frame(Time_emmeans$emmeans)

ggplot(
  Time_plot,
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
    width = 0.10,
    linewidth = 1.1
  ) +
  
  geom_point(
    size = 5.5
  ) +
  
  # Significance bracket between Time Periods 1 and 2
  annotate(
    "segment",
    x = 1,
    xend = 2,
    y = 26,
    yend = 26,
    linewidth = 0.8
  ) +
  
  annotate(
    "segment",
    x = 1,
    xend = 1,
    y = 25.5,
    yend = 26,
    linewidth = 0.8
  ) +
  
  annotate(
    "segment",
    x = 2,
    xend = 2,
    y = 25.5,
    yend = 26,
    linewidth = 0.8
  ) +
  
  annotate(
    "text",
    x = 1.5,
    y = 26.5,
    label = "p = 0.033",
    size = 5
  ) +
  
  scale_colour_manual(
    values = c(
      "1" = "#E3B23C",
      "2" = "#F28E2B",
      "3" = "#4E79A7"
    )
  ) +
  
  scale_x_discrete(
    labels = c(
      "1" = "Time Period 1",
      "2" = "Time Period 2",
      "3" = "Time Period 3"
    )
  ) +
  
  scale_y_continuous(
    limits = c(0, 28),
    breaks = seq(0, 25, 5),
    expand = expansion(mult = c(0, 0))
  ) +
  
  labs(
    x = NULL,
    y = "Predicted visits per\n20-minute survey"
  ) +
  
  theme_classic(
    base_size = 18
  ) +
  
  theme(
    legend.position = "none",
    axis.title.y = element_text(
      face = "bold",
      margin = margin(r = 12)
    ),
    axis.text = element_text(
      colour = "black"
    ),
    axis.line = element_line(
      linewidth = 0.8
    ),
    axis.ticks = element_line(
      linewidth = 0.8
    )
  )

sum(pollinators$Visits, na.rm = TRUE)

# Total visitation records
n_records <- nrow(pollinators)

# Total visits
total_visits <- sum(pollinators$Visits, na.rm = TRUE)

# Plant species
n_plant_species <- pollinators %>%
  filter(Group != "None") %>%
  distinct(Plant.ID) %>%
  nrow()

# Pollinator taxa (highest taxonomic resolution used in the networks)
n_pollinator_taxa <- species_network_data %>%
  distinct(Pollinator.Label) %>%
  nrow()

n_records
total_visits
n_plant_species
n_pollinator_taxa
head(n_pollinator_taxa)
head(pollinators)
species_network_data %>%
  distinct(Pollinator.Label) %>%
  arrange(Pollinator.Label)

library(dplyr)

n_surveys <- pollinators %>%
  distinct(Date, Site, Plant.ID, Time.Period, Replication) %>%
  nrow()

n_surveys

emmeans(
  ModelSpecies_nb,
  ~ Plant.Species,
  type = "response",
  offset = 0
)

packageVersion("bipartite")
help.search("nullmodel", package = "bipartite")


richness_data <- pollinators %>%
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
    TaxonRichness = n_distinct(Group[Visits > 0]),
    .groups = "drop"
  )

summary(pollinators$Humidity....RH.)
sd(pollinators$Humidity....RH., na.rm = TRUE)
range(pollinators$Humidity....RH., na.rm = TRUE)
