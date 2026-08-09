#######################################
#Pollinator Group Temperature Responses
#######################################

# Aim 5:
# Determine whether pollinator groups respond differently
# to temperature.
#
# Negative-binomial GLMM with survey Date as a random intercept.


source("Overall Survey and Data.R")


#Prep data

# Mining Bee and Wasp are excluded because there are too few
# observations to estimate reliable group-specific responses.

analysis_data_temp <- analysis_data %>%
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
      
      TRUE ~
        Date
    )
  ) %>%
  
  filter(
    !(Group %in% c(
      "Mining Bee",
      "Wasp"
    )),
    
    !is.na(Visits),
    !is.na(Group),
    
    !is.na(
      Temperature..C.degrees.
    ),
    
    !is.na(Site),
    
    !is.na(
      Time.Period
    ),
    
    !is.na(
      Plant.ID
    ),
    
    !is.na(Date)
  ) %>%
  
  mutate(
    Group = droplevels(
      factor(
        Group
      )
    ),
    
    Site = factor(
      Site
    ),
    
    Time.Period = factor(
      Time.Period
    ),
    
    Plant.ID = droplevels(
      factor(
        Plant.ID
      )
    ),
    
    Date = droplevels(
      factor(
        Date
      )
    )
  )


#Check group and date sample sizes

table(
  analysis_data_temp$Group
)

table(
  analysis_data_temp$Date
)


#Check temperature ranges by pollinator group

temperature_ranges <- analysis_data_temp %>%
  group_by(
    Group
  ) %>%
  
  summarise(
    Observations = n(),
    
    Distinct.Temperatures = n_distinct(
      Temperature..C.degrees.
    ),
    
    Minimum.Temperature = min(
      Temperature..C.degrees.,
      na.rm = TRUE
    ),
    
    Maximum.Temperature = max(
      Temperature..C.degrees.,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

temperature_ranges


#Mean-centre temperature

mean_temperature <- mean(
  analysis_data_temp$
    Temperature..C.degrees.,
  na.rm = TRUE
)


analysis_data_temp <- analysis_data_temp %>%
  mutate(
    Temperature_c =
      Temperature..C.degrees. -
      mean_temperature
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
  
  data =
    analysis_data_temp
)


summary(
  ModelTempGroup_linear
)


TempGroup_linear_Anova <- car::Anova(
  ModelTempGroup_linear,
  type = 2
)

TempGroup_linear_Anova


#Prep data for quadratic analysis

# Sweat bees are excluded from the quadratic analysis because
# their observations span too few distinct temperatures to
# estimate a group-specific quadratic curve reliably.

analysis_data_quad <- analysis_data_temp %>%
  filter(
    Group != "Sweat Bee"
  ) %>%
  droplevels()


table(
  analysis_data_quad$Group
)

table(
  analysis_data_quad$Date
)


#Re-centre temperature using the reduced dataset

mean_temperature_quad <- mean(
  analysis_data_quad$
    Temperature..C.degrees.,
  na.rm = TRUE
)


analysis_data_quad <- analysis_data_quad %>%
  mutate(
    Temperature_c =
      Temperature..C.degrees. -
      mean_temperature_quad
  )


#Linear model on the same quadratic-analysis dataset

ModelTempGroup_linear_quad_data <- glmmTMB(
  Visits ~
    Temperature_c * Group +
    Site +
    Time.Period +
    Plant.ID +
    (1 | Date),
  
  family = nbinom2,
  
  data =
    analysis_data_quad
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
  
  data =
    analysis_data_quad
)


summary(
  ModelTempGroup_quad
)


TempGroup_quad_Anova <- car::Anova(
  ModelTempGroup_quad,
  type = 2
)

TempGroup_quad_Anova

# Temperature x Group was significant.
# Temperature^2 x Group was not significant.


#Does a group-specific quadratic response improve fit?

anova(
  ModelTempGroup_linear_quad_data,
  ModelTempGroup_quad
)


AIC(
  ModelTempGroup_linear_quad_data,
  ModelTempGroup_quad
)

# The full group-specific quadratic model did not significantly
# improve fit.
#
# The AIC difference was <2, indicating similar support.


#Simplified Model

# Shared quadratic temperature response, but group-specific
# linear temperature slopes.

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
  
  data =
    analysis_data_quad
)


summary(
  ModelTempGroup_final
)


TempGroup_final_Anova <- car::Anova(
  ModelTempGroup_final,
  type = 2
)

TempGroup_final_Anova


#Compare simplified model with full quadratic interaction

anova(
  ModelTempGroup_final,
  ModelTempGroup_quad
)


AIC(
  ModelTempGroup_final,
  ModelTempGroup_quad
)

# The simplified model has comparable support while avoiding
# unsupported group-specific quadratic temperature terms.


#Test overall linear temperature x group interaction

# This reduced model removes group-specific linear temperature
# slopes while retaining the shared quadratic response.

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
  
  data =
    analysis_data_quad
)


anova(
  ModelTempGroup_no_linear_interaction,
  ModelTempGroup_final
)


AIC(
  ModelTempGroup_no_linear_interaction,
  ModelTempGroup_final
)

# A significant likelihood-ratio test indicates that the
# temperature response differs among pollinator groups.


#Check Date random effect

VarCorr(
  ModelTempGroup_final
)


ModelTempGroup_final$sdr$pdHess


#Compare with model without Date

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
  
  data =
    analysis_data_quad
)


AIC(
  ModelTempGroup_no_date,
  ModelTempGroup_final
)

# Retain Date if the random-effect variance, model stability
# and AIC support its inclusion.


#Model diagnostics

TempGroup_residuals <- simulateResiduals(
  fittedModel =
    ModelTempGroup_final,
  
  n = 1000
)


plot(
  TempGroup_residuals
)


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

# Only prefer nbinom1 if it converges cleanly,
# has lower AIC and has satisfactory diagnostics.


#Generate predictions within each group's observed temperature range

group_temperature_ranges <- analysis_data_quad %>%
  group_by(
    Group
  ) %>%
  
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
  
  group_by(
    Group
  ) %>%
  
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
      Temperature -
      mean_temperature_quad,
    
    Group = factor(
      Group,
      levels = levels(
        analysis_data_quad$Group
      )
    ),
    
    Site = factor(
      levels(
        analysis_data_quad$Site
      )[1],
      
      levels = levels(
        analysis_data_quad$Site
      )
    ),
    
    Time.Period = factor(
      levels(
        analysis_data_quad$
          Time.Period
      )[1],
      
      levels = levels(
        analysis_data_quad$
          Time.Period
      )
    ),
    
    Plant.ID = factor(
      levels(
        analysis_data_quad$
          Plant.ID
      )[1],
      
      levels = levels(
        analysis_data_quad$
          Plant.ID
      )
    ),
    
    # A valid Date is required in newdata.
    # re.form = NA excludes the date-specific random effect.
    
    Date = factor(
      levels(
        analysis_data_quad$Date
      )[1],
      
      levels = levels(
        analysis_data_quad$Date
      )
    )
  )


prediction_results <- predict(
  ModelTempGroup_final,
  
  newdata =
    Temp_predictions_final_df,
  
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
        1.96 *
        prediction_results$se.fit
    ),
    
    upper_CL = exp(
      prediction_results$fit +
        1.96 *
        prediction_results$se.fit
    )
  )


head(
  Temp_predictions_final_df
)


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
      x =
        Temperature..C.degrees.,
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
    values =
      pollinator_colours
  ) +
  
  scale_fill_manual(
    values =
      pollinator_colours
  ) +
  
  labs(
    x = "Temperature (°C)",
    y = "Predicted number of visits",
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


#Final Aim 5 outputs

temperature_ranges

TempGroup_linear_Anova

TempGroup_quad_Anova

TempGroup_final_Anova

anova(
  ModelTempGroup_no_linear_interaction,
  ModelTempGroup_final
)

AIC(
  ModelTempGroup_no_linear_interaction,
  ModelTempGroup_final
)

VarCorr(
  ModelTempGroup_final
)

ModelTempGroup_final$sdr$pdHess

TempGroup_uniformity

TempGroup_dispersion

TempGroup_zero_inflation

Temperature_plot
