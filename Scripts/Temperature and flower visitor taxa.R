# TEMPERATURE AND FLOWER-VISITOR TAXA
# Do different flower-visitor taxa respond differently to temperature?


# Packages

library(dplyr)
library(tidyr)
library(ggplot2)
library(glmmTMB)
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


#Prep data

analysis_data_temp <- pollinators %>%
  mutate(
    Group = trimws(
      as.character(Group)
    ),
    
    Date = as.Date(
      Date,
      format = "%d/%m/%Y"
    ),
    
    #Fix dates
    
    Date = case_when(
      Date == as.Date("2026-06-06") ~
        as.Date("2026-06-05"),
      
      Date == as.Date("2026-07-03") ~
        as.Date("2026-07-02"),
      
      TRUE ~ Date
    )
  ) %>%
  
  filter(
    !(Group %in% c(
      "None",
      "Mining Bee",
      "Wasp"
    )),
    
    !is.na(Visits),
    !is.na(Group),
    !is.na(Temperature..C.degrees.),
    !is.na(Site),
    !is.na(Time.Period),
    !is.na(Plant.ID),
    !is.na(Date)
  ) %>%
  
  mutate(
    Group = droplevels(
      factor(Group)
    ),
    
    Site = factor(Site),
    
    Time.Period = factor(
      Time.Period
    ),
    
    Plant.ID = droplevels(
      factor(Plant.ID)
    ),
    
    Date = factor(Date)
  )


#Check group and date sample sizes

table(
  analysis_data_temp$Group
)

table(
  analysis_data_temp$Date
)


#Prepare survey x flower-visitor group dataset

#This analysis requires one row for each survey x flower-visitor group.
#Multiple records from the same flower-visitor group within the same
#20-minute survey are therefore summed together.

#Flower-visitor groups included in the final quadratic temperature analysis

groups_quad <- c(
  "Bumblebee",
  "Butterfly",
  "Honeybee",
  "Hoverfly",
  "Non-syrphid Fly"
)


#Create one row for each survey and retain the temperature
#recorded during that survey

survey_conditions <- analysis_data_temp %>%
  filter(
    Group %in% groups_quad
  ) %>%
  
  group_by(
    Date,
    Site,
    Plant.ID,
    Time.Period,
    Block,
    Replication
  ) %>%
  
  summarise(
    Temperature..C.degrees. =
      first(
        Temperature..C.degrees.
      ),
    
    .groups = "drop"
  )


#Sum all visits from the same flower-visitor group within
#each 20-minute survey

survey_group_visits <- analysis_data_temp %>%
  filter(
    Group %in% groups_quad
  ) %>%
  
  group_by(
    Date,
    Site,
    Plant.ID,
    Time.Period,
    Block,
    Replication,
    Group
  ) %>%
  
  summarise(
    Visits = sum(
      Visits,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


#Create every survey x flower-visitor group combination

#This adds explicit zeroes where a flower-visitor group was
#not recorded during a particular survey.

analysis_data_quad2 <- survey_conditions %>%
  tidyr::crossing(
    Group = groups_quad
  ) %>%
  
  left_join(
    survey_group_visits,
    by = c(
      "Date",
      "Site",
      "Plant.ID",
      "Time.Period",
      "Block",
      "Replication",
      "Group"
    )
  ) %>%
  
  mutate(
    Visits = replace_na(
      Visits,
      0
    ),
    
    Group = factor(
      Group,
      levels = groups_quad
    ),
    
    Date = factor(Date),
    Site = factor(Site),
    Plant.ID = factor(Plant.ID),
    Time.Period = factor(Time.Period)
  )


#Check final analysis structure

nrow(
  analysis_data_quad2
)

table(
  analysis_data_quad2$Group
)

table(
  analysis_data_quad2$Visits == 0
)


#Check that there is now only one row for each
#survey x flower-visitor group combination

analysis_data_quad2 %>%
  count(
    Date,
    Site,
    Plant.ID,
    Time.Period,
    Block,
    Replication,
    Group
  ) %>%
  
  count(
    n,
    name = "number_of_combinations"
  )


#Mean-centre temperature

mean_temperature_quad2 <- mean(
  analysis_data_quad2$Temperature..C.degrees.,
  na.rm = TRUE
)

analysis_data_quad2 <- analysis_data_quad2 %>%
  mutate(
    Temperature_c =
      Temperature..C.degrees. -
      mean_temperature_quad2
  )


#Check temperature range

summary(
  analysis_data_quad2$Temperature..C.degrees.
)


#Full group-specific quadratic model

#This model allows both the linear and quadratic effects of
#temperature to vary among flower-visitor groups.

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
  data = analysis_data_quad2
)

summary(
  ModelTempGroup_quad
)

TempGroup_quad_Anova <- car::Anova(
  ModelTempGroup_quad,
  type = 2
)

TempGroup_quad_Anova

#Temperature x Group was significant.
#Temperature^2 x Group was not significant.


#Simplified model

#Shared quadratic temperature response, but group-specific
#linear temperature slopes.

ModelTempGroup_with_date <- glmmTMB(
  Visits ~
    Temperature_c +
    I(Temperature_c^2) +
    Group +
    Temperature_c:Group +
    Time.Period +
    Plant.ID +
    (1 | Date),
  
  family = nbinom2,
  data = analysis_data_quad2
)


#Compare the simplified model with the full quadratic interaction

anova(
  ModelTempGroup_with_date,
  ModelTempGroup_quad
)

AIC(
  ModelTempGroup_with_date,
  ModelTempGroup_quad
)

#Full Temperature^2 x Group interaction did not significantly improve fit.
#Simplified model was retained.


#Check random date effect

VarCorr(
  ModelTempGroup_with_date
)

ModelTempGroup_with_date$sdr$pdHess


#Same model without Date random effect

ModelTempGroup_no_date <- glmmTMB(
  Visits ~
    Temperature_c +
    I(Temperature_c^2) +
    Group +
    Temperature_c:Group +
    Time.Period +
    Plant.ID,
  
  family = nbinom2,
  data = analysis_data_quad2
)


AIC(
  ModelTempGroup_with_date,
  ModelTempGroup_no_date
)


#Date was initially included as a random effect.
#Its estimated variance was effectively zero and removing
#Date reduced AIC from 1140.41 to 1138.41.
#The simpler model without Date was therefore retained.


#Final temperature x flower-visitor group model

#Site is excluded because Plant ID is nested within site.

ModelTempGroup_final <- glmmTMB(
  Visits ~
    Temperature_c +
    I(Temperature_c^2) +
    Group +
    Temperature_c:Group +
    Time.Period +
    Plant.ID,
  
  family = nbinom2,
  data = analysis_data_quad2
)


#Check that the correct model has been fitted

formula(
  ModelTempGroup_final
)

ModelTempGroup_final$sdr$pdHess

AIC(
  ModelTempGroup_final
)

summary(
  ModelTempGroup_final
)


#Overall fixed-effect tests

TempGroup_final_Anova <- car::Anova(
  ModelTempGroup_final,
  type = 2
)

TempGroup_final_Anova


#Test the overall linear temperature x group interaction

#This reduced model removes group-specific linear temperature slopes
#while retaining the shared quadratic temperature response.

ModelTempGroup_no_linear_interaction <- glmmTMB(
  Visits ~
    Temperature_c +
    I(Temperature_c^2) +
    Group +
    Time.Period +
    Plant.ID,
  
  family = nbinom2,
  data = analysis_data_quad2
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
#temperature effect differs among flower-visitor groups.


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

#nbinom2 had a much lower AIC than nbinom1.
#Both converged but nbinom2 was retained.


#Model diagnostics

set.seed(123)

TempGroup_residuals <- simulateResiduals(
  fittedModel = ModelTempGroup_final,
  n = 1000
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


#Four-panel diagnostic figure

par(
  mfrow = c(2, 2),
  mar = c(5, 5, 3, 2),
  cex = 0.8
)


#Panel 1: QQ residual diagnostic

plotQQunif(
  TempGroup_residuals
)


#Panel 2: Residuals against fitted values

plotResiduals(
  TempGroup_residuals,
  form = fitted(
    ModelTempGroup_final
  )
)


#Panel 3: Dispersion diagnostic

testDispersion(
  TempGroup_residuals,
  plot = TRUE
)


#Panel 4: Zero-inflation diagnostic

testZeroInflation(
  TempGroup_residuals,
  plot = TRUE
)


#Reset plotting layout afterwards

par(
  mfrow = c(1, 1)
)


#Group-specific temperature slopes from the FINAL model

temp_gradients_final <- emtrends(
  ModelTempGroup_final,
  ~ Group,
  var = "Temperature_c"
)

temp_gradient_results <- summary(
  temp_gradients_final,
  infer = c(
    TRUE,
    TRUE
  )
)

temp_gradient_results


#Generate predictions from the FINAL temperature x flower-visitor group model

#Calculate the observed temperature range for each flower-visitor group

group_temperature_ranges_final <- analysis_data_quad2 %>%
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


group_temperature_ranges_final


#Create a sequence of temperatures within the observed
#temperature range of each flower-visitor group

Temp_predictions_final_df <- group_temperature_ranges_final %>%
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
      mean_temperature_quad2,
    
    Group = factor(
      Group,
      levels = levels(
        analysis_data_quad2$Group
      )
    ),
    
    Time.Period = factor(
      levels(
        analysis_data_quad2$Time.Period
      )[1],
      levels = levels(
        analysis_data_quad2$Time.Period
      )
    ),
    
    Plant.ID = factor(
      levels(
        analysis_data_quad2$Plant.ID
      )[1],
      levels = levels(
        analysis_data_quad2$Plant.ID
      )
    )
  )


#Generate predicted visitation and 95% confidence intervals

prediction_results_final <- predict(
  ModelTempGroup_final,
  newdata = Temp_predictions_final_df,
  type = "link",
  se.fit = TRUE
)


Temp_predictions_final_df <- Temp_predictions_final_df %>%
  mutate(
    predicted_visits = exp(
      prediction_results_final$fit
    ),
    
    lower_CL = exp(
      prediction_results_final$fit -
        1.96 *
        prediction_results_final$se.fit
    ),
    
    upper_CL = exp(
      prediction_results_final$fit +
        1.96 *
        prediction_results_final$se.fit
    )
  )


#Check predictions

head(
  Temp_predictions_final_df
)


#Plot temperature responses

flower_visitor_colours <- c(
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
    data = analysis_data_quad2,
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
    values = flower_visitor_colours
  ) +
  
  scale_fill_manual(
    values = flower_visitor_colours
  ) +
  
  labs(
    x = "Temperature (°C)",
    y = "Predicted number of visits per 20-minute survey"
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    strip.text = element_text(
      face = "bold"
    ),
    
    legend.position = "none"
  )


Temperature_plot


#Final outputs

summary(
  ModelTempGroup_final
)

TempGroup_final_Anova

AIC(
  ModelTempGroup_final
)

ModelTempGroup_final$sdr$pdHess

temp_gradient_results

TempGroup_uniformity
TempGroup_dispersion
TempGroup_zero_inflation

Temperature_plot
