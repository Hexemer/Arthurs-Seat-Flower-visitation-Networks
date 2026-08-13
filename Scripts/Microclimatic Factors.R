# MICROCLIMATIC FACTORS
# Which environmental variables best predict flower visitation?


# Packages

library(dplyr)
library(tidyr)
library(glmmTMB)
library(car)
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


#Fix dates

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
    
    Abundance = first(
      Abundance
    ),
    
    Minutes = first(
      Minutes
    ),
    
    Temperature..C.degrees. = first(
      Temperature..C.degrees.
    ),
    
    Solar.Power..W.m.2. = first(
      Solar.Power..W.m.2.
    ),
    
    Humidity....RH. = first(
      Humidity....RH.
    ),
    
    Wind.Speed..m.s. = first(
      Wind.Speed..m.s.
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


#Prep data

survey_environment2 <- survey_data %>%
  filter(
    !(Plant.ID %in% c(
      "6",
      "13"
    )),
    
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
    Plant.ID = droplevels(
      factor(Plant.ID)
    ),
    
    Date = factor(
      Date
    ),
    
    Time.Period = factor(
      Time.Period
    ),
    
    Temp_c =
      Temperature..C.degrees. -
      mean(
        Temperature..C.degrees.,
        na.rm = TRUE
      ),
    
    Solar_c =
      Solar.Power..W.m.2. -
      mean(
        Solar.Power..W.m.2.,
        na.rm = TRUE
      ),
    
    Humidity_c =
      Humidity....RH. -
      mean(
        Humidity....RH.,
        na.rm = TRUE
      ),
    
    Wind_c =
      Wind.Speed..m.s. -
      mean(
        Wind.Speed..m.s.,
        na.rm = TRUE
      )
  )


#Check final data

nrow(
  survey_environment2
)

table(
  survey_environment2$Date
)

table(
  survey_environment2$Plant.ID
)

summary(
  survey_environment2$Number.of.Flowers
)


# Null model for environmental predictors

#Null model contains the same controls as the final model,
#but excludes all environmental predictors.

Environment_null <- glmmTMB(
  Visits ~
    Plant.ID +
    offset(
      log(Number.of.Flowers)
    ) +
    (1 | Date),
  
  family = nbinom2,
  data = survey_environment2
)

summary(
  Environment_null
)


# Linear abiotic model

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

summary(
  Environment_linear
)

car::Anova(
  Environment_linear,
  type = 2
)


# Test quadratic temperature response

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


# Test quadratic solar-radiation response

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


# Test quadratic humidity response

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


# Final environmental model

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

summary(
  Environment_final2
)


#Overall fixed-effect tests

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

VarCorr(
  Environment_final2
)


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


# Model diagnostics

set.seed(123)

Environment_diagnostics <- simulateResiduals(
  fittedModel = Environment_final2,
  n = 1000
)

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


#Four-panel diagnostic figure

par(
  mfrow = c(2, 2),
  mar = c(5, 5, 3, 2),
  cex = 0.8
)


#Panel 1: QQ residual diagnostic

plotQQunif(
  Environment_diagnostics
)


#Panel 2: Residuals against fitted values

plotResiduals(
  Environment_diagnostics,
  form = fitted(
    Environment_final2
  )
)


#Panel 3: Dispersion diagnostic

testDispersion(
  Environment_diagnostics,
  plot = TRUE
)


#Panel 4: Zero-inflation diagnostic

testZeroInflation(
  Environment_diagnostics,
  plot = TRUE
)


#Reset plotting layout

par(
  mfrow = c(1, 1)
)


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


#Compare null model with final environmental model

anova(
  Environment_null,
  Environment_final2
)

AIC(
  Environment_null,
  Environment_final2
)


# Final outputs

summary(
  Environment_final2
)

Environment_final_Anova

AIC(
  Environment_final2
)

VarCorr(
  Environment_final2
)

Environment_final2$sdr$pdHess

Environment_dispersion
Environment_zero_inflation
Environment_uniformity

anova(
  Environment_final2,
  Environment_final_time
)

AIC(
  Environment_final2,
  Environment_final_time
)

anova(
  Environment_null,
  Environment_final2
)

AIC(
  Environment_null,
  Environment_final2
)

