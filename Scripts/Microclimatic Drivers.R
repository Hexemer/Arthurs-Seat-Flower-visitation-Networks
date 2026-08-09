##################################
#Microclimatic Drivers
##################################

# Aim 4:
# Determine the extent to which microclimatic variables
# correlate with pollinator visitation and whether they
# explain the temporal pattern in visitation.


source("Overall Survey and Data.R")


#Create survey-level environmental dataset

survey_environment <- pollinators %>%
  
  mutate(
    Date = as.Date(
      Date,
      format = "%d/%m/%Y"
    ),
    
    # Correct confirmed date-entry errors
    
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
    
    Solar.Power..W.m.2. = first(
      Solar.Power..W.m.2.
    ),
    
    Temperature..C.degrees. = first(
      Temperature..C.degrees.
    ),
    
    Humidity....RH. = first(
      Humidity....RH.
    ),
    
    Wind.Speed..m.s. = first(
      Wind.Speed..m.s.
    ),
    
    .groups = "drop"
  )


#Prepare environmental analysis dataset

survey_environment2 <- survey_environment %>%
  
  filter(
    !(Plant.ID %in% c(
      "6",
      "13"
    )),
    
    !is.na(Visits),
    
    !is.na(
      Number.of.Flowers
    ),
    
    Number.of.Flowers > 0,
    
    !is.na(
      Temperature..C.degrees.
    ),
    
    !is.na(
      Solar.Power..W.m.2.
    ),
    
    !is.na(
      Humidity....RH.
    ),
    
    !is.na(
      Wind.Speed..m.s.
    ),
    
    !is.na(Date),
    
    !is.na(
      Plant.ID
    ),
    
    !is.na(
      Time.Period
    )
  ) %>%
  
  mutate(
    Plant.ID = droplevels(
      factor(
        Plant.ID
      )
    ),
    
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


#Check final environmental dataset

nrow(
  survey_environment2
)

table(
  survey_environment2$Date
)

table(
  survey_environment2$Plant.ID
)

table(
  survey_environment2$Time.Period
)

summary(
  survey_environment2$Number.of.Flowers
)

summary(
  survey_environment2$
    Temperature..C.degrees.
)

summary(
  survey_environment2$
    Solar.Power..W.m.2.
)

summary(
  survey_environment2$
    Humidity....RH.
)

summary(
  survey_environment2$
    Wind.Speed..m.s.
)


#Linear environmental model

# Plant.ID controls for differences in visitation among plant taxa.
# Number of flowers is included as an offset so visitation is
# modelled relative to floral availability.
# Date is included as a random intercept to account for
# observations collected on the same survey day.

Environment_linear <- glmmTMB(
  Visits ~
    Temp_c +
    Solar_c +
    Humidity_c +
    Wind_c +
    Plant.ID +
    offset(
      log(
        Number.of.Flowers
      )
    ) +
    (1 | Date),
  
  family = nbinom2,
  
  data =
    survey_environment2
)

summary(
  Environment_linear
)

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
      log(
        Number.of.Flowers
      )
    ) +
    (1 | Date),
  
  family = nbinom2,
  
  data =
    survey_environment2
)


anova(
  Environment_linear,
  Environment_temp_quad
)

AIC(
  Environment_linear,
  Environment_temp_quad
)

# A significant likelihood-ratio test and/or lower AIC
# supports retaining the quadratic temperature term.


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
      log(
        Number.of.Flowers
      )
    ) +
    (1 | Date),
  
  family = nbinom2,
  
  data =
    survey_environment2
)


anova(
  Environment_linear,
  Environment_solar_quad
)

AIC(
  Environment_linear,
  Environment_solar_quad
)

# If this does not improve model fit,
# solar radiation is retained as a linear effect.


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
      log(
        Number.of.Flowers
      )
    ) +
    (1 | Date),
  
  family = nbinom2,
  
  data =
    survey_environment2
)


anova(
  Environment_linear,
  Environment_humidity_quad
)

AIC(
  Environment_linear,
  Environment_humidity_quad
)

# If this does not improve model fit,
# humidity is retained as a linear effect.


#Final environmental model

# Temperature is retained as a quadratic response.
# Solar radiation, humidity and wind speed are retained
# as linear effects.

Environment_final2 <- glmmTMB(
  Visits ~
    Temp_c +
    I(Temp_c^2) +
    Solar_c +
    Humidity_c +
    Wind_c +
    Plant.ID +
    offset(
      log(
        Number.of.Flowers
      )
    ) +
    (1 | Date),
  
  family = nbinom2,
  
  data =
    survey_environment2
)


summary(
  Environment_final2
)


Environment_final_Anova <- car::Anova(
  Environment_final2,
  type = 2
)

Environment_final_Anova


#Does time period explain anything beyond microclimate?

# This model is identical to the final environmental model
# except that Time.Period is added.

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
      log(
        Number.of.Flowers
      )
    ) +
    (1 | Date),
  
  family = nbinom2,
  
  data =
    survey_environment2
)


anova(
  Environment_final2,
  Environment_final_time
)

AIC(
  Environment_final2,
  Environment_final_time
)

# If adding Time.Period does not significantly improve fit,
# the apparent temporal pattern in visitation is largely
# accounted for by the measured environmental variables.


#Check Date random effect

VarCorr(
  Environment_final2
)


#Fit identical model without Date

Environment_final_no_date <- glmmTMB(
  Visits ~
    Temp_c +
    I(Temp_c^2) +
    Solar_c +
    Humidity_c +
    Wind_c +
    Plant.ID +
    offset(
      log(
        Number.of.Flowers
      )
    ),
  
  family = nbinom2,
  
  data =
    survey_environment2
)


AIC(
  Environment_final_no_date,
  Environment_final2
)

# Use the Date variance estimate, convergence and AIC
# to assess whether the random intercept should be retained.


#Model diagnostics

Environment_diagnostics <- simulateResiduals(
  fittedModel =
    Environment_final2,
  
  n = 1000
)


plot(
  Environment_diagnostics
)


Environment_uniformity <- testUniformity(
  Environment_diagnostics
)

Environment_dispersion <- testDispersion(
  Environment_diagnostics
)

Environment_zero_inflation <- testZeroInflation(
  Environment_diagnostics
)


Environment_uniformity

Environment_dispersion

Environment_zero_inflation


#Check convergence

Environment_final2$sdr$pdHess

# TRUE supports successful model convergence.


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

# Prefer the distribution with the lower AIC only if
# convergence and DHARMa diagnostics are satisfactory.


#Environmental descriptive statistics

environment_summary <- survey_environment2 %>%
  summarise(
    Temperature_Mean = mean(
      Temperature..C.degrees.,
      na.rm = TRUE
    ),
    
    Temperature_SD = sd(
      Temperature..C.degrees.,
      na.rm = TRUE
    ),
    
    Temperature_Min = min(
      Temperature..C.degrees.,
      na.rm = TRUE
    ),
    
    Temperature_Max = max(
      Temperature..C.degrees.,
      na.rm = TRUE
    ),
    
    Solar_Mean = mean(
      Solar.Power..W.m.2.,
      na.rm = TRUE
    ),
    
    Solar_SD = sd(
      Solar.Power..W.m.2.,
      na.rm = TRUE
    ),
    
    Humidity_Mean = mean(
      Humidity....RH.,
      na.rm = TRUE
    ),
    
    Humidity_SD = sd(
      Humidity....RH.,
      na.rm = TRUE
    ),
    
    Wind_Mean = mean(
      Wind.Speed..m.s.,
      na.rm = TRUE
    ),
    
    Wind_SD = sd(
      Wind.Speed..m.s.,
      na.rm = TRUE
    )
  )

environment_summary


#Final Aim 4 outputs

summary(
  Environment_final2
)

Environment_final_Anova

anova(
  Environment_final2,
  Environment_final_time
)

AIC(
  Environment_final2,
  Environment_final_time
)

VarCorr(
  Environment_final2
)

AIC(
  Environment_final_no_date,
  Environment_final2
)

Environment_final2$sdr$pdHess

Environment_uniformity

Environment_dispersion

Environment_zero_inflation

environment_summary
