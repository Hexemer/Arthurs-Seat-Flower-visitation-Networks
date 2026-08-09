##################################
#Visitation Time and Site
##################################

getwd()
source("Overall Survey and Data.R")
list.files()


#Prepare dataset

# Plant IDs 6 and 13 received no pollinator visits in any survey.
# They are retained in descriptive summaries but excluded from
# inferential count models because they provide no within-plant
# variation in visitation.

survey_model <- survey_data %>%
  filter(
    !(Plant.ID %in% c("6", "13"))
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


#Check final analysis structure

nrow(
  survey_model
)

table(
  survey_model$Date
)

table(
  survey_model$Plant.ID
)

table(
  survey_model$Time.Period
)

table(
  survey_model$Site
)


#Block or no block

# Block was determined as the first two weeks and the second two weeks.
# This was determined using a negative binomial GLM.

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

# Addition of Block was not significant.


#Time-period and site model

# Inclusion of Date and Plant.ID as random effects.

ModelTimeSite <- glmmTMB(
  Visits ~
    Time.Period +
    Site +
    (1 | Date) +
    (1 | Plant.ID),
  family = nbinom2,
  data = survey_model
)


#Model diagnostics for nbinom2

simulationOutput_TimeSite_nbinom2 <- simulateResiduals(
  fittedModel = ModelTimeSite,
  n = 1000
)

plot(
  simulationOutput_TimeSite_nbinom2
)

Uniformity_TimeSite_nbinom2 <- testUniformity(
  simulationOutput_TimeSite_nbinom2
)

Dispersion_TimeSite_nbinom2 <- testDispersion(
  simulationOutput_TimeSite_nbinom2
)

ZeroInflation_TimeSite_nbinom2 <- testZeroInflation(
  simulationOutput_TimeSite_nbinom2
)

Uniformity_TimeSite_nbinom2
Dispersion_TimeSite_nbinom2
ZeroInflation_TimeSite_nbinom2


#Compare nbinom1 and nbinom2

# nbinom1: variance increases linearly with the mean.
# nbinom2: variance increases quadratically with the mean.

ModelTimeSite_nbinom1 <- update(
  ModelTimeSite,
  family = nbinom1
)

ModelTimeSite_nbinom1$sdr$pdHess

AIC(
  ModelTimeSite,
  ModelTimeSite_nbinom1
)


#Model diagnostics for nbinom1

simulationOutput_TimeSite_nbinom1 <- simulateResiduals(
  fittedModel = ModelTimeSite_nbinom1,
  n = 1000
)

plot(
  simulationOutput_TimeSite_nbinom1
)

Uniformity_TimeSite_nbinom1 <- testUniformity(
  simulationOutput_TimeSite_nbinom1
)

Dispersion_TimeSite_nbinom1 <- testDispersion(
  simulationOutput_TimeSite_nbinom1
)

ZeroInflation_TimeSite_nbinom1 <- testZeroInflation(
  simulationOutput_TimeSite_nbinom1
)

Uniformity_TimeSite_nbinom1
Dispersion_TimeSite_nbinom1
ZeroInflation_TimeSite_nbinom1

# nbinom1 was chosen as it converges cleanly;
# pdHess is TRUE; diagnostics are satisfactory;
# AIC is meaningfully lower.


#Final model

ModelTimeSite_final <- ModelTimeSite_nbinom1


#Test overall fixed effects

TimeSite_Anova <- car::Anova(
  ModelTimeSite_final,
  type = 2
)

TimeSite_Anova

# Time.Period:
# Does visitation vary throughout the day?

# Site:
# Does visitation differ between the two sites?


#Estimated visitation by time period

# emmeans tells you what visitation the model predicts for each group
# after accounting for the other variables in the model.

Time_emmeans <- emmeans(
  ModelTimeSite_final,
  pairwise ~ Time.Period,
  type = "response",
  adjust = "tukey"
)

Time_emmeans

# Predicted visits per survey for each period

Time_emmeans$emmeans

# Tukey-adjusted pairwise comparisons

Time_emmeans$contrasts


#Estimated visitation by site

Site_emmeans <- emmeans(
  ModelTimeSite_final,
  pairwise ~ Site,
  type = "response"
)

Site_emmeans

# Predicted visits per survey for each site

Site_emmeans$emmeans

# Site comparison as a response-scale ratio

Site_emmeans$contrasts


#Descriptive summaries

survey_data %>%
  group_by(
    Time.Period
  ) %>%
  summarise(
    Surveys = n(),
    Total_Visits = sum(
      Visits
    ),
    Mean_Visits = mean(
      Visits
    ),
    Median_Visits = median(
      Visits
    ),
    Zero_Visit_Surveys = sum(
      Visits == 0
    ),
    .groups = "drop"
  )

survey_data %>%
  group_by(
    Site
  ) %>%
  summarise(
    Surveys = n(),
    Total_Visits = sum(
      Visits
    ),
    Mean_Visits = mean(
      Visits
    ),
    Median_Visits = median(
      Visits
    ),
    Zero_Visit_Surveys = sum(
      Visits == 0
    ),
    .groups = "drop"
  )

survey_data %>%
  group_by(
    Plant.ID,
    Site
  ) %>%
  summarise(
    Surveys = n(),
    Total_Visits = sum(
      Visits
    ),
    Mean_Visits = mean(
      Visits
    ),
    .groups = "drop"
  )


#Check random effects and convergence

VarCorr(
  ModelTimeSite_final
)

ModelTimeSite_final$sdr$pdHess

# TRUE indicates a positive-definite Hessian and supports
# successful model convergence.


#Check the value of Date as a random effect

ModelTimeSite_no_date <- glmmTMB(
  Visits ~
    Time.Period +
    Site +
    (1 | Plant.ID),
  family = nbinom1,
  data = survey_model
)

AIC(
  ModelTimeSite_no_date,
  ModelTimeSite_final
)

# AIC was better in the model which included Date as a random effect.


#Final outputs

summary(
  ModelTimeSite_final
)

TimeSite_Anova

Time_emmeans$emmeans
Time_emmeans$contrasts

Site_emmeans$emmeans
Site_emmeans$contrasts

VarCorr(
  ModelTimeSite_final
)

AIC(
  ModelTimeSite_final
)

ModelTimeSite_final$sdr$pdHess

plot(
  simulationOutput_TimeSite_nbinom1
)

Uniformity_TimeSite_nbinom1
Dispersion_TimeSite_nbinom1
ZeroInflation_TimeSite_nbinom1

