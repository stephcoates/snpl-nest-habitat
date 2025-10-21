# Modeling nests vs. random sites with logistic regression
# 1. define the null model for all data, summer subset, and fall subset
# 2. test use of linear or quadratic terms for covariates
# 3. create sets of global models with uncorrelated variables
# 4. use stepwise function to test submodels of each global model, compare AIC
# 5. check model fit, VIF; remove poor fitting models and check for high VIF
# 6. compare remaining candidate models

# load packages and functions and data
source('R/packages.R')
source('R/functions_glm.R')

all_df <- read_csv('data/glm_alldat.csv')
summer_df <- read_csv('data/glm_breedingdat.csv')
fall_df <- read_csv('data/glm_falldat.csv')

# 1. Null models ----
nullmod.noR <- glmmTMB(Outcome ~1,
                       data = all_df, 
                       family = binomial) #AIC = 236.3
nullmod.nesting <- glmmTMB(Outcome ~1, 
                           data = summer_df,
                           family = binomial) #AIC = 112.9
nullmod.fall <- glmmTMB(Outcome ~1, 
                        data = fall_df,
                        family = binomial) #AIC = 112.9
# specify null model with random effect
# nullmod <- glmer(Outcome ~(1|PairID),
#                data = all_df,
#                family = binomial) 
# did not end up using because model fit was not improved with random effect

# 2. Test linear and quadratic terms of each covariate ----

# Functional groups and overall percent cover covariates
# variables to test
all_variables <- c("Pioneer_Stabilizers", "Pioneer_Builders", "Secondary_Stabilizers", 
                   "Half", "DeadVegManual", "Viewshed_sc", "NonVeg", "Two", "TwoVeg")

# check we are only using variables that exist in the dataset
all_variables <- intersect(all_variables, names(all_df))

# pass variables to the function
aic_results_all <- run_glmer_models(all_variables, all_df)
aic_results_summer <- run_glmer_models(all_variables, summer_df)
aic_results_fall <- run_glmer_models(all_variables, fall_df)

# inspect
print(aic_results_all) 
# AIC values and term to use in models
# Quadratic: half^2(187.0), Two^2(219.1), TwoVeg(223.9), deadvegmanual^2(224.8), 
#             PioneerStabilizers^2(230.2) 
# Linear: NonVeg (215.6), Secondary Stabilizers(235.1), Viewshed(231.5),
# Cusp:  pioneer builders^2(229.9)
print(aic_results_summer)
# Quadratic: Half, Two, TWoVeg, Secondary Stabilizers, DeadVegManual
# Linear: NonVeg, Viewshed_sc, Pioneer Stabilizer
# Cusp: Pioneer Builder
print(aic_results_fall)
# Quadratic: Half, Two, TwoVeg, Secondary Stabilizers, DeadVegManual, Pioneer Stabilizers
# Linear: NonVeg, Viewshed_sc
# Cusp: Pioneer Builder

# Individual species
# covariates not tested in plant functional group section
all_variables.overall <- c("SR", "BP", "IP", "BS", "VE", "SG", "DG", "WO", "SL")

# Ensure we are only using variables that exist in the dataset
all_variables.overall <- intersect(all_variables.overall, names(all_df))

# Now pass only these variables to the function
aic_results_all.overall <- run_glmer_models(all_variables.overall, all_df)
aic_results_summer.overall <- run_glmer_models(all_variables.overall, summer_df)
aic_results_fall.overall <- run_glmer_models(all_variables.overall, fall_df)

print(aic_results_all.overall) 
# use quadratic for SR(231.6), BS(230.1), 
# use linear for IP(222.0), WO (207.8), SL (236.0)
print(aic_results_summer.overall)
# quadratic: SR
# linear: WO, IP, SL, BS
print(aic_results_fall.overall)
# quadratic: 
# linear: WO, IP, SR, SL
# cusp: BS

# 3. Global Models and stepwise selection to get candidate models ----
# 3a. All data: Functional Groups ----
# groups species by the functional group they belong to. no individual species
# linear vs. quadratic terms decided from the results of the function; add 
# random effect back in later. 

# Check vif for covariates in global model
fullmod.lm <-lm(Outcome ~ poly(Half,2) + poly(Pioneer_Stabilizers,2)
                + poly(Two,2) + Viewshed_sc + NonVeg + Secondary_Stabilizers
                + Pioneer_Builders + poly(DeadVegManual,2) + poly(TwoVeg,2),
                data = all_df)
vif(fullmod.lm)
# TwoVeg and Two have high VIF; TwoVeg may be slightly of more interest because
# Non-veg cover is also included here (Non-veg + TwoVeg = Two). After taking out 
# Two, TwoVeg still has VIF > 5 (5.115). Same thing if Two is left instead;
# VIF > 5 (5.511). Both removed because they exceeded cutoff value of 5.

#stepwise function to select candidate models
fullmod <- glmmTMB(Outcome ~ Object_Type + poly(Half,2) + poly(Pioneer_Stabilizers,2)
                   + Viewshed_sc + NonVeg + Secondary_Stabilizers
                   + Pioneer_Builders + poly(DeadVegManual,2),
                   data = all_df, family = binomial)
fullmod.step <- step(nullmod.noR, list(lower=formula(nullmod.noR), 
                                       upper=formula(fullmod)), direction = "both")
summary(fullmod.step) 

# exploring the top ranked models according to stepwise function 
fullmod.1 <- glmmTMB(Outcome ~ Object_Type + poly(Half,2) + Pioneer_Builders 
                     + poly(Pioneer_Stabilizers,2), 
                     data = all_df, family = binomial)
summary(fullmod.1) #104.6

fullmod.2 <- glmmTMB(Outcome ~ Object_Type + poly(Half,2) + Pioneer_Builders 
                     + poly(Pioneer_Stabilizers,2) + Viewshed_sc, 
                     data = all_df, family = binomial)
summary(fullmod.2) #105.5

fullmod.3 <- glmmTMB(Outcome ~ Object_Type + poly(Half,2) + Pioneer_Builders 
                     + poly(Pioneer_Stabilizers,2) + Secondary_Stabilizers, 
                      data = all_df, family = binomial)
summary(fullmod.3) #106.0

fullmod.4 <- glmmTMB(Outcome ~ Object_Type + poly(Half,2) + Pioneer_Builders 
                     + poly(Pioneer_Stabilizers,2) + NonVeg, 
                      data = all_df, family = binomial)
summary(fullmod.4) #106.2

fullmod.5 <- glmmTMB(Outcome ~ Object_Type + poly(Half,2) + Pioneer_Builders 
                     + poly(Pioneer_Stabilizers,2) + poly(DeadVegManual,2), 
                      data = all_df, family = binomial)
summary(fullmod.5) #107.4

fullmod.6 <- glmmTMB(Outcome ~ Object_Type + Pioneer_Builders 
                     + poly(Pioneer_Stabilizers,2), 
                     data = all_df, family = binomial)
summary(fullmod.6) #119.3


# create list of models
models.full.fg <- list(
  fullmod.1 = fullmod.1,
  fullmod.2 = fullmod.2,
  fullmod.3 = fullmod.3,
  fullmod.4 = fullmod.4,
  fullmod.5 = fullmod.5
  #fullmod.6 = fullmod.6
)

# check all for goodness of fit
fit_output <- run_all_dharma_diagnostics(models.full.fg)
fit_output$summary
# all look ok here; fullmod.6 has p-value of 0.07 (on the cusp of being an issue)

# create candidate model table
cand.mods.full.fg <- create_candidate_model_table(models.full.fg)

# Selected top model; lowest AIC, most parsimonious
finalmod.all.fg <- fullmod.1

# Test adding random effect into selected model:
# point ID is nested within the pair ID 
# including PairID accounts for resampling the same territory. 
# can't use PointID because it perfectly predicts the outcome; nest is always nest

# selected model
REmod.all.fg <- glmmTMB(Outcome ~ (1|PairID) + Object_Type + poly(Half, 2)
                        + poly(Pioneer_Stabilizers,2) + Pioneer_Builders,
                        family = binomial, data = all_df)
summary(REmod.all.fg)
# adding PairID random effect does not change coefficient estimates, explains no 
# additional variance

# Random effect conclusions:
# The effect is negligible; does not explain the variation in outcome. 
# Having a random effect also creates convergence problems because it explains 
# no variance.

# 3b. All data: Individual Species or Functional Groups ----
# includes individual species or functional groups, whichever create the best-
# fitting model
# linear vs. quadratic terms decided from the results of the function; add 
# random effect back in later. 

# From correlation matrix we know already:
# BS is perfectly correlated with Pioneer Stabilizers; eliminated BS
# SR and Pioneer Builders has a .97 correlation coefficient; eliminated SR
# WO and NonVeg have 0.9 corr coef; eliminated NonVeg
# VE and Pioneer Stabilizers has 0.77 corr coeff; eliminated VE
# Two and TwoVeg were a known issue from before; removed both

# Check vif for covariates in global model after removing perfectly correlated variables
fullmod.o.lm <-lm(Outcome ~ poly(Half,2)
                + Viewshed_sc + poly(DeadVegManual,2)
                + poly(Pioneer_Stabilizers,2) + Pioneer_Builders
                + IP + WO + SL,
                data = all_df)
vif(fullmod.o.lm)
# Secondary Stabilizers was causing issues because the combo of individual species
# included here equals the Secondary stabilizers. Removed secondary stabilizers.
# After removal, no issues w/ VIF, all <5


#stepwise function to select candidate models
fullmod.o <- glmmTMB(Outcome ~ Object_Type + poly(Half,2)
                   + Viewshed_sc + poly(DeadVegManual,2)
                   + poly(Pioneer_Stabilizers,2) + Pioneer_Builders
                   + IP + WO + SL,
                   data = all_df, family = binomial)
fullmod.o.step <- step(nullmod.noR, list(lower=formula(nullmod.noR), 
                                       upper=formula(fullmod.o)), direction = "both")
summary(fullmod.o.step) 

# exploring the top ranked models according to stepwise function 
fullmod.o.1 <- glmmTMB(Outcome ~ Object_Type + poly(Half,2) + Pioneer_Builders 
                     + poly(Pioneer_Stabilizers,2) + IP + WO, 
                     data = all_df, family = binomial)
summary(fullmod.o.1) #98.0

fullmod.o.2 <- glmmTMB(Outcome ~ Object_Type + poly(Half,2) + Pioneer_Builders 
                     + poly(Pioneer_Stabilizers,2) + IP + WO + SL, 
                     data = all_df, family = binomial)
summary(fullmod.o.2) #99.6; adding covariate does not improve fit by >2AIC

fullmod.o.3 <- glmmTMB(Outcome ~ Object_Type + poly(Half,2) + Pioneer_Builders 
                       + poly(Pioneer_Stabilizers,2) + IP, 
                     data = all_df, family = binomial)
summary(fullmod.o.3) #99.8; this is most parsimonious model w/in 2 AIC of top

fullmod.o.4 <- glmmTMB(Outcome ~ Object_Type + poly(Half,2) + Pioneer_Builders 
                     + poly(Pioneer_Stabilizers,2) + IP + WO + Viewshed_sc, 
                     data = all_df, family = binomial)
summary(fullmod.o.4) #99.9; adding covariate does not improve fit by >2AIC

fullmod.o.5 <- glmmTMB(Outcome ~ Object_Type + poly(Half,2) + Pioneer_Builders 
                     + poly(Pioneer_Stabilizers,2) + IP + WO 
                     + poly(DeadVegManual,2), 
                     data = all_df, family = binomial)
summary(fullmod.o.5) # 100.17; adding covariate does not improve fit by >2AIC

fullmod.o.6 <- glmmTMB(Outcome ~ Object_Type + poly(Half,2) + Pioneer_Builders 
                     + poly(Pioneer_Stabilizers,2) + WO, 
                     data = all_df, family = binomial)
summary(fullmod.o.6) # 102.55; removing IP increases AIC by >2

fullmod.o.7 <- glmmTMB(Outcome ~ Object_Type + Pioneer_Builders 
                       + poly(Pioneer_Stabilizers,2) + IP + WO, 
                       data = all_df, family = binomial)
summary(fullmod.o.7) # 111.3
# the most parsimonious model w/in 2 AIC of top-ranked model is: fullmod.o.3


# create list of models
models.full.o <- list(
  fullmod.o.1 = fullmod.o.1,
  fullmod.o.2 = fullmod.o.2,
  fullmod.o.3 = fullmod.o.3,
  fullmod.o.4 = fullmod.o.4,
  fullmod.o.5 = fullmod.o.5
  #fullmod.o.6 = fullmod.o.6, 
  #fullmod.o.7 = fullmod.o.7
)

# check all for goodness of fit
fit_output <- run_all_dharma_diagnostics(models.full.o)
fit_output$summary
# all look ok here

# create candidate model table
cand.mods.full.o <- create_candidate_model_table(models.full.o)

# Selected top model; most parsimonious w/in 2 AIC of top-ranked model
finalmod.all.o <- fullmod.o.3

# Test adding random effect into selected model:
# point ID is nested within the pair ID 
# including PairID accounts for resampling the same territory. 
# can't use PointID because it perfectly predicts the outcome; nest is always nest

# selected model
REmod.all.o <- glmmTMB(Outcome ~ (1|PairID) + Object_Type + poly(Half, 2)
                        + poly(Pioneer_Stabilizers,2) + Pioneer_Builders + IP,
                        family = binomial, data = all_df)
summary(REmod.all.o)
# adding PairID random effect does not change coefficient estimates or p-value, 
# explains no additional variance

# Random effect conclusions:
# The effect is negligible; does not explain the variation in outcome. 
# Having a random effect also creates convergence problems because it explains 
# no variance.


# 3c. Summer subset: Functional Groups ----

# From correlation matrix: 
# Pioneer Builders and Secondary Stabilizers are correlated with Half
# Two and TwoVeg are correlated with Half and PionBuild and SecondStab
# Half is correlated w/ Two and TwoVeg

# Check vif for covariates in global model 1
summermod.lm <-lm(Outcome ~ poly(Half,2) + Pioneer_Stabilizers 
                  + Viewshed_sc + NonVeg 
                  + poly(DeadVegManual,2),
                data = summer_df)
vif(summermod.lm) # all ok

# Check vif for covariates in global model 2
summermod.lm <-lm(Outcome ~  Pioneer_Stabilizers 
                  + Viewshed_sc + NonVeg + poly(Secondary_Stabilizers,2)
                  + Pioneer_Builders + poly(DeadVegManual,2),
                  data = summer_df)
vif(summermod.lm) # all ok

# Check vif for covariates in global model 3
summermod.lm <-lm(Outcome ~ poly(Two,2) + Pioneer_Stabilizers 
                  + Viewshed_sc + NonVeg 
                  + poly(DeadVegManual,2),
                  data = summer_df)
vif(summermod.lm) #all ok


#Global model 1. stepwise function to select candidate models. 
summermod.g1 <- glmmTMB(Outcome ~ Object_Type + poly(Half,2) + Pioneer_Stabilizers 
                        + Viewshed_sc + NonVeg 
                        + poly(DeadVegManual,2),
                   data = summer_df, family = binomial)
summermod.step.g1 <- step(nullmod.nesting, list(lower=formula(nullmod.nesting), 
                                       upper=formula(summermod.g1)), direction = "both")
summary(summermod.step.g1) # AIC 61.5

#Global model 2. stepwise function to select candidate models. 
summermod.g2 <- glmmTMB(Outcome ~ Object_Type + Pioneer_Stabilizers 
                        + Viewshed_sc + NonVeg + poly(Secondary_Stabilizers,2)
                        + Pioneer_Builders + poly(DeadVegManual,2),
                        data = summer_df, family = binomial)
summermod.step.g2 <- step(nullmod.nesting, list(lower=formula(nullmod.nesting), 
                                                upper=formula(summermod.g2)), direction = "both")
summary(summermod.step.g2) # AIC 61.5

#Global model 3. stepwise function to select candidate models. 
summermod.g3 <- glmmTMB(Outcome ~ Object_Type + poly(Two,2) + Pioneer_Stabilizers 
                        + Viewshed_sc + NonVeg 
                        + poly(DeadVegManual,2),
                        data = summer_df, family = binomial)
summermod.step.g3 <- step(nullmod.nesting, list(lower=formula(nullmod.nesting), 
                                                upper=formula(summermod.g3)), direction = "both")
summary(summermod.step.g3) 

# exploring the top ranked models according to stepwise function 
summermod.1 <- glmmTMB(Outcome ~ Object_Type + poly(Half,2) + Pioneer_Builders,
                     data = summer_df, family = binomial)
summary(summermod.1)
# AIC 50.6

summermod.2 <- glmmTMB(Outcome ~ Object_Type + poly(Half,2) + Pioneer_Builders
                       + Viewshed_sc,
                       data = summer_df, family = binomial)
summary(summermod.2)
# AIC 52.3

summermod.3 <- glmmTMB(Outcome ~ Object_Type + poly(Half,2) + Pioneer_Builders
                       + Secondary_Stabilizers,
                       data = summer_df, family = binomial)
summary(summermod.3)
# AIC 52.4

summermod.4 <- glmmTMB(Outcome ~ Object_Type + poly(Half,2) + Pioneer_Builders
                       + NonVeg,
                       data = summer_df, family = binomial)
summary(summermod.4)
# AIC 52.5

summermod.5 <- glmmTMB(Outcome ~ Object_Type + poly(Half,2) + Pioneer_Builders
                       + poly(Pioneer_Stabilizers, 2),
                       data = summer_df, family = binomial)
summary(summermod.5)
# AIC 53.1

summermod.6 <- glmmTMB(Outcome ~ Object_Type + poly(Half,2) + Pioneer_Builders
                       + poly(DeadVegManual, 2),
                       data = summer_df, family = binomial)
summary(summermod.6)
# AIC 54.6

summermod.7 <- glmmTMB(Outcome ~ poly(Half,2) + Pioneer_Builders,
                       data = summer_df, family = binomial)
summary(summermod.7)
# AIC 56.9

# create list of models
models.summer.fg <- list(
  summermod.1 = summermod.1,
  summermod.2 = summermod.2,
  summermod.3 = summermod.3,
  summermod.4 = summermod.4,
  summermod.5 = summermod.5
  #summermod.6 = summermod.6
)

# check all for goodness of fit
fit_output <- run_all_dharma_diagnostics(models.summer.fg)
fit_output$summary
# all look ok here

# create candidate model table
cand.mods.summer.fg <- create_candidate_model_table(models.summer.fg)

# Selected top model; lowest AIC, most parsimonious
finalmod.summer.fg <- summermod.1
summary(finalmod.summer.fg)

# Test adding random effect into selected model:
# point ID is nested within the pair ID 
# including PairID accounts for resampling the same territory. 
# can't use PointID because it perfectly predicts the outcome; nest is always nest

# selected model
REmod.summer.fg <- glmmTMB(Outcome ~ (1|PairID) + Object_Type + poly(Half, 2)
                       + Pioneer_Builders,
                       family = binomial, data = summer_df)
summary(REmod.summer.fg)
# adding PairID random effect does not change coefficient estimates or p-value, 
# explains no additional variance

# Random effect conclusions:
# The effect is negligible; does not explain the variation in outcome. 

# 3d. Summer subset: Individual Species or Functional Groups ----
# includes individual species or functional groups, whichever create the best-
# fitting model
# linear vs. quadratic terms decided from the results of the function; add 
# random effect back in later. 

# From correlation matrix we know already:
# BS is perfectly correlated with Pioneer Stabilizers; eliminated BS
# SR and Pioneer Builders has a .97 correlation coefficient; eliminated SR
# WO and NonVeg have 0.8 corr coef; trying with both
# Two and TwoVeg were a known issue from before; removed both

# Check vif for covariates in global model after removing perfectly correlated variables
summermod.o.lm <-lm(Outcome ~ poly(Half,2) + poly(Secondary_Stabilizers,2) 
                  + Viewshed_sc + poly(DeadVegManual,2)
                  + poly(Pioneer_Stabilizers,2) + Pioneer_Builders
                  + IP + WO + SL + NonVeg,
                  data = summer_df)
vif(summermod.o.lm)
# no issues w/ VIF, all <5


#stepwise function to select candidate models
summermod.o <- glmmTMB(Outcome ~ Object_Type + poly(Half,2) 
                     + poly(Secondary_Stabilizers,2) + poly(DeadVegManual,2)
                     + poly(Pioneer_Stabilizers,2) + Pioneer_Builders
                     + Viewshed_sc + IP + WO + SL + NonVeg,
                     data = summer_df, family = binomial)
summermod.o.step <- step(nullmod.nesting, list(lower=formula(nullmod.nesting), 
                                         upper=formula(summermod.o)), direction = "both")
summary(summermod.o.step) 

# exploring the top ranked models according to stepwise function 
summermod.o.1 <- glmmTMB(Outcome ~ Object_Type + poly(Half,2) + Pioneer_Builders, 
                       data = summer_df, family = binomial)
summary(summermod.o.1) #50.63

summermod.o.2 <- glmmTMB(Outcome ~ Object_Type + poly(Half,2) + Pioneer_Builders 
                       + WO, 
                       data = summer_df, family = binomial)
summary(summermod.o.2) # 52.1

summermod.o.3 <- glmmTMB(Outcome ~ Object_Type + poly(Half,2) + Pioneer_Builders 
                       + Viewshed_sc, 
                       data = summer_df, family = binomial)
summary(summermod.o.3) # 52.3

summermod.o.4 <- glmmTMB(Outcome ~ Object_Type + poly(Half,2) + Pioneer_Builders 
                        + IP, 
                       data = summer_df, family = binomial)
summary(summermod.o.4) # 52.5

summermod.o.5 <- glmmTMB(Outcome ~ Object_Type + poly(Half,2) + Pioneer_Builders 
                       + NonVeg, 
                       data = summer_df, family = binomial)
summary(summermod.o.5) # 52.5

summermod.o.6 <- glmmTMB(Outcome ~ Object_Type + poly(Half,2) + Pioneer_Builders 
                       + poly(Secondary_Stabilizers,2), 
                       data = summer_df, family = binomial)
summary(summermod.o.6) # 52.6

summermod.o.7 <- glmmTMB(Outcome ~ Object_Type + poly(Half,2) + Pioneer_Builders 
                       + SL, 
                       data = summer_df, family = binomial)
summary(summermod.o.7) # 52.6

summermod.o.8 <- glmmTMB(Outcome ~ Object_Type + poly(Half,2) + Pioneer_Builders 
                         + poly(Pioneer_Stabilizers,2), 
                         data = summer_df, family = binomial)
summary(summermod.o.8) # 53.1

summermod.o.9 <- glmmTMB(Outcome ~ Object_Type + poly(Half,2) + Pioneer_Builders 
                         + poly(DeadVegManual,2), 
                         data = summer_df, family = binomial)
summary(summermod.o.9) # 54.6

# the most parsimonious model w/in 2 AIC of top-ranked model is: summermod.o.1

# create list of models
models.summer.o <- list(
  summermod.o.1 = summermod.o.1,
  summermod.o.2 = summermod.o.2,
  summermod.o.3 = summermod.o.3,
  summermod.o.4 = summermod.o.4,
  summermod.o.5 = summermod.o.5,
  summermod.o.6 = summermod.o.6,
  summermod.o.7 = summermod.o.7,
  summermod.o.8 = summermod.o.8,
  summermod.o.9 = summermod.o.9
)

# check all for goodness of fit
fit_output <- run_all_dharma_diagnostics(models.summer.o)
fit_output$summary
# all look ok here

# create candidate model table
cand.mods.summer.o <- create_candidate_model_table(models.summer.o)

# Selected top model; lowest AIC, most parsimonious
finalmod.summer.o <- summermod.o.1
summary(finalmod.summer.o)

# Test adding random effect into selected model:
# point ID is nested within the pair ID 
# including PairID accounts for resampling the same territory. 
# can't use PointID because it perfectly predicts the outcome; nest is always nest

# selected model
REmod.summer.o <- glmmTMB(Outcome ~ (1|PairID) + Object_Type + poly(Half, 2)
                           + Pioneer_Builders,
                           family = binomial, data = summer_df)
summary(REmod.summer.o)
# adding PairID random effect does not change coefficient estimates or p-value, 
# explains no additional variance

# Random effect conclusions:
# The effect is negligible; does not explain the variation in outcome. 
 

# 3e. Fall subset: Functional Groups ----
# groups species by the functional group they belong to. no individual species
# linear vs. quadratic terms decided from the results of the function; add 
# random effect back in later. 

# From correlation matrix:
# Pioneer Stabilizers and DeadVegManual are correlated 0.8; DeadVegManual has
# lower AIC; try both but not in same model
# Pionner Stabilizers and Secondary Stabilizers are no longer correlated with Half

# Check vif for covariates in global model
fallmod.lm <-lm(Outcome ~ poly(Half,2) #+ poly(Two,2) + poly(TwoVeg,2) 
                + poly(Pioneer_Stabilizers,2) + poly(Secondary_Stabilizers,2)
                + Viewshed_sc + NonVeg +
                + Pioneer_Builders + poly(DeadVegManual,2),
                data = fall_df)
vif(fallmod.lm)
# TwoVeg and Two have very high VIF; removed Two and re-ran. Two and 
# PioneerStabilizers are correlated, but having PioneerStabs gives lower AIC.
# Removed TwoVeg. All VIF ok after

#stepwise function to select candidate models
fallmod <- glmmTMB(Outcome ~ Object_Type + poly(Half,2) 
                   + poly(Secondary_Stabilizers,2) + poly(Pioneer_Stabilizers,2)
                     + Viewshed_sc + NonVeg 
                     + Pioneer_Builders, #+ poly(DeadVegManual,2),
                     data = fall_df, family = binomial)
fallmod.step <- step(nullmod.fall, list(lower=formula(nullmod.fall), 
                                             upper=formula(fallmod)), direction = "both")
summary(fallmod.step) 
# Tried with either Pioneer Stabs (AIC 48) or DeadVegManual (54); Kept PiSt

# exploring the top ranked models according to stepwise function 
fallmod.1 <- glmmTMB(Outcome ~ Object_Type + poly(Half,2) + Pioneer_Builders
                     + poly(Pioneer_Stabilizers,2) + poly(Secondary_Stabilizers,2)
                     + Viewshed_sc,
                       data = fall_df, family = binomial)
summary(fallmod.1)
# AIC 48.8

fallmod.2 <- glmmTMB(Outcome ~ Object_Type + poly(Half,2) + Pioneer_Builders
                     + poly(Pioneer_Stabilizers,2) + poly(Secondary_Stabilizers,2)
                     + Viewshed_sc + NonVeg,
                       data = fall_df, family = binomial)
summary(fallmod.2)
# AIC 49.6

fallmod.3 <- glmmTMB(Outcome ~ Object_Type + poly(Half,2) + Pioneer_Builders
                     + poly(Pioneer_Stabilizers,2)
                     + Viewshed_sc,
                       data = fall_df, family = binomial)
summary(fallmod.3)
# AIC 56.3

fallmod.4 <- glmmTMB(Outcome ~ Object_Type + poly(Half,2) + Pioneer_Builders,
                     data = fall_df, family = binomial)
summary(fallmod.4)
# AIC 77.9; for comparison with selected summer model


# create list of models
models.fall.fg <- list(
  fallmod.1 = fallmod.1,
  fallmod.2 = fallmod.2,
  fallmod.3 = fallmod.3
  #fallmod.4 = fallmod.4,
  #fallmod.5 = fallmod.5,
  #fallmod.6 = fallmod.6
)

# check all for goodness of fit
fit_output <- run_all_dharma_diagnostics(models.fall.fg)
fit_output$summary
# all look ok here

# create candidate model table
cand.mods.fall.fg <- create_candidate_model_table(models.fall.fg)

# Selected top model; lowest AIC, most parsimonious
finalmod.fall.fg <- fallmod.1
summary(finalmod.fall.fg)

# Test adding random effect into selected model:
# point ID is nested within the pair ID 
# including PairID accounts for resampling the same territory. 
# can't use PointID because it perfectly predicts the outcome; nest is always nest

# selected model
REmod.fall.o <- glmmTMB(Outcome ~ (1|PairID) + Object_Type + poly(Half,2) 
                        + Pioneer_Builders + poly(Pioneer_Stabilizers,2)
                        + poly(Secondary_Stabilizers,2) + Viewshed_sc,
                        family = binomial, data = fall_df)
summary(REmod.fall.o)
# adding PairID random effect does not change coefficient estimates, 
# has convergence issues, and explains no additional variance

# Random effect conclusions:
# The effect is negligible; does not explain the variation in outcome. 


# 3f. FALL subset: Individual Species or Functional Groups ----
# includes individual species or functional groups, whichever create the best-
# fitting model
# linear vs. quadratic terms decided from the results of the function; add 
# random effect back in later. 

print(aic_results_fall)
# Quadratic: Half, Two, TwoVeg, Secondary Stabilizers, DeadVegManual, Pioneer Stabilizers
# Linear: NonVeg, Viewshed_sc
# Cusp: Pioneer Builder
print(aic_results_fall.overall)
# quadratic: 
# linear: WO, IP, SR, SL
# cusp: BS


# From correlation matrix we know already:
# BS is perfectly correlated with Pioneer Stabilizers; eliminated BS
# SR and Pioneer Builders has a .96 correlation coefficient; eliminated SR
# WO and NonVeg have 0.93 corr coef; eliminated NonVeg
# Two and TwoVeg were a known issue from before; removed both
# From correlation matrix:
# Pioneer Stabilizers and DeadVegManual are correlated 0.8; removed DeadVeg

# Check vif for covariates in global model after removing correlated variables
fallmod.o.lm <-lm(Outcome ~ poly(Half,2) + poly(Secondary_Stabilizers,2) 
                    + Viewshed_sc + poly(Pioneer_Stabilizers,2) + Pioneer_Builders
                    + IP + WO + SL,
                    data = fall_df)
vif(fallmod.o.lm)
# no issues w/ VIF, all <5


#stepwise function to select candidate models
fallmod.o <- glmmTMB(Outcome ~ Object_Type + poly(Half,2) + poly(Secondary_Stabilizers,2) 
                     + Viewshed_sc + poly(Pioneer_Stabilizers,2) + Pioneer_Builders
                     + IP + WO + SL,
                       data = fall_df, family = binomial)
fallmod.o.step <- step(nullmod.fall, list(lower=formula(nullmod.fall), 
                                               upper=formula(fallmod.o)), direction = "both")
summary(fallmod.o.step) 

# exploring the top ranked models according to stepwise function 
fallmod.o.1 <- glmmTMB(Outcome ~ Object_Type + poly(Half,2) + Pioneer_Builders
                       + poly(Pioneer_Stabilizers,2) + IP, 
                         data = fall_df, family = binomial)
summary(fallmod.o.1) # 47.8

fallmod.o.2 <- glmmTMB(Outcome ~ Object_Type + poly(Half,2) + Pioneer_Builders
                       + poly(Pioneer_Stabilizers,2) + IP + WO, 
                         data = fall_df, family = binomial)
summary(fallmod.o.2) # 47.9

fallmod.o.3 <- glmmTMB(Outcome ~ Object_Type + poly(Half,2) + Pioneer_Builders
                       + poly(Pioneer_Stabilizers,2) + IP 
                       + poly(Secondary_Stabilizers,2),
                         data = fall_df, family = binomial)
summary(fallmod.o.3) # 48.2

fallmod.o.4 <- glmmTMB(Outcome ~ Object_Type + poly(Half,2) + Pioneer_Builders
                       + poly(Pioneer_Stabilizers,2) + IP + SL, 
                         data = fall_df, family = binomial)
summary(fallmod.o.4) # 48.4

fallmod.o.5 <- glmmTMB(Outcome ~ Object_Type + poly(Half,2) + Pioneer_Builders
                       + poly(Pioneer_Stabilizers,2) + IP + Viewshed_sc, 
                         data = fall_df, family = binomial)
summary(fallmod.o.5) # 49.1

fallmod.o.6 <- glmmTMB(Outcome ~ Object_Type + poly(Half,2) + Pioneer_Builders
                       + poly(Pioneer_Stabilizers,2), 
                         data = fall_df, family = binomial)
summary(fallmod.o.6) # 57.2

# the most parsimonious model w/in 2 AIC of top-ranked model is: summermod.o.1

# create list of models
models.fall.o <- list(
  fallmod.o.1 = fallmod.o.1,
  fallmod.o.2 = fallmod.o.2,
  fallmod.o.3 = fallmod.o.3,
  fallmod.o.4 = fallmod.o.4,
  fallmod.o.5 = fallmod.o.5,
  fallmod.o.6 = fallmod.o.6
)

# check all for goodness of fit
fit_output <- run_all_dharma_diagnostics(models.fall.o)
fit_output$summary
# all look ok here

# create candidate model table
cand.mods.fall.o <- create_candidate_model_table(models.fall.o)

# Selected top model; lowest AIC, most parsimonious
finalmod.fall.o <- fallmod.o.1
summary(finalmod.fall.o)

# Test adding random effect into selected model:
# point ID is nested within the pair ID 
# including PairID accounts for resampling the same territory. 
# can't use PointID because it perfectly predicts the outcome; nest is always nest

# selected model
REmod.fall.o <- glmmTMB(Outcome ~ (1|PairID) + Object_Type + poly(Half,2) 
                        + Pioneer_Builders + poly(Pioneer_Stabilizers,2) + IP,
                          family = binomial, data = fall_df)
summary(REmod.fall.o)
# adding PairID random effect does not change coefficient estimates or p-value, 
# explains no additional variance

# Random effect conclusions:
# The effect is negligible; does not explain the variation in outcome. 

# 4. Candidate model table ----
# Show all models within 2 AIC of top-ranked model, the next best models after
# that until the first one where support drops off
# include the null model

# ALL DATA
cand.mods.full.fg
cand.mods.full.o

# BREEDING SEASON
cand.mods.summer.fg
cand.mods.summer.o

# FALL
cand.mods.fall.fg
cand.mods.fall.o

# Combine them into one table
cand.mods.all <- dplyr::bind_rows(
  cand.mods.full.fg,
  cand.mods.full.o,
  cand.mods.summer.fg,
  cand.mods.summer.o,
  cand.mods.fall.fg,
  cand.mods.fall.o
)

write_csv(cand.mods.all, "tables/Candidate_Models_AllSubsets.csv")
