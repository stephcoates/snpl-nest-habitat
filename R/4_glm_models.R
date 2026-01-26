# Modeling nests vs. random sites with logistic regression
# 1. define the null model for all data, summer subset, and fall subset
# 2. test use of linear or quadratic terms for covariates
# 3. Check VIF, create global models and test candidate models for GoF and RE
# 3.a. all data: summarized microhabitat, functional groups, individual species
# 3.b. breeding season data: summarized microhabitat, functional groups, 
# individual species
# 3.c. fall data: summarized microhabitat, functional groups, individual species
# 4. Candidate model table (Tables S1)
# 5. Model results table (Tables S2)

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

# 2. Test linear and quadratic terms of each covariate ----

# Functional groups and summarized microhabitat covariates
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
# use linear for IP(222.0), WO (207.8), SL (236.0), DG, VE, BP, SG
print(aic_results_summer.overall)
# quadratic: SR, DG, SL
# linear: WO, IP, SL, BS, VE, BP, SG
print(aic_results_fall.overall)
# quadratic: BS
# linear: WO, IP, SR, SL


# 3. Global Models and stepwise selection to get candidate models ----
# 3a. ALL data: Summarized Microhabitat, Functional Groups, and Species ----
# linear vs. quadratic terms decided from the results of the function; add 
# random effect back in later. 

# Check vif for covariates in global model - SUMMARIZED MICROHABITAT
fullmod.lm <-lm(Outcome ~ Viewshed_sc + poly(Half,2) + poly(Two,2) + poly(TwoVeg,2)
                + NonVeg + poly(DeadVegManual,2),
                data = all_df)
vif(fullmod.lm)
# Two and TwoVeg have high VIF; don't use in same model

# Check vif for covariates in global model - FUNCTIONAL GROUPS
fullmod.lm <-lm(Outcome ~ Viewshed_sc + poly(Pioneer_Stabilizers,2)
                + NonVeg + Secondary_Stabilizers + Pioneer_Builders,
                data = all_df)
vif(fullmod.lm)
# All ok

# Check vif for covariates in global model - INDIVIDUAL SPECIES
fullmod.lm <-lm(Outcome ~ Viewshed_sc + poly(BS,2) + poly(SR,2)
                + WO + IP + SL,
                data = all_df)
vif(fullmod.lm)
# All ok

# SUMMARIZED MICROHABITAT - stepwise function to select candidate models
fullmod <- glmmTMB(Outcome ~ Object_Type + Viewshed_sc + poly(Half,2) 
                   + poly(TwoVeg,2)
                   + NonVeg + poly(DeadVegManual,2),
                   data = all_df, family = binomial)
fullmod.step <- step(nullmod.noR, list(lower=formula(nullmod.noR), 
                                       upper=formula(fullmod)), direction = "both")
summary(fullmod.step) 
#using TwoVeg: AIC = 118.5; Outcome ~ Object_Type + poly(TwoVeg, 2) + poly(Half, 2)
#using Two: AIC = 119.5; Outcome ~ Object_Type + poly(Two, 2) + poly(Half, 2)


# FUNCTIONAL GROUPS - stepwise function to select candidate models
fullmod <- glmmTMB(Outcome ~ Object_Type + Viewshed_sc + poly(Pioneer_Stabilizers,2) 
                   + Secondary_Stabilizers + NonVeg + Pioneer_Builders,
                   data = all_df, family = binomial)
fullmod.step <- step(nullmod.noR, list(lower=formula(nullmod.noR), 
                                       upper=formula(fullmod)), direction = "both")
summary(fullmod.step) 
#using PioneerStabs^2: AIC = 119.34; Outcome ~ Object_Type + poly(Pioneer_Stabilizers, 2) + Pioneer_Builders
#using PioneerStabs: AIC = 128.57; Outcome ~ Object_Type + Pioneer_Builders + Pioneer_Stabilizers

# INDIVIDUAL SPECIES - stepwise function to select candidate models
fullmod <- glmmTMB(Outcome ~ Object_Type + Viewshed_sc + poly(SR,2) + poly(BS,2)
                   + IP + WO + SL,
                   data = all_df, family = binomial)
fullmod.step <- step(nullmod.noR, list(lower=formula(nullmod.noR), 
                                       upper=formula(fullmod)), direction = "both")
summary(fullmod.step) 
# AIC = 119.67; Outcome ~ Object_Type + poly(BS, 2) + IP + poly(SR, 2) + WO

# 3.b. Candidate Models for All Data: Summarized Microhabitat, Functional Groups,
# and Individual Species

# SUMMARIZED MICROHABITAT - exploring the top ranked models according to stepwise function 
fullmod.1.sm <- glmmTMB(Outcome ~ Object_Type + poly(TwoVeg,2) +  poly(Half,2), 
                     data = all_df, family = binomial)
summary(fullmod.1.sm) #118.5

fullmod.2.sm <- glmmTMB(Outcome ~ Object_Type + poly(TwoVeg,2) +  poly(Half,2) 
                     + Viewshed_sc, 
                     data = all_df, family = binomial)
summary(fullmod.2.sm) # 118.9

fullmod.3.sm <- glmmTMB(Outcome ~ Object_Type + poly(TwoVeg,2) +  poly(Half,2)
                     + NonVeg, 
                     data = all_df, family = binomial)
summary(fullmod.3.sm) # 120.1

fullmod.4.sm <- glmmTMB(Outcome ~ Object_Type + poly(TwoVeg,2) +  poly(Half,2)
                     + poly(DeadVegManual,2), 
                     data = all_df, family = binomial)
summary(fullmod.4.sm) # 122.4


# create list of models
models.all.sm <- list(
  fullmod.1.sm = fullmod.1.sm,
  fullmod.2.sm = fullmod.2.sm,
  fullmod.3.sm = fullmod.3.sm,
  fullmod.4.sm = fullmod.4.sm
)

# FUNCTIONAL GROUPS - exploring the top ranked models according to stepwise function 
fullmod.1.fg <- glmmTMB(Outcome ~ Object_Type + poly(Pioneer_Stabilizers, 2) 
                     + Pioneer_Builders, 
                     data = all_df, family = binomial)
summary(fullmod.1.fg) # 119.34

fullmod.2.fg <- glmmTMB(Outcome ~ Object_Type + poly(Pioneer_Stabilizers, 2) 
                        + Pioneer_Builders + NonVeg, 
                     data = all_df, family = binomial)
summary(fullmod.2.fg) # 119.83

fullmod.3.fg <- glmmTMB(Outcome ~ Object_Type + poly(Pioneer_Stabilizers, 2) 
                        + Pioneer_Builders + Viewshed_sc, 
                     data = all_df, family = binomial)
summary(fullmod.3.fg) # 119.94

fullmod.4.fg <- glmmTMB(Outcome ~ Object_Type + poly(Pioneer_Stabilizers, 2) 
                        + Pioneer_Builders + Secondary_Stabilizers, 
                     data = all_df, family = binomial)
summary(fullmod.4.fg) # 121.34

fullmod.5.fg <- glmmTMB(Outcome ~ Object_Type + poly(Pioneer_Stabilizers, 2), 
                        data = all_df, family = binomial)
summary(fullmod.5.fg) # 133.74



# create list of models
models.all.fg <- list(
  fullmod.1.fg = fullmod.1.fg,
  fullmod.2.fg = fullmod.2.fg,
  fullmod.3.fg = fullmod.3.fg,
  fullmod.4.fg = fullmod.4.fg,
  fullmod.5.fg = fullmod.5.fg
)

# INDIVIDUAL SPECIES - exploring the top ranked models according to stepwise function 
fullmod.1.is <- glmmTMB(Outcome ~ Object_Type + poly(BS, 2) + IP + poly(SR, 2) 
                        + WO, 
                        data = all_df, family = binomial)
summary(fullmod.1.is) # 119.67

fullmod.2.is <- glmmTMB(Outcome ~ Object_Type + poly(BS, 2) + IP + poly(SR, 2),
                        data = all_df, family = binomial)
summary(fullmod.2.is) # 120.64

fullmod.3.is <- glmmTMB(Outcome ~ Object_Type + poly(BS, 2) + IP + poly(SR, 2) 
                        + WO + Viewshed_sc, 
                        data = all_df, family = binomial)
summary(fullmod.3.is) # 120.83

fullmod.4.is <- glmmTMB(Outcome ~ Object_Type + poly(BS, 2) + IP + poly(SR, 2) 
                        + WO + SL,  
                        data = all_df, family = binomial)
summary(fullmod.4.is) # 121.65

fullmod.5.is <- glmmTMB(Outcome ~ Object_Type + poly(BS, 2) + IP
                        + WO,  
                        data = all_df, family = binomial)
summary(fullmod.5.is) # 124.83


# create list of models
models.all.is <- list(
  fullmod.1.is = fullmod.1.is,
  fullmod.2.is = fullmod.2.is,
  fullmod.3.is = fullmod.3.is,
  fullmod.4.is = fullmod.4.is,
  fullmod.5.is = fullmod.5.is
)

# check all for goodness of fit
fit_output <- run_all_dharma_diagnostics(models.all.sm)
fit_output$summary
# all look ok here
fit_output <- run_all_dharma_diagnostics(models.all.fg)
fit_output$summary
# all ok, but p values range from 0.07 to 0.13
fit_output <- run_all_dharma_diagnostics(models.all.is)
fit_output$summary
# all ok

# create candidate model table
models.all.all <- list(
  fullmod.1.sm = fullmod.1.sm,
  fullmod.2.sm = fullmod.2.sm,
  fullmod.3.sm = fullmod.3.sm,
  fullmod.4.sm = fullmod.4.sm,
  fullmod.1.fg = fullmod.1.fg,
  fullmod.2.fg = fullmod.2.fg,
  fullmod.3.fg = fullmod.3.fg,
  fullmod.4.fg = fullmod.4.fg,
  fullmod.5.fg = fullmod.5.fg,
  fullmod.1.is = fullmod.1.is,
  fullmod.2.is = fullmod.2.is,
  fullmod.3.is = fullmod.3.is,
  fullmod.4.is = fullmod.4.is,
  fullmod.5.is = fullmod.5.is
)

cand.mods.full.sm <- create_candidate_model_table(models.all.sm)
cand.mods.full.fg <- create_candidate_model_table(models.all.fg)
cand.mods.full.is <- create_candidate_model_table(models.all.is)
cand.mods.full.all <- create_candidate_model_table(models.all.all)

# Selected top model; lowest AIC, most parsimonious
finalmod.all.sm <- fullmod.1.sm
finalmod.all.fg <- fullmod.1.fg
finalmod.all.is <- fullmod.2.is

# Test adding random effect into selected model:
# point ID is nested within the pair ID 
# including PairID accounts for resampling the same territory. 
# can't use PointID because it perfectly predicts the outcome; nest is always nest

# selected models
REmod.all.sm <- glmmTMB(Outcome ~ (1|PairID) + Object_Type + poly(TwoVeg,2)
                        + poly(Half,2),
                        family = binomial, data = all_df)
summary(REmod.all.sm)
# adding PairID random effect does not change coefficient estimates, explains no 
# additional variance (AIC = 120.5)

REmod.all.fg <- glmmTMB(Outcome ~ (1|PairID) + Object_Type
                        + poly(Pioneer_Stabilizers, 2) + Pioneer_Builders,
                        family = binomial, data = all_df)
summary(REmod.all.fg)
# adding PairID random effect does not change coefficient estimates, explains no 
# additional variance (AIC = 121.3)

REmod.all.is <- glmmTMB(Outcome ~ (1|PairID) + Object_Type
                        + poly(BS, 2) + IP + poly(SR, 2),
                        family = binomial, data = all_df)
summary(REmod.all.is)
# adding PairID random effect does not change coefficient estimates, explains no 
# additional variance (AIC = 122.6)

# Random effect conclusions:
# The effect is negligible; does not explain the variation in outcome. 
# Having a random effect also creates convergence problems because it explains 
# no variance.

# 3b. BREEDING subset: Summarized Microhabitat, Functional Groups, and Species ----
# linear vs. quadratic terms decided from the results of the function; add 
# random effect back in later. 

# Check vif for covariates in global model - SUMMARIZED MICROHABITAT
summermod.lm <-lm(Outcome ~ Viewshed_sc + poly(Half,2) + poly(Two,2) + poly(TwoVeg,2)
                + NonVeg + poly(DeadVegManual,2),
                data = summer_df)
vif(summermod.lm)
# Two and TwoVeg have high VIF, don't use in same model

# Check vif for covariates in global model - FUNCTIONAL GROUPS
summermod.lm <-lm(Outcome ~ Viewshed_sc + poly(Secondary_Stabilizers,2)
                + NonVeg + Pioneer_Stabilizers + Pioneer_Builders,
                data = summer_df)
vif(summermod.lm)
# all ok

# Check vif for covariates in global model - INDIVIDUAL SPECIES
summermod.lm <-lm(Outcome ~ Viewshed_sc + BS + poly(SR,2)
                + WO + IP + SL,
                data = summer_df)
vif(summermod.lm)
# All ok

# SUMMARIZED MICROHABITAT - stepwise function to select candidate models
summermod <- glmmTMB(Outcome ~ Object_Type + Viewshed_sc + poly(Half,2) 
                     + poly(TwoVeg,2)
                     + NonVeg + poly(DeadVegManual,2),
                   data = summer_df, family = binomial)
summermod.step <- step(nullmod.nesting, list(lower=formula(nullmod.nesting), 
                                       upper=formula(summermod)), direction = "both")
summary(summermod.step) 
#using TwoVeg: AIC = 53.78; Outcome ~ Object_Type + poly(TwoVeg, 2) + poly(Half, 2)
#using Two: AIC = 54.02; Outcome ~ Object_Type + poly(Two, 2) + poly(Half, 2)


# FUNCTIONAL GROUPS - stepwise function to select candidate models
summermod <- glmmTMB(Outcome ~ Object_Type + Viewshed_sc + Secondary_Stabilizers
                   + NonVeg + Pioneer_Stabilizers + Pioneer_Builders,
                   data = summer_df, family = binomial)
summermod.step <- step(nullmod.nesting, list(lower=formula(nullmod.nesting), 
                                       upper=formula(summermod)), direction = "both")
summary(summermod.step) 
#using SecondStabs^2: AIC = 62.42; Outcome ~ Object_Type + poly(Secondary_Stabilizers, 2) + Pioneer_Builders
#using SecondStabs: AIC = 62.28; Outcome ~ Object_Type + Pioneer_Stabilizers + Pioneer_Builders

# INDIVIDUAL SPECIES - stepwise function to select candidate models
summermod <- glmmTMB(Outcome ~ Object_Type + Viewshed_sc + poly(SR,2) + BS
                   + IP + WO + SL,
                   data = summer_df, family = binomial)
summermod.step <- step(nullmod.nesting, list(lower=formula(nullmod.nesting), 
                                       upper=formula(summermod)), direction = "both")
summary(summermod.step) 
# AIC = 61.89; Outcome ~ Object_Type + BS + poly(SR, 2)

# 3.b. Candidate Models for SUMMER Data: Summarized Microhabitat, Functional Groups,
# and Individual Species

# SUMMARIZED MICROHABITAT - exploring the top ranked models according to stepwise function 
summermod.1.sm <- glmmTMB(Outcome ~ Object_Type + poly(TwoVeg, 2) + poly(Half, 2), 
                        data = summer_df, family = binomial)
summary(summermod.1.sm) #53.78

summermod.2.sm <- glmmTMB(Outcome ~ Object_Type + poly(TwoVeg, 2) + poly(Half, 2)
                          + NonVeg, 
                        data = summer_df, family = binomial)
summary(summermod.2.sm) # 55.291

summermod.3.sm <- glmmTMB(Outcome ~ Object_Type + poly(TwoVeg,2) +  poly(Half,2)
                        + Viewshed_sc, 
                        data = summer_df, family = binomial)
summary(summermod.3.sm) # 55.774

summermod.4.sm <- glmmTMB(Outcome ~ Object_Type + poly(TwoVeg,2) +  poly(Half,2)
                        + poly(DeadVegManual,2), 
                        data = summer_df, family = binomial)
summary(summermod.4.sm) # 56.173


# create list of models
models.summer.sm <- list(
  summermod.1.sm = summermod.1.sm,
  summermod.2.sm = summermod.2.sm,
  summermod.3.sm = summermod.3.sm,
  summermod.4.sm = summermod.4.sm
)

# FUNCTIONAL GROUPS - exploring the top ranked models according to stepwise function 
summermod.1.fg <- glmmTMB(Outcome ~ Object_Type + Pioneer_Stabilizers 
                        + Pioneer_Builders, 
                        data = summer_df, family = binomial)
summary(summermod.1.fg) # 64.28

summermod.2.fg <- glmmTMB(Outcome ~ Object_Type + Pioneer_Stabilizers, 
                        data = summer_df, family = binomial)
summary(summermod.2.fg) # 64.417

summermod.3.fg <- glmmTMB(Outcome ~ Object_Type + Pioneer_Stabilizers 
                        + Pioneer_Builders + NonVeg, 
                        data = summer_df, family = binomial)
summary(summermod.3.fg) # 64.742

summermod.4.fg <- glmmTMB(Outcome ~ Object_Type + Pioneer_Stabilizers 
                        + Pioneer_Builders  + Secondary_Stabilizers, 
                        data = summer_df, family = binomial)
summary(summermod.4.fg) # 65.718

summermod.5.fg <- glmmTMB(Outcome ~ Object_Type + Pioneer_Stabilizers
                        + Pioneer_Builders + Viewshed_sc,
                        data = summer_df, family = binomial)
summary(summermod.5.fg) # 66.269

summermod.6.fg <- glmmTMB(Outcome ~ Object_Type + Pioneer_Builders,
                          data = summer_df, family = binomial)
summary(summermod.6.fg) # 68.187


# create list of models
models.summer.fg <- list(
  summermod.1.fg = summermod.1.fg,
  summermod.2.fg = summermod.2.fg,
  summermod.3.fg = summermod.3.fg,
  summermod.4.fg = summermod.4.fg,
  summermod.5.fg = summermod.5.fg,
  summermod.6.fg = summermod.6.fg
)

# INDIVIDUAL SPECIES - exploring the top ranked models according to stepwise function 
summermod.1.is <- glmmTMB(Outcome ~ Object_Type + BS + poly(SR, 2), 
                        data = summer_df, family = binomial)
summary(summermod.1.is) # 61.89

summermod.2.is <- glmmTMB(Outcome ~ Object_Type + BS + poly(SR, 2) + WO,
                        data = summer_df, family = binomial)
summary(summermod.2.is) # 63.070

summermod.3.is <- glmmTMB(Outcome ~ Object_Type + BS + poly(SR, 2) + SL, 
                        data = summer_df, family = binomial)
summary(summermod.3.is) # 63.245

summermod.4.is <- glmmTMB(Outcome ~ Object_Type + BS + poly(SR, 2) + IP,  
                        data = summer_df, family = binomial)
summary(summermod.4.is) # 63.67

summermod.5.is <- glmmTMB(Outcome ~ Object_Type + BS + poly(SR, 2) + Viewshed_sc,  
                        data = summer_df, family = binomial)
summary(summermod.5.is) # 63.871

summermod.6.is <- glmmTMB(Outcome ~ Object_Type + BS + Viewshed_sc,  
                          data = summer_df, family = binomial)
summary(summermod.6.is) # 64.417


# create list of models
models.summer.is <- list(
  summermod.1.is = summermod.1.is,
  summermod.2.is = summermod.2.is,
  summermod.3.is = summermod.3.is,
  summermod.4.is = summermod.4.is,
  summermod.5.is = summermod.5.is,
  summermod.6.is = summermod.6.is
)

# check all for goodness of fit
fit_output <- run_all_dharma_diagnostics(models.summer.sm)
fit_output$summary
# all look ok here
fit_output <- run_all_dharma_diagnostics(models.summer.fg)
fit_output$summary
# all ok
fit_output <- run_all_dharma_diagnostics(models.summer.is)
fit_output$summary
# all ok

# create overall candidate model table for summer subset
models.summer.all <- list(
  summermod.1.sm = summermod.1.sm,
  summermod.2.sm = summermod.2.sm,
  summermod.3.sm = summermod.3.sm,
  summermod.4.sm = summermod.4.sm,
  summermod.1.fg = summermod.1.fg,
  summermod.2.fg = summermod.2.fg,
  summermod.3.fg = summermod.3.fg,
  summermod.4.fg = summermod.4.fg,
  summermod.5.fg = summermod.5.fg,
  summermod.6.fg = summermod.6.fg,
  summermod.1.is = summermod.1.is,
  summermod.2.is = summermod.2.is,
  summermod.3.is = summermod.3.is,
  summermod.4.is = summermod.4.is,
  summermod.5.is = summermod.5.is,
  summermod.6.is = summermod.6.is
)

cand.mods.summer.sm <- create_candidate_model_table(models.summer.sm)
cand.mods.summer.fg <- create_candidate_model_table(models.summer.fg)
cand.mods.summer.is <- create_candidate_model_table(models.summer.is)
cand.mods.summer.all <- create_candidate_model_table(models.summer.all)

# Selected top model; lowest AIC, most parsimonious
finalmod.summer.sm <- summermod.1.sm
finalmod.summer.fg <- summermod.2.fg
finalmod.summer.is <- summermod.1.is

# Test adding random effect into selected model:
# point ID is nested within the pair ID 
# including PairID accounts for resampling the same territory. 
# can't use PointID because it perfectly predicts the outcome; nest is always nest

# selected models
REmod.summer.sm <- glmmTMB(Outcome ~ (1|PairID) + Object_Type + poly(TwoVeg,2)
                        + poly(Half,2),
                        family = binomial, data = summer_df)
summary(REmod.summer.sm)
# adding PairID random effect does not change coefficient estimates, explains no 
# additional variance (AIC = 55.8)

REmod.summer.fg <- glmmTMB(Outcome ~ (1|PairID) + Object_Type 
                           + Pioneer_Stabilizers,
                        family = binomial, data = summer_df)
summary(REmod.summer.fg)
# adding PairID random effect does not change coefficient estimates, explains no 
# additional variance (AIC = 66.4)

REmod.summer.is <- glmmTMB(Outcome ~ (1|PairID) + Object_Type
                        + BS + poly(SR, 2),
                        family = binomial, data = summer_df)
summary(REmod.summer.is)
# adding PairID random effect does not change coefficient estimates, explains no 
# additional variance (AIC = 63.9)

# Random effect conclusions:
# The effect is negligible; does not explain the variation in outcome. 
# Having a random effect also creates convergence problems because it explains 
# no variance.

# 3c. FALL subset: Summarized Microhabitat, Functional Groups, and Species ----
# linear vs. quadratic terms decided from the results of the function; add 
# random effect back in later. 

print(aic_results_fall)
# Quadratic: Half, Two, TwoVeg, Secondary Stabilizers, DeadVegManual, Pioneer Stabilizers
# Linear: NonVeg, Viewshed_sc
# Cusp: Pioneer Builder
print(aic_results_fall.overall)
# quadratic: BS
# linear: WO, IP, SR, SL


# Check vif for covariates in global model - SUMMARIZED MICROHABITAT
fallmod.lm <-lm(Outcome ~ Viewshed_sc + poly(Half,2) + poly(Two,2) + poly(TwoVeg,2)
                  + NonVeg + poly(DeadVegManual,2),
                  data = fall_df)
vif(fallmod.lm)
# Two and TwoVeg have high VIF, don't use in same model

# Check vif for covariates in global model - FUNCTIONAL GROUPS
fallmod.lm <-lm(Outcome ~ Viewshed_sc + poly(Secondary_Stabilizers,2)
                  + NonVeg + poly(Pioneer_Stabilizers,2) + Pioneer_Builders,
                  data = fall_df)
vif(fallmod.lm)
# all ok

# Check vif for covariates in global model - INDIVIDUAL SPECIES
fallmod.lm <-lm(Outcome ~ Viewshed_sc + SR + poly(BS,2)
                  + WO + IP + SL,
                  data = fall_df)
vif(fallmod.lm)
# All ok

# SUMMARIZED MICROHABITAT - stepwise function to select candidate models
fallmod <- glmmTMB(Outcome ~ Object_Type + Viewshed_sc + poly(Half,2) 
                   + poly(TwoVeg,2)
                   + NonVeg + poly(DeadVegManual,2),
                     data = fall_df, family = binomial)
fallmod.step <- step(nullmod.fall, list(lower=formula(nullmod.fall), 
                                             upper=formula(fallmod)), direction = "both")
summary(fallmod.step) 
#using TwoVeg: AIC = 70.33; Object_Type + poly(TwoVeg, 2) + Viewshed_sc + poly(Half, 2)
#using Two: AIC = 70.61; Object_Type + poly(Two, 2) + Viewshed_sc + poly(Half, 2)


# FUNCTIONAL GROUPS - stepwise function to select candidate models
fallmod <- glmmTMB(Outcome ~ Object_Type + Viewshed_sc + poly(Secondary_Stabilizers,2)
                   + NonVeg + poly(Pioneer_Stabilizers,2) + Pioneer_Builders,
                     data = fall_df, family = binomial)
fallmod.step <- step(nullmod.fall, list(lower=formula(nullmod.fall), 
                                             upper=formula(fallmod)), direction = "both")
summary(fallmod.step) 
# AIC = 64.06; Outcome ~ Object_Type + poly(Pioneer_Stabilizers, 2) 
# + Pioneer_Builders + Viewshed_sc

# INDIVIDUAL SPECIES - stepwise function to select candidate models
fallmod <- glmmTMB(Outcome ~ Object_Type + Viewshed_sc + SR + poly(BS,2)
                   + WO + IP + SL,
                     data = fall_df, family = binomial)
fallmod.step <- step(nullmod.fall, list(lower=formula(nullmod.fall), 
                                             upper=formula(fallmod)), direction = "both")
summary(fallmod.step) 
# AIC = 61.89; Outcome ~ Object_Type + BS + poly(SR, 2)

# Candidate Models for FALL Data: Summarized Microhabitat, Functional Groups,
# and Individual Species

# SUMMARIZED MICROHABITAT - exploring the top ranked models according to stepwise function 
fallmod.1.sm <- glmmTMB(Outcome ~ Object_Type + poly(TwoVeg, 2) + poly(Half, 2)
                          + Viewshed_sc, 
                          data = fall_df, family = binomial)
summary(fallmod.1.sm) # 70.33

fallmod.2.sm <- glmmTMB(Outcome ~ Object_Type + poly(TwoVeg, 2) 
                        + Viewshed_sc, 
                          data = fall_df, family = binomial)
summary(fallmod.2.sm) # 70.53

fallmod.3.sm <- glmmTMB(Outcome ~ Object_Type + poly(TwoVeg,2) + poly(Half,2), 
                          data = fall_df, family = binomial)
summary(fallmod.3.sm) # 70.75

fallmod.4.sm <- glmmTMB(Outcome ~ Object_Type + poly(TwoVeg, 2) + poly(Half, 2)
                        + Viewshed_sc + NonVeg, 
                          data = fall_df, family = binomial)
summary(fallmod.4.sm) # 71.880

fallmod.5.sm <- glmmTMB(Outcome ~ Object_Type + poly(TwoVeg, 2) + poly(Half, 2)
                        + Viewshed_sc + poly(DeadVegManual,2), 
                        data = fall_df, family = binomial)
summary(fallmod.5.sm) # 72.545


# create list of models
models.fall.sm <- list(
  fallmod.1.sm = fallmod.1.sm,
  fallmod.2.sm = fallmod.2.sm,
  fallmod.3.sm = fallmod.3.sm,
  fallmod.4.sm = fallmod.4.sm,
  fallmod.5.sm = fallmod.5.sm
)

# FUNCTIONAL GROUPS - exploring the top ranked models according to stepwise function 
fallmod.1.fg <- glmmTMB(Outcome ~ Object_Type + poly(Pioneer_Stabilizers, 2) 
                          + Pioneer_Builders + Viewshed_sc, 
                          data = fall_df, family = binomial)
summary(fallmod.1.fg) # 64.06

fallmod.2.fg <- glmmTMB(Outcome ~ Object_Type + poly(Pioneer_Stabilizers, 2) 
                        + Pioneer_Builders + Viewshed_sc + NonVeg, 
                          data = fall_df, family = binomial)
summary(fallmod.2.fg) # 65.655

fallmod.3.fg <- glmmTMB(Outcome ~ Object_Type + poly(Pioneer_Stabilizers, 2) 
                        + Pioneer_Builders, 
                          data = fall_df, family = binomial)
summary(fallmod.3.fg) # 65.728

fallmod.4.fg <- glmmTMB(Outcome ~ Object_Type + poly(Pioneer_Stabilizers, 2) 
                        + Pioneer_Builders + Viewshed_sc + poly(Secondary_Stabilizers,2), 
                          data = fall_df, family = binomial)
summary(fallmod.4.fg) # 67.398


# create list of models
models.fall.fg <- list(
  fallmod.1.fg = fallmod.1.fg,
  fallmod.2.fg = fallmod.2.fg,
  fallmod.3.fg = fallmod.3.fg,
  fallmod.4.fg = fallmod.4.fg
)

# INDIVIDUAL SPECIES - exploring the top ranked models according to stepwise function 
fallmod.1.is <- glmmTMB(Outcome ~ Object_Type + IP + poly(BS, 2) + SR 
                        + Viewshed_sc, 
                          data = fall_df, family = binomial)
summary(fallmod.1.is) # 63.219

fallmod.2.is <- glmmTMB(Outcome ~ Object_Type + IP + poly(BS, 2) + SR,
                          data = fall_df, family = binomial)
summary(fallmod.2.is) # 63.456

fallmod.3.is <- glmmTMB(Outcome ~ Object_Type + IP + poly(BS, 2) + Viewshed_sc, 
                          data = fall_df, family = binomial)
summary(fallmod.3.is) # 64.98

fallmod.4.is <- glmmTMB(Outcome ~ Object_Type + IP + poly(BS, 2) + SR 
                        + Viewshed_sc + WO,  
                          data = fall_df, family = binomial)
summary(fallmod.4.is) # 65.021

fallmod.5.is <- glmmTMB(Outcome ~ Object_Type + IP + poly(BS, 2) + SR 
                        + Viewshed_sc + SL,  
                          data = fall_df, family = binomial)
summary(fallmod.5.is) # 65.210

fallmod.6.is <- glmmTMB(Outcome ~ Object_Type + poly(BS, 2) + SR + Viewshed_sc,  
                        data = fall_df, family = binomial)
summary(fallmod.6.is) # 68.808


# create list of models
models.fall.is <- list(
  fallmod.1.is = fallmod.1.is,
  fallmod.2.is = fallmod.2.is,
  fallmod.3.is = fallmod.3.is,
  fallmod.4.is = fallmod.4.is,
  fallmod.5.is = fallmod.5.is,
  fallmod.6.is = fallmod.6.is
)

# check all for goodness of fit
fit_output <- run_all_dharma_diagnostics(models.fall.sm)
fit_output$summary
# all look ok here
fit_output <- run_all_dharma_diagnostics(models.fall.fg)
fit_output$summary
# all ok
fit_output <- run_all_dharma_diagnostics(models.fall.is)
fit_output$summary
# all ok

# create candidate model table

models.fall.all <- list(
  fallmod.1.sm = fallmod.1.sm,
  fallmod.2.sm = fallmod.2.sm,
  fallmod.3.sm = fallmod.3.sm,
  fallmod.4.sm = fallmod.4.sm,
  fallmod.5.sm = fallmod.5.sm,
  fallmod.1.fg = fallmod.1.fg,
  fallmod.2.fg = fallmod.2.fg,
  fallmod.3.fg = fallmod.3.fg,
  fallmod.4.fg = fallmod.4.fg,
  fallmod.1.is = fallmod.1.is,
  fallmod.2.is = fallmod.2.is,
  fallmod.3.is = fallmod.3.is,
  fallmod.4.is = fallmod.4.is,
  fallmod.5.is = fallmod.5.is,
  fallmod.6.is = fallmod.6.is
)
  
cand.mods.fall.sm <- create_candidate_model_table(models.fall.sm)
cand.mods.fall.fg <- create_candidate_model_table(models.fall.fg)
cand.mods.fall.is <- create_candidate_model_table(models.fall.is)
cand.mods.fall.all <- create_candidate_model_table(models.fall.all)

# Selected top model; lowest AIC, most parsimonious
finalmod.fall.sm.1 <- fallmod.2.sm
finalmod.fall.sm.2 <- fallmod.3.sm
finalmod.fall.fg <- fallmod.3.fg
finalmod.fall.is.1 <- fallmod.2.is
finalmod.fall.is.2 <- fallmod.3.is

# Test adding random effect into selected model:
# point ID is nested within the pair ID 
# including PairID accounts for resampling the same territory. 
# can't use PointID because it perfectly predicts the outcome; nest is always nest

# selected models
REmod.fall.sm <- glmmTMB(Outcome ~ (1|PairID) + Object_Type + poly(TwoVeg,2)
                           + Viewshed,
                           family = binomial, data = fall_df)
summary(REmod.fall.sm)
# adding PairID random effect does not change coefficient estimates, explains no 
# additional variance (AIC = 72.5)

REmod.fall.sm <- glmmTMB(Outcome ~ (1|PairID) + Object_Type + poly(TwoVeg,2)
                         + poly(Half,2),
                         family = binomial, data = fall_df)
summary(REmod.fall.sm)
# adding PairID random effect does not change coefficient estimates, explains no 
# additional variance (AIC = 72.8)

REmod.fall.fg <- glmmTMB(Outcome ~ (1|PairID) + Object_Type 
                         + poly(Pioneer_Stabilizers, 2) + Pioneer_Builders,
                           family = binomial, data =fall_df)
summary(REmod.fall.fg)
# adding PairID random effect does not change coefficient estimates, explains no 
# additional variance (AIC = 67.7)

REmod.fall.is <- glmmTMB(Outcome ~ (1|PairID) + Object_Type + IP + poly(BS, 2) + SR,
                           family = binomial, data = fall_df)
summary(REmod.fall.is)
# adding PairID random effect does not change coefficient estimates, explains no 
# additional variance (AIC = 65.5)

REmod.fall.is <- glmmTMB(Outcome ~ (1|PairID) + Object_Type + IP + poly(BS, 2) 
                         + Viewshed_sc,
                         family = binomial, data = fall_df)
summary(REmod.fall.is)
# adding PairID random effect does not change coefficient estimates, explains no 
# additional variance (AIC = 67.0)

# Random effect conclusions:
# The effect is negligible; does not explain the variation in outcome. 
# Having a random effect also creates convergence problems because it explains 
# no variance.


# 4. Candidate model table ----
# Show all models within 2 AIC of top-ranked model, the next best models after
# that until the first one where support drops off
# include the null model

# ALL DATA
cand.mods.full.sm
cand.mods.full.fg
cand.mods.full.is

# BREEDING SEASON
cand.mods.summer.sm
cand.mods.summer.fg
cand.mods.summer.is

# FALL
cand.mods.fall.sm
cand.mods.fall.fg
cand.mods.fall.is

# Combine them into one table
cand.mods.all.3groups <- dplyr::bind_rows(
  cand.mods.full.all,
  cand.mods.summer.all,
  cand.mods.fall.all
)

write_csv(cand.mods.all.3groups, "tables/Candidate_Models_AllSubsets_3Groups.csv")

# 5. Create model results tables ----
# results for all "selected" models 

#create list of selected models - all data
selectmods.alldata <- list(
  fullmod.1.sm = fullmod.1.sm,
  fullmod.1.fg = fullmod.1.fg,
  fullmod.2.is = fullmod.2.is
)

modresults.alldata <- tidy_mod_table(selectmods.alldata)
write_csv(modresults.alldata, "tables/ModelResultsTable_AllData.csv")

#create list of selected models - summer data
selectmods.summerdata <- list(
  summermod.1.sm = summermod.1.sm,
  summermod.2.fg = summermod.2.fg,
  summermod.1.is = summermod.1.is
)

modresults.summerdata <- tidy_mod_table(selectmods.summerdata)
write_csv(modresults.summerdata, "tables/ModelResultsTable_SummerData.csv")

#create list of selected models - fall data
selectmods.falldata <- list(
  fallmod.2.sm = fallmod.2.sm,
  fallmod.1.fg = fallmod.1.fg,
  fallmod.3.fg = fallmod.3.fg,
  fallmod.2.is = fallmod.2.is,
  fallmod.3.is = fallmod.3.is
)

modresults.falldata <- tidy_mod_table(selectmods.falldata)
write_csv(modresults.falldata, "tables/ModelResultsTable_FallData.csv")

# Interpretation help: get outputs for making odds ratio comparisons between object types ----
# Get pairwise contrasts on the log-odds scale
emm <- emmeans(fullmod.1.sm, ~Object_Type)
pw <- contrast(
  emm,
  method = "pairwise",
  adjust = "none"
)

# Convert emmeans pairwise contrasts to a summary table with OR, 95% CI, and p-values
# Turn into a tidy table and compute OR + 95% CI. 
tab <- as.data.frame(pw) %>%
  mutate(
    # Odds ratio and CI from log-odds (estimate) and SE
    OR      = exp(estimate),
    CI_low  = exp(estimate - 1.96 * SE),
    CI_high = exp(estimate + 1.96 * SE),
    
    # Clean p-value formatting 
    p_value = ifelse(p.value < 1e-4, "<0.0001", sprintf("%.4f", p.value)),
    
    # Make a nicer comparison label
    Comparison = contrast %>%
      str_replace_all("\\(", "") %>%
      str_replace_all("\\)", "") %>%
      str_replace_all("\\s+", " ") %>%
      str_replace(" - ", " vs ")
  ) %>%
  select(Comparison, OR, CI_low, CI_high, p_value) %>%
  mutate(
    OR = round(OR, 2),
    CI_low = round(CI_low, 2),
    CI_high = round(CI_high, 2),
    `95% CI` = paste0(CI_low, "–", CI_high)
  ) %>%
  select(Comparison, OR, `95% CI`, p_value)

# make final table
tab_final <- bind_rows(
  tab %>% filter(Comparison != "None vs Vegetation"),
  tab %>% filter(Comparison == "None vs Vegetation") %>%
    mutate(
      Comparison = "Vegetation vs None",
      OR = round(1 / OR, 2),
      # Recompute CI by swapping/inverting the bounds:
      `95% CI` = {
        # pull original CI numbers back out
        parts <- str_split_fixed(`95% CI`, "–", 2)
        lo <- as.numeric(parts[,1]); hi <- as.numeric(parts[,2])
        paste0(round(1 / hi, 2), "–", round(1 / lo, 2))
      }
    ) %>%
    select(Comparison, OR, `95% CI`, p_value)
) %>%
  mutate(Comparison = factor(
    Comparison,
    levels = c(
      "Non-veg Wood or Shell vs Vegetation",
      "Non-veg Wood or Shell vs None",
      "Vegetation vs None"
    )
  )) %>%
  arrange(Comparison) %>%
  mutate(Comparison = as.character(Comparison))

tab_final
