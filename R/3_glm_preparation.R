# Prepare data and examine covariates for logistic regression

# load packages and functions
source('R/packages.R')

alldata <- read_csv('data/SNPL_nest_habitat.csv') %>% 
  filter(!is.na(ID)) 
pairingsdata <- read_csv('data/nest_random_pairs.csv') %>% 
  filter(!is.na(`Collected Fall`)) 

# Process data ----
# add the nest-random point pairings
alldata1 <- alldata %>%
  left_join(pairingsdata, by = c("ID" = "rSX RP")) %>%
  mutate(PairID = case_when(
    Type == "Nest" ~ paste("sx", ID),
    Type == "Random Point" ~ paste("sx", `SX Nest`),
    TRUE ~ NA_character_
  )) %>%
  select(-`SX Nest`, -`Collected Nesting`, -`Collected Fall`) 

# data prep for analysis: convert to percentages, create object column, 
# add object type, scale viewshed, remove points that were in the water
reg1data <- alldata1 %>%
  mutate(across(Half:MI, ~ . * 100),  # convert to percentages
         across(Pioneer_Builders:NonVeg, ~ . *100), # convert to percentages
         Object = ifelse(grepl("None", Spec, ignore.case = TRUE), "Absent", "Present"), 
         Object_Type = case_when(  # Add Object Type column based on conditions
           Spec %in% c("Woody", "Shell") ~ "Non-veg Wood or Shell",
           Spec == "None" ~ "None",
           !(Spec %in% c("Woody", "Shell", "None")) ~ "Vegetation"),
         Viewshed_sc = Viewshed/360, 
         SeasonN = ifelse(Season=="Nesting",1,2),
         Outcome = ifelse(Type=="Nest", 1, 0) 
  ) %>%
  mutate(
    PointID = paste0(PairID, Type) 
  ) %>%
  filter(Spec != "Water") %>% 
  select(-Dead, -Unspec, -Spec, -`Collection Period`, -Fate, - DeadVeg)  

# create seasonal data sets: breeding season and post-breeding fall; filtering to 
# just the points that were sampled in both seasons.
reg1data_nesting <- reg1data %>%
  add_count(PairID, name = "pair_count") %>%  
  filter(pair_count == 4) %>%
  dplyr::select(-pair_count) %>%  
  filter(Season == "Nesting") 

reg1data_fall <- reg1data %>%
  add_count(PairID, name = "pair_count") %>%  
  filter(pair_count == 4) %>% 
  dplyr::select(-pair_count) %>%  
  filter(Season == "Fall")

# Assess covariates ----

# Remove cover types that occurred very rarely; 
# summarize how often there is 0% cover of a cover type across all samples
covertypes_df1 <- reg1data %>% dplyr::select(Half:MI, Viewshed_sc, Viewshed:NonVeg, Outcome) 
covertypes_df1.nest <- reg1data_nesting %>% dplyr::select(Half:MI, Viewshed_sc, Viewshed:NonVeg, Outcome) 
covertypes_df1.fall <- reg1data_fall %>% dplyr::select(Half:MI, Viewshed_sc, Viewshed:NonVeg, Outcome) 

# all data
summary_of_zeros1 <- covertypes_df1 %>%
  group_by(Outcome) %>%
  summarize(across(everything(), ~mean(. == 0) * 100))

# breeding subset
summary_of_zeros1.nest <- covertypes_df1.nest %>%
  group_by(Outcome) %>%
  summarize(across(everything(), ~mean(. == 0) * 100)) #%>%

# fall subset
summary_of_zeros1.fall <- covertypes_df1.fall %>%
  group_by(Outcome) %>%
  summarize(across(everything(), ~mean(. == 0) * 100)) #%>%

# eliminate variables that do not occur at >90% of nests and random sites
covertypes_remove <- summary_of_zeros1 %>%
  dplyr::select(where(~ all(. > 85, na.rm = TRUE))) %>%
  names()
covertypes_remove.nest <- summary_of_zeros1.nest %>%
  dplyr::select(where(~ all(. > 85, na.rm = TRUE))) %>%
  names()
covertypes_remove.fall <- summary_of_zeros1.fall %>%
  dplyr::select(where(~ all(. > 85, na.rm = TRUE))) %>%
  names()

# remove variables that we aren't interested in for this analysis
all_df <- reg1data %>% dplyr::select(-c(all_of(covertypes_remove), MI, Sand, Live))
summer_df <- reg1data_nesting %>% dplyr::select(-c(all_of(covertypes_remove), MI, Sand, Live))
fall_df <- reg1data_fall %>% dplyr::select(-c(all_of(covertypes_remove), MI, Sand, Live))

write.csv(all_df, "data/glm_alldat.csv")
write.csv(summer_df, "data/glm_breedingdat.csv")
write.csv(fall_df, "data/glm_falldat.csv")

# Test covariates for correlation ----
# covariates with corellation coeffecient >0.7 not included in same model

# all data
hab_corr_df1 <- all_df %>% dplyr::select(Half:SL, Viewshed_sc, Pioneer_Builders:NonVeg)
corr_matrix1 <- cor(hab_corr_df1)
# based on correlation coefficient >0.7, these can not be included in same model:
# Two and BS, Secondary_Stabilizers, Pioneer_Stabilizers
# SR and Pioneer_Builders 
# BS and Two, Pioneer_Stabilizer
# WO and NonVeg
# Pioneer_Builder and SR
# Secondary_Stabilizer and Two
# Pioneer_Stabilizer and Two, BS
# NonVeg and WO

# breeding subset
hab_corr_df1.nest <- summer_df %>% dplyr::select(Half:SL, Viewshed_sc, Pioneer_Builders:NonVeg)
corr_matrix1.nest <- cor(hab_corr_df1.nest)
# based on correlation coefficient >0.7, these can not be included in same model:
# Half and Two, BS, VE, Secondary Stabilizers, Pioneer Stabilizers
# Two and Half, BS, VE, Secondary Stabilizers, Pioneer Stabilizers
# SR and Pioneer Builders
# BS and Half, Two
# WO and non veg
# Pioneer Builders and SR
# Secondary Stabilizers and Half, Two
# Pioneer Stabilizers and Half, Two
# NonVeg and WO

# fall subset
hab_corr_df1.fall <- fall_df %>% dplyr::select(Half:SL, Viewshed_sc, Pioneer_Builders:NonVeg)
corr_matrix1.fall <- cor(hab_corr_df1.fall)
# based on correlation coefficient >0.7, these can not be included in same model:
# Two and DeadVegManual, BS, Secondary Stabilizers, Pioneer Stabilizers
# DeadVegManual and Two, BS, VE, Pioneer Stabilizers
# SR and Pioneer builders
# BS and Two, DeadVegManual
# SG and Secondary Stabilizers
# WO and non veg
# Pioneer Builders and SR
# Secondary Stabilizers and Two, SG
# Pioneer Stabilizers and Two, DeadVegManual
# NonVeg and WO