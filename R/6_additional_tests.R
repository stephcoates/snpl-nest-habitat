## Follow up tests based on glm results 
# looking into viewshed and iceplant: why are there differences for these 
# in breeding season versus fall models?

# load packages and functions and data
source('R/packages.R')
source('R/functions_glm.R')

# make an exploratory dataframe
ex_df <- read.csv("data/glm_alldat.csv")
  
# Compare ice plant cover at nests and random points in breeding season and fall
IP_tab <- ex_df %>%
  group_by(Season, Type) %>%
  summarise(
    mean_IP = mean(IP, na.rm = TRUE),
    min_IP = min(IP, na.rm = TRUE),
    max_IP = max(IP, na.rm = TRUE),
    n = sum(!is.na(IP)),
    .groups = "drop"
  )
# test for statistical differences:
ex_df$Season <- factor(ex_df$Season)
ex_df$Type   <- factor(ex_df$Type)

anova_mod <- aov(IP ~ Season * Type, data = ex_df)
summary(anova_mod)

emmeans(anova_mod, pairwise ~ Season * Type)
# these tests tell us that nests and random points have different IP cover (specifically
# nests and random points within both seasons), but IP cover is not statistically
# different between seasons

# Compare viewshed at nests and random points in breeding season and fall
View_tab <- ex_df %>%
  group_by(Season, Type) %>%
  summarise(
    mean_View = mean(Viewshed, na.rm = TRUE),
    min_View = min(Viewshed, na.rm = TRUE),
    max_View = max(Viewshed, na.rm = TRUE),
    n = sum(!is.na(Viewshed)),
    .groups = "drop"
  )

# test for statistical difference
anova_mod_v <- aov(Viewshed ~ Season * Type, data = ex_df)
summary(anova_mod_v)

emmeans(anova_mod_v, pairwise ~ Season * Type)
# same thing here, type (nest or random point) is statistically different, this 
# time only for Breeding season nests and fall random points, but Viewshed is not
# statistically different by season.

## Compare changes for other variables in the fall vs. summer models ----
# IP and Viewshed showing up in fall models but not summer models isn't really 
# explained by seasonal differences in either (because it's not statistically different)

# Compare TwoVeg at nests and random points in breeding season and fall
Two_tab <- ex_df %>%
  group_by(Season, Type) %>%
  summarise(
    mean_TwoVeg = mean(TwoVeg, na.rm = TRUE),
    min_TwoVeg = min(TwoVeg, na.rm = TRUE),
    max_TwoVeg = max(TwoVeg, na.rm = TRUE),
    n = sum(!is.na(TwoVeg)),
    .groups = "drop"
  )
# table tells us % cover of veg in 2m is higher on average at random points than nests.

# test for statistical difference
anova_mod_2 <- aov(TwoVeg ~ Season * Type, data = ex_df)
summary(anova_mod_2)

emmeans(anova_mod_2, pairwise ~ Season * Type)
# again it's by type (nest vs. random); no statistical significance with season

# Compare Half at nests and random points in breeding season and fall
Half_tab <- ex_df %>%
  group_by(Season, Type) %>%
  summarise(
    mean_Half = mean(Half, na.rm = TRUE),
    min_Half = min(Half, na.rm = TRUE),
    max_Half = max(Half, na.rm = TRUE),
    n = sum(!is.na(Half)),
    .groups = "drop"
  )
# table tells us % cover of veg in 0.5m is lower on average at random points than nests.

# test for statistical difference
anova_mod_.5 <- aov(Half ~ Season * Type, data = ex_df)
summary(anova_mod_.5)

emmeans(anova_mod_.5, pairwise ~ Season * Type)
# again it's by type (nest vs. random); no statistical significance with season

