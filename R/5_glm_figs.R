# Plotting selected models from logistic regression
# 1.a. All data - Summarized Microhabitat
# 1.b. All data - Functional Groups
# 1.c. All data - Individual species
# 2. Summer subset - 1 best model
# 3. Fall subset - 3 way tie

# load packages and functions and data
source('R/packages.R')
source('R/functions_glm.R')

all_df <- read_csv('data/glm_alldat.csv')
summer_df <- read_csv('data/glm_breedingdat.csv')
fall_df <- read_csv('data/glm_falldat.csv')

# prepare data by renaming the object type factor so that it's formatted for Wader Study
all_df$Object_Type <- factor(all_df$Object_Type)
levels(all_df$Object_Type)[levels(all_df$Object_Type) == "Non-veg Wood or Shell"] <- "Non-veg wood or shell"

fall_df$Object_Type <- factor(fall_df$Object_Type)
levels(fall_df$Object_Type)[levels(fall_df$Object_Type) == "Non-veg Wood or Shell"] <- "Non-veg wood or shell"

summer_df$Object_Type <- factor(summer_df$Object_Type)
levels(summer_df$Object_Type)[levels(summer_df$Object_Type) == "Non-veg Wood or Shell"] <- "Non-veg wood or shell"


# 1.a. All data - Summarized Microhabitat ----

# rename selected model and use glm for "predict" function later
f.mod.all.sm <- glm(Outcome ~ Object_Type + poly(TwoVeg,2) +  poly(Half,2), 
                        data = all_df, family = binomial)

# Plot covariate - TwoVeg
range <- all_df %>%
  filter(Outcome == 1) %>% # to limit to nest site values in plot
  summarise(min_TwoVeg = min(TwoVeg, na.rm = TRUE), max_TwoVeg = max(TwoVeg, na.rm = TRUE))

# Create a data frame with a sequence of values for variable of interest
# and separate rows for each level of `Object_Type`
new_data <- expand.grid(
  TwoVeg = seq(min(all_df$TwoVeg, na.rm = TRUE), max(all_df$TwoVeg, na.rm = TRUE), length.out = 100),
  Half = mean(all_df$Half, na.rm = TRUE),
  Object_Type = c("None", "Vegetation", "Non-veg wood or shell")
)

# Get predictions with standard errors
predictions <- predict(f.mod.all.sm, 
                       new_data, 
                       type = "link", 
                       se.fit = TRUE) 

# Calculate the predicted probabilities and confidence intervals
new_data$predicted_prob <- plogis(predictions$fit)
new_data$lower_ci <- plogis(predictions$fit - 1.96 * predictions$se.fit)
new_data$upper_ci <- plogis(predictions$fit + 1.96 * predictions$se.fit)

# Plot using ggplot2 with confidence intervals as a ribbon
TwoVeg_cover <- ggplot(new_data, aes(x = TwoVeg, y = predicted_prob, color = Object_Type)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = Object_Type), color = NA, alpha = 0.3) +
  labs(
    title = "% Cover Vegetation in 2m Radius",
    x = "Percent Cover",
    y = "Probability of nest",
    color = "Object Type",  # Change legend title for line color
    fill = "Object Type"    # Change legend title for ribbon fill
  ) +
  theme_classic() +
  xlim(range$min_TwoVeg, range$max_TwoVeg)

# Display the plot
print(TwoVeg_cover) 

# Plot covariate - Half
range <- all_df %>%
  filter(Outcome == 1) %>% # to limit to nest site values in plot
  summarise(min_Half = min(Half, na.rm = TRUE), max_Half = max(Half, na.rm = TRUE))

# Create a data frame with a sequence of values for variable of interest
# and separate rows for each level of `Object_Type`
new_data <- expand.grid(
  Half = seq(min(all_df$Half, na.rm = TRUE), max(all_df$Half, na.rm = TRUE), length.out = 100),
  TwoVeg = mean(all_df$TwoVeg, na.rm = TRUE),
  Object_Type = c("None", "Vegetation", "Non-veg wood or shell")
)

# Get predictions with standard errors
predictions <- predict(f.mod.all.sm, 
                       new_data, 
                       type = "link", 
                       se.fit = TRUE) 

# Calculate the predicted probabilities and confidence intervals
new_data$predicted_prob <- plogis(predictions$fit)
new_data$lower_ci <- plogis(predictions$fit - 1.96 * predictions$se.fit)
new_data$upper_ci <- plogis(predictions$fit + 1.96 * predictions$se.fit)

# Plot using ggplot2 with confidence intervals as a ribbon
Half_cover <- ggplot(new_data, aes(x = Half, y = predicted_prob, color = Object_Type)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = Object_Type), color = NA, alpha = 0.3) +
  labs(
    title = "% Cover Vegetation in 0.5m Radius",
    x = "Percent Cover",
    y = "Probability of nest",
    color = "Object Type",  # Change legend title for line color
    fill = "Object Type"    # Change legend title for ribbon fill
  ) +
  theme_classic() +
  xlim(range$min_Half, range$max_Half)

# Display the plot
print(Half_cover) 


# 1.b. All data - Functional Groups ----

# rename selected model and use glm for "predict" function later
f.mod.all.fg <- glm(Outcome ~ Object_Type + poly(Pioneer_Stabilizers,2) 
                    + Pioneer_Builders, 
                    data = all_df, family = binomial)

# Plot covariate - Pioneer Stabilizers
range <- all_df %>%
  filter(Outcome == 1) %>% # to limit to nest site values in plot
  summarise(min_PS = min(Pioneer_Stabilizers, na.rm = TRUE), 
            max_PS = max(Pioneer_Stabilizers, na.rm = TRUE))

# Create a data frame with a sequence of values for variable of interest
# and separate rows for each level of `Object_Type`
new_data <- expand.grid(
  Pioneer_Stabilizers = seq(min(all_df$Pioneer_Stabilizers, na.rm = TRUE), 
                            max(all_df$Pioneer_Stabilizers, na.rm = TRUE), 
                            length.out = 100),
  Pioneer_Builders = mean(all_df$Pioneer_Builders, na.rm = TRUE),
  Object_Type = c("None", "Vegetation", "Non-veg wood or shell")
)

# Get predictions with standard errors
predictions <- predict(f.mod.all.fg, 
                       new_data, 
                       type = "link", 
                       se.fit = TRUE) 

# Calculate the predicted probabilities and confidence intervals
new_data$predicted_prob <- plogis(predictions$fit)
new_data$lower_ci <- plogis(predictions$fit - 1.96 * predictions$se.fit)
new_data$upper_ci <- plogis(predictions$fit + 1.96 * predictions$se.fit)

# Plot using ggplot2 with confidence intervals as a ribbon
PS_cover <- ggplot(new_data, aes(x = Pioneer_Stabilizers, y = predicted_prob, color = Object_Type)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = Object_Type), color = NA, alpha = 0.3) +
  labs(
    title = "% Cover Pioneer Stabilizers in 2m Radius",
    x = "Percent Cover",
    y = "Probability of nest",
    color = "Object Type",  # Change legend title for line color
    fill = "Object Type"    # Change legend title for ribbon fill
  ) +
  theme_classic() +
  xlim(range$min_PS, range$max_PS)

# Display the plot
print(PS_cover) 

# Plot covariate - Pioneer Builders
range <- all_df %>%
  filter(Outcome == 1) %>% # to limit to nest site values in plot
  summarise(min_PB = min(Pioneer_Builders, na.rm = TRUE), 
            max_PB = max(Pioneer_Builders, na.rm = TRUE))

# Create a data frame with a sequence of values for variable of interest
# and separate rows for each level of `Object_Type`
new_data <- expand.grid(
  Pioneer_Builders = seq(min(all_df$Pioneer_Builders, na.rm = TRUE), 
                            max(all_df$Pioneer_Builders, na.rm = TRUE), 
                            length.out = 100),
  Pioneer_Stabilizers = mean(all_df$Pioneer_Stabilizers, na.rm = TRUE),
  Object_Type = c("None", "Vegetation", "Non-veg wood or shell")
)

# Get predictions with standard errors
predictions <- predict(f.mod.all.fg, 
                       new_data, 
                       type = "link", 
                       se.fit = TRUE) 

# Calculate the predicted probabilities and confidence intervals
new_data$predicted_prob <- plogis(predictions$fit)
new_data$lower_ci <- plogis(predictions$fit - 1.96 * predictions$se.fit)
new_data$upper_ci <- plogis(predictions$fit + 1.96 * predictions$se.fit)

# Plot using ggplot2 with confidence intervals as a ribbon
PB_cover <- ggplot(new_data, aes(x = Pioneer_Builders, y = predicted_prob, color = Object_Type)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = Object_Type), color = NA, alpha = 0.3) +
  labs(
    title = "% Cover Pioneer Builders in 2m Radius",
    x = "Percent Cover",
    y = "Probability of nest",
    color = "Object Type",  # Change legend title for line color
    fill = "Object Type"    # Change legend title for ribbon fill
  ) +
  theme_classic() +
  xlim(range$min_PB, range$max_PB)

# Display the plot
print(PB_cover) 













# 1.c. All data - Individual Species ----

# rename selected model and use glm for "predict" function later
f.mod.all.is <- glm(Outcome ~ Object_Type + poly(BS,2) + IP + poly(SR,2), 
                    data = all_df, family = binomial)

# Plot covariate - BS
range <- all_df %>%
  filter(Outcome == 1) %>% # to limit to nest site values in plot
  summarise(min_BS = min(BS, na.rm = TRUE), 
            max_BS = max(BS, na.rm = TRUE))

# Create a data frame with a sequence of values for variable of interest
# and separate rows for each level of `Object_Type`
new_data <- expand.grid(
  BS = seq(min(all_df$BS, na.rm = TRUE), max(all_df$BS, na.rm = TRUE), 
           length.out = 100),
  SR = mean(all_df$SR, na.rm = TRUE),
  IP = mean(all_df$IP, na.rm = TRUE),
  Object_Type = c("None", "Vegetation", "Non-veg wood or shell")
)

# Get predictions with standard errors
predictions <- predict(f.mod.all.is, 
                       new_data, 
                       type = "link", 
                       se.fit = TRUE) 

# Calculate the predicted probabilities and confidence intervals
new_data$predicted_prob <- plogis(predictions$fit)
new_data$lower_ci <- plogis(predictions$fit - 1.96 * predictions$se.fit)
new_data$upper_ci <- plogis(predictions$fit + 1.96 * predictions$se.fit)

# Plot using ggplot2 with confidence intervals as a ribbon
BS_cover <- ggplot(new_data, aes(x = BS, y = predicted_prob, color = Object_Type)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = Object_Type), color = NA, alpha = 0.3) +
  labs(
    title = "% Cover Beach Bur Sage in 2m Radius",
    x = "Percent Cover",
    y = "Probability of nest",
    color = "Object Type",  # Change legend title for line color
    fill = "Object Type"    # Change legend title for ribbon fill
  ) +
  theme_classic() +
  xlim(range$min_BS, range$max_BS)

# Display the plot
print(BS_cover) 


# Plot covariate - SR
range <- all_df %>%
  filter(Outcome == 1) %>% # to limit to nest site values in plot
  summarise(min_SR = min(SR, na.rm = TRUE), 
            max_SR = max(SR, na.rm = TRUE))

# Create a data frame with a sequence of values for variable of interest
# and separate rows for each level of `Object_Type`
new_data <- expand.grid(
  SR = seq(min(all_df$SR, na.rm = TRUE), 
                         max(all_df$SR, na.rm = TRUE), 
                         length.out = 100),
  BS = mean(all_df$BS, na.rm = TRUE),
  IP = mean(all_df$IP, na.rm = TRUE),
  Object_Type = c("None", "Vegetation", "Non-veg wood or shell")
)

# Get predictions with standard errors
predictions <- predict(f.mod.all.is, 
                       new_data, 
                       type = "link", 
                       se.fit = TRUE) 

# Calculate the predicted probabilities and confidence intervals
new_data$predicted_prob <- plogis(predictions$fit)
new_data$lower_ci <- plogis(predictions$fit - 1.96 * predictions$se.fit)
new_data$upper_ci <- plogis(predictions$fit + 1.96 * predictions$se.fit)

# Plot using ggplot2 with confidence intervals as a ribbon
SR_cover <- ggplot(new_data, aes(x = SR, y = predicted_prob, color = Object_Type)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = Object_Type), color = NA, alpha = 0.3) +
  labs(
    title = "% Cover Sea Rocket in 2m Radius",
    x = "Percent Cover",
    y = "Probability of nest",
    color = "Object Type",  # Change legend title for line color
    fill = "Object Type"    # Change legend title for ribbon fill
  ) +
  theme_classic() +
  xlim(range$min_SR, range$max_SR)

# Display the plot
print(SR_cover) 


# Plot covariate - IP
range <- all_df %>%
  filter(Outcome == 1) %>% # to limit to nest site values in plot
  summarise(min_IP = min(IP, na.rm = TRUE), 
            max_IP = max(IP, na.rm = TRUE))

# Create a data frame with a sequence of values for variable of interest
# and separate rows for each level of `Object_Type`
new_data <- expand.grid(
  IP = seq(min(all_df$IP, na.rm = TRUE), 
           max(all_df$IP, na.rm = TRUE), 
           length.out = 100),
  BS = mean(all_df$BS, na.rm = TRUE),
  SR = mean(all_df$SR, na.rm = TRUE),
  Object_Type = c("None", "Vegetation", "Non-veg wood or shell")
)

# Get predictions with standard errors
predictions <- predict(f.mod.all.is, 
                       new_data, 
                       type = "link", 
                       se.fit = TRUE) 

# Calculate the predicted probabilities and confidence intervals
new_data$predicted_prob <- plogis(predictions$fit)
new_data$lower_ci <- plogis(predictions$fit - 1.96 * predictions$se.fit)
new_data$upper_ci <- plogis(predictions$fit + 1.96 * predictions$se.fit)

# Plot using ggplot2 with confidence intervals as a ribbon
IP_cover <- ggplot(new_data, aes(x = IP, y = predicted_prob, color = Object_Type)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = Object_Type), color = NA, alpha = 0.3) +
  labs(
    title = "% Cover Ice Plant in 2m Radius",
    x = "Percent Cover",
    y = "Probability of nest",
    color = "Object Type",  # Change legend title for line color
    fill = "Object Type"    # Change legend title for ribbon fill
  ) +
  theme_classic() +
  xlim(range$min_IP, range$max_IP)

# Display the plot
print(IP_cover) 












# 1.d. Facet Plot for All Data ----

# Remove legends and apply correct titles/axis labels
TwoVegp <- TwoVeg_cover + theme(legend.position = "none",
                                plot.title = element_text(size = 11),
                                axis.title.x = element_text(size = 10), 
                                axis.title.y = element_text(size = 10)) +
  scale_color_grey() + scale_fill_grey() +
  font_theme + 
  labs(title = "(a) Vegetative cover in 2 m radius", x = "Percent cover in 2 m radius", 
       y = "Probability of nest")

Halfp <- Half_cover + theme(legend.position = "none") +
  scale_color_grey() + scale_fill_grey() +
  font_theme + 
  labs(title = "(b) Vegetative cover in 0.5 m radius", 
       x = "Percent cover in 0.5 m radius")

PSp <- PS_cover + theme(legend.position = "none", 
                        plot.title = element_text(size = 11), 
                        axis.title.x = element_text(size = 10), 
                        axis.title.y = element_text(size = 10)) +
  scale_color_grey() + scale_fill_grey() +
  labs(title = "(a) Pioneer stabilizers", x = "Percent cover in 2 m radius", 
       y = "Probability of nest")

PBp <- PB_cover + theme(legend.position = "none") +
  scale_color_grey() + scale_fill_grey() +
  font_theme + 
  labs(title = "(b) Pioneer builders", x = "Percent cover in 2 m radius", 
       y = "Probability of nest")

BSp <- BS_cover + theme(legend.position = "none", 
                        plot.title = element_text(size = 11), 
                        axis.title.x = element_text(size = 10), 
                        axis.title.y = element_text(size = 10)) +
  scale_color_grey() + scale_fill_grey() +
  font_theme + 
  labs(title = "(a) Beach bur sage", x = "Percent cover in 2 m radius", 
       y = "Probability of nest")

SRp <- SR_cover + theme(legend.position = "none") +
  scale_color_grey() + scale_fill_grey() +
  font_theme + 
  labs(title = "(b) Sea rocket", x = "Percent cover in 2 m radius", 
       y = "Probability of nest")

IPp <- IP_cover + theme(legend.position = "none") +
  scale_color_grey() + scale_fill_grey() +
  font_theme + 
  labs(title = "(c) Dead ice plant", x = "Percent cover in 2 m radius", 
       y = "Probability of nest")


# Extract single legend before removing it from the plot
legend <- get_legend(
  BSp +
    theme(
      legend.position = "right",
      legend.direction = "horizontal",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 10)
    ) +
    guides(color = guide_legend(nrow=1))
)

# Make an empty plot for spacing
empty_plot <- ggplot() + theme_void()

# Create label plots
label_top    <- make_row_label("Summarized Microhabitat")
label_middle <- make_row_label("Functional Groups")
label_bottom <- make_row_label("Individual Species")

# Example top row: label + plots
top_row <- plot_grid(
  label_top,
  plot_grid(TwoVegp, Halfp, empty_plot, nrow = 1, rel_widths = c(1,1,1)),
  ncol = 2,
  rel_widths = c(0.05, 0.95)  # narrow label, wide plots
)

middle_row <- plot_grid(
  label_middle,
  plot_grid(PSp, PBp, empty_plot, nrow = 1, rel_widths = c(1,1,1)),
  ncol = 2,
  rel_widths = c(0.05, 0.95)
)

bottom_row <- plot_grid(
  label_bottom,
  plot_grid(BSp, SRp, IPp, nrow = 1, rel_widths = c(1,1,1)),
  ncol = 2,
  rel_widths = c(0.05, 0.95)
)

# Combine all rows with legend
legend_row <- plot_grid(legend, nrow = 1)
facet_alldata_3groups <- plot_grid(top_row, middle_row, bottom_row, legend_row, ncol = 1, rel_heights = c(1,1,1,0.2))

# Display
print(final_plot)

ggsave("fig/facet_alldata_3groups.jpg", plot = facet_alldata_3groups, width = 170, height = 170, units = "mm", dpi = 300) #adjust size as needed























# 2.a. Summer subset - Summarized Microhabitat ----

# rename selected model and use glm for "predict" function later
f.mod.summer.sm <- glm(Outcome ~ Object_Type + poly(TwoVeg,2) +  poly(Half,2), 
                    data = summer_df, family = binomial)

# Plot covariate - TwoVeg
range <- summer_df %>%
  filter(Outcome == 1) %>% # to limit to nest site values in plot
  summarise(min_TwoVeg = min(TwoVeg, na.rm = TRUE), max_TwoVeg = max(TwoVeg, na.rm = TRUE))

# Create a data frame with a sequence of values for variable of interest
# and separate rows for each level of `Object_Type`
new_data <- expand.grid(
  TwoVeg = seq(min(summer_df$TwoVeg, na.rm = TRUE), max(summer_df$TwoVeg, na.rm = TRUE), length.out = 100),
  Half = mean(summer_df$Half, na.rm = TRUE),
  Object_Type = c("None", "Vegetation", "Non-veg wood or shell")
)

# Get predictions with standard errors
predictions <- predict(f.mod.summer.sm, 
                       new_data, 
                       type = "link", 
                       se.fit = TRUE) 

# Calculate the predicted probabilities and confidence intervals
new_data$predicted_prob <- plogis(predictions$fit)
new_data$lower_ci <- plogis(predictions$fit - 1.96 * predictions$se.fit)
new_data$upper_ci <- plogis(predictions$fit + 1.96 * predictions$se.fit)

# Plot using ggplot2 with confidence intervals as a ribbon
TwoVeg_cover <- ggplot(new_data, aes(x = TwoVeg, y = predicted_prob, color = Object_Type)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = Object_Type), color = NA, alpha = 0.3) +
  labs(
    title = "% Cover Vegetation in 2m Radius",
    x = "Percent Cover",
    y = "Probability of nest",
    color = "Object Type",  # Change legend title for line color
    fill = "Object Type"    # Change legend title for ribbon fill
  ) +
  theme_classic() +
  xlim(range$min_TwoVeg, range$max_TwoVeg)

# Display the plot
print(TwoVeg_cover) 

# Plot covariate - Half
range <- summer_df %>%
  filter(Outcome == 1) %>% # to limit to nest site values in plot
  summarise(min_Half = min(Half, na.rm = TRUE), max_Half = max(Half, na.rm = TRUE))

# Create a data frame with a sequence of values for variable of interest
# and separate rows for each level of `Object_Type`
new_data <- expand.grid(
  Half = seq(min(summer_df$Half, na.rm = TRUE), max(summer_df$Half, na.rm = TRUE), length.out = 100),
  TwoVeg = mean(summer_df$TwoVeg, na.rm = TRUE),
  Object_Type = c("None", "Vegetation", "Non-veg wood or shell")
)

# Get predictions with standard errors
predictions <- predict(f.mod.summer.sm, 
                       new_data, 
                       type = "link", 
                       se.fit = TRUE) 

# Calculate the predicted probabilities and confidence intervals
new_data$predicted_prob <- plogis(predictions$fit)
new_data$lower_ci <- plogis(predictions$fit - 1.96 * predictions$se.fit)
new_data$upper_ci <- plogis(predictions$fit + 1.96 * predictions$se.fit)

# Plot using ggplot2 with confidence intervals as a ribbon
Half_cover <- ggplot(new_data, aes(x = Half, y = predicted_prob, color = Object_Type)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = Object_Type), color = NA, alpha = 0.3) +
  labs(
    title = "% Cover Vegetation in 0.5m Radius",
    x = "Percent Cover",
    y = "Probability of nest",
    color = "Object Type",  # Change legend title for line color
    fill = "Object Type"    # Change legend title for ribbon fill
  ) +
  theme_classic() +
  xlim(range$min_Half, range$max_Half)

# Display the plot
print(Half_cover) 


# 2.b. Summer subset - Functional Groups ----

# rename selected model and use glm for "predict" function later
f.mod.summer.fg <- glm(Outcome ~ Object_Type + Pioneer_Stabilizers, 
                    data = summer_df, family = binomial)

# Plot covariate - Pioneer Stabilizers
range <- summer_df %>%
  filter(Outcome == 1) %>% # to limit to nest site values in plot
  summarise(min_PS = min(Pioneer_Stabilizers, na.rm = TRUE), 
            max_PS = max(Pioneer_Stabilizers, na.rm = TRUE))

# Create a data frame with a sequence of values for variable of interest
# and separate rows for each level of `Object_Type`
new_data <- expand.grid(
  Pioneer_Stabilizers = seq(min(summer_df$Pioneer_Stabilizers, na.rm = TRUE), 
                            max(summer_df$Pioneer_Stabilizers, na.rm = TRUE), 
                            length.out = 100),
  Object_Type = c("None", "Vegetation", "Non-veg wood or shell")
)

# Get predictions with standard errors
predictions <- predict(f.mod.summer.fg, 
                       new_data, 
                       type = "link", 
                       se.fit = TRUE) 

# Calculate the predicted probabilities and confidence intervals
new_data$predicted_prob <- plogis(predictions$fit)
new_data$lower_ci <- plogis(predictions$fit - 1.96 * predictions$se.fit)
new_data$upper_ci <- plogis(predictions$fit + 1.96 * predictions$se.fit)

# Plot using ggplot2 with confidence intervals as a ribbon
PS_cover <- ggplot(new_data, aes(x = Pioneer_Stabilizers, y = predicted_prob, color = Object_Type)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = Object_Type), color = NA, alpha = 0.3) +
  labs(
    title = "% Cover Pioneer Stabilizers in 2m Radius",
    x = "Percent Cover",
    y = "Probability of nest",
    color = "Object Type",  # Change legend title for line color
    fill = "Object Type"    # Change legend title for ribbon fill
  ) +
  theme_classic() +
  xlim(range$min_PS, range$max_PS)

# Display the plot
print(PS_cover) 

# 2.c. Summer subset - Individual Species ----

# rename selected model and use glm for "predict" function later
f.mod.summer.is <- glm(Outcome ~ Object_Type + BS + poly(SR,2), 
                    data = summer_df, family = binomial)

# Plot covariate - BS
range <- summer_df %>%
  filter(Outcome == 1) %>% # to limit to nest site values in plot
  summarise(min_BS = min(BS, na.rm = TRUE), 
            max_BS = max(BS, na.rm = TRUE))

# Create a data frame with a sequence of values for variable of interest
# and separate rows for each level of `Object_Type`
new_data <- expand.grid(
  BS = seq(min(summer_df$BS, na.rm = TRUE), max(summer_df$BS, na.rm = TRUE), 
           length.out = 100),
  SR = mean(summer_df$SR, na.rm = TRUE),
  Object_Type = c("None", "Vegetation", "Non-veg wood or shell")
)

# Get predictions with standard errors
predictions <- predict(f.mod.summer.is, 
                       new_data, 
                       type = "link", 
                       se.fit = TRUE) 

# Calculate the predicted probabilities and confidence intervals
new_data$predicted_prob <- plogis(predictions$fit)
new_data$lower_ci <- plogis(predictions$fit - 1.96 * predictions$se.fit)
new_data$upper_ci <- plogis(predictions$fit + 1.96 * predictions$se.fit)

# Plot using ggplot2 with confidence intervals as a ribbon
BS_cover <- ggplot(new_data, aes(x = BS, y = predicted_prob, color = Object_Type)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = Object_Type), color = NA, alpha = 0.3) +
  labs(
    title = "% Cover Beach Bur Sage in 2m Radius",
    x = "Percent Cover",
    y = "Probability of nest",
    color = "Object Type",  # Change legend title for line color
    fill = "Object Type"    # Change legend title for ribbon fill
  ) +
  theme_classic() +
  xlim(range$min_BS, range$max_BS)

# Display the plot
print(BS_cover) 


# Plot covariate - SR
range <- summer_df %>%
  filter(Outcome == 1) %>% # to limit to nest site values in plot
  summarise(min_SR = min(SR, na.rm = TRUE), 
            max_SR = max(SR, na.rm = TRUE))

# Create a data frame with a sequence of values for variable of interest
# and separate rows for each level of `Object_Type`
new_data <- expand.grid(
  SR = seq(min(summer_df$SR, na.rm = TRUE), 
           max(summer_df$SR, na.rm = TRUE), 
           length.out = 100),
  BS = mean(summer_df$BS, na.rm = TRUE),
  Object_Type = c("None", "Vegetation", "Non-veg wood or shell")
)

# Get predictions with standard errors
predictions <- predict(f.mod.summer.is, 
                       new_data, 
                       type = "link", 
                       se.fit = TRUE) 

# Calculate the predicted probabilities and confidence intervals
new_data$predicted_prob <- plogis(predictions$fit)
new_data$lower_ci <- plogis(predictions$fit - 1.96 * predictions$se.fit)
new_data$upper_ci <- plogis(predictions$fit + 1.96 * predictions$se.fit)

# Plot using ggplot2 with confidence intervals as a ribbon
SR_cover <- ggplot(new_data, aes(x = SR, y = predicted_prob, color = Object_Type)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = Object_Type), color = NA, alpha = 0.3) +
  labs(
    title = "% Cover Sea Rocket in 2m Radius",
    x = "Percent Cover",
    y = "Probability of nest",
    color = "Object Type",  # Change legend title for line color
    fill = "Object Type"    # Change legend title for ribbon fill
  ) +
  theme_classic() +
  xlim(range$min_SR, range$max_SR)

# Display the plot
print(SR_cover) 

# 2.d. Facet Plot for Summer subset ----

# Remove legends and apply correct titles/axis labels
TwoVegp <- TwoVeg_cover + theme(legend.position = "none",
                                plot.title = element_text(size = 11),
                                axis.title.x = element_text(size = 10), 
                                axis.title.y = element_text(size = 10)) +
  scale_color_grey() + scale_fill_grey() +
  font_theme + 
  labs(title = "(a) Vegetative cover in 2 m radius", x = "Percent cover in 2 m radius", 
       y = "Probability of nest")

Halfp <- Half_cover + theme(legend.position = "none") +
  scale_color_grey() + scale_fill_grey() +
  font_theme + 
  labs(title = "(b) Vegetative cover in 0.5 m radius", 
       x = "Percent cover in 0.5 m radius")

PSp <- PS_cover + theme(legend.position = "none", 
                        plot.title = element_text(size = 11), 
                        axis.title.x = element_text(size = 10), 
                        axis.title.y = element_text(size = 10)) +
  scale_color_grey() + scale_fill_grey() +
  labs(title = "(a) Pioneer stabilizers", x = "Percent cover in 2 m radius", 
       y = "Probability of nest")

BSp <- BS_cover + theme(legend.position = "none", 
                        plot.title = element_text(size = 11), 
                        axis.title.x = element_text(size = 10), 
                        axis.title.y = element_text(size = 10)) +
  scale_color_grey() + scale_fill_grey() +
  font_theme + 
  labs(title = "(a) Beach bur sage", x = "Percent cover in 2 m radius", 
       y = "Probability of nest")

SRp <- SR_cover + theme(legend.position = "none") +
  scale_color_grey() + scale_fill_grey() +
  font_theme + 
  labs(title = "(b) Sea rocket", x = "Percent cover in 2 m radius", 
       y = "Probability of nest")

# Extract single legend before removing it from the plot
legend <- get_legend(
  BSp +
    theme(
      legend.position = "right",
      legend.direction = "horizontal",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 10)
    ) +
    guides(color = guide_legend(nrow=1))
)

# Make an empty plot for spacing
empty_plot <- ggplot() + theme_void()

# Create label plots
label_top    <- make_row_label("Summarized Microhabitat")
label_middle <- make_row_label("Functional Groups")
label_bottom <- make_row_label("Individual Species")

# Example top row: label + plots
top_row <- plot_grid(
  label_top,
  plot_grid(TwoVegp, Halfp, nrow = 1, rel_widths = c(1,1)),
  ncol = 2,
  rel_widths = c(0.05, 0.95)  # narrow label, wide plots
)

middle_row <- plot_grid(
  label_middle,
  plot_grid(PSp, empty_plot, nrow = 1, rel_widths = c(1,1)),
  ncol = 2,
  rel_widths = c(0.05, 0.95)
)

bottom_row <- plot_grid(
  label_bottom,
  plot_grid(BSp, SRp, nrow = 1, rel_widths = c(1,1)),
  ncol = 2,
  rel_widths = c(0.05, 0.95)
)

# Combine all rows with legend
legend_row <- plot_grid(legend, nrow = 1)
facet_summerdata_3groups <- plot_grid(top_row, middle_row, bottom_row, legend_row, ncol = 1, rel_heights = c(1,1,1,0.2))

# Display
print(facet_summerdata_3groups)

ggsave("fig/facet_summerdata_3groups.jpg", plot = facet_summerdata_3groups, width = 170, height = 170, units = "mm", dpi = 300) #adjust size as needed
























# 3.a. Fall data - Summarized Microhabitat ----

# rename selected model and use glm for "predict" function later
f.mod.fall.sm <- glm(Outcome ~ Object_Type + poly(TwoVeg,2) +  Viewshed_sc, 
                    data = fall_df, family = binomial)

# Plot covariate - TwoVeg
range <- fall_df %>%
  filter(Outcome == 1) %>% # to limit to nest site values in plot
  summarise(min_TwoVeg = min(TwoVeg, na.rm = TRUE), max_TwoVeg = max(TwoVeg, na.rm = TRUE))

# Create a data frame with a sequence of values for variable of interest
# and separate rows for each level of `Object_Type`
new_data <- expand.grid(
  TwoVeg = seq(min(fall_df$TwoVeg, na.rm = TRUE), 
               max(fall_df$TwoVeg, na.rm = TRUE), length.out = 100),
  Viewshed_sc = mean(fall_df$Viewshed_sc, na.rm = TRUE),
  Object_Type = c("None", "Vegetation", "Non-veg wood or shell")
)

# Get predictions with standard errors
predictions <- predict(f.mod.fall.sm, 
                       new_data, 
                       type = "link", 
                       se.fit = TRUE) 

# Calculate the predicted probabilities and confidence intervals
new_data$predicted_prob <- plogis(predictions$fit)
new_data$lower_ci <- plogis(predictions$fit - 1.96 * predictions$se.fit)
new_data$upper_ci <- plogis(predictions$fit + 1.96 * predictions$se.fit)

# Plot using ggplot2 with confidence intervals as a ribbon
TwoVeg_cover <- ggplot(new_data, aes(x = TwoVeg, y = predicted_prob, color = Object_Type)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = Object_Type), color = NA, alpha = 0.3) +
  labs(
    title = "% Cover Vegetation in 2m Radius",
    x = "Percent Cover",
    y = "Probability of nest",
    color = "Object Type",  # Change legend title for line color
    fill = "Object Type"    # Change legend title for ribbon fill
  ) +
  theme_classic() +
  xlim(range$min_TwoVeg, range$max_TwoVeg)

# Display the plot
print(TwoVeg_cover) 

# Plot covariate - Viewshed_sc
range <- fall_df %>%
  filter(Outcome == 1) %>% # to limit to nest site values in plot
  summarise(min_View = min(Viewshed_sc, na.rm = TRUE), max_View = max(Viewshed_sc, na.rm = TRUE))

# Create a data frame with a sequence of values for variable of interest
# and separate rows for each level of `Object_Type`
new_data <- expand.grid(
  Viewshed_sc = seq(min(fall_df$Viewshed_sc, na.rm = TRUE), max(fall_df$Viewshed_sc, na.rm = TRUE), length.out = 100),
  TwoVeg = mean(fall_df$TwoVeg, na.rm = TRUE),
  Object_Type = c("None", "Vegetation", "Non-veg wood or shell")
)

# Get predictions with standard errors
predictions <- predict(f.mod.fall.sm, 
                       new_data, 
                       type = "link", 
                       se.fit = TRUE) 

# Calculate the predicted probabilities and confidence intervals
new_data$predicted_prob <- plogis(predictions$fit)
new_data$lower_ci <- plogis(predictions$fit - 1.96 * predictions$se.fit)
new_data$upper_ci <- plogis(predictions$fit + 1.96 * predictions$se.fit)

# Plot using ggplot2 with confidence intervals as a ribbon
View_sm <- ggplot(new_data, aes(x = Viewshed_sc, y = predicted_prob, color = Object_Type)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = Object_Type), color = NA, alpha = 0.3) +
  labs(
    title = "Viewshed from nest",
    x = "Proportion unobscured",
    y = "Probability of nest",
    color = "Object Type",  # Change legend title for line color
    fill = "Object Type"    # Change legend title for ribbon fill
  ) +
  theme_classic() +
  xlim(range$min_View, range$max_View)

# Display the plot
print(View_sm) 


# 3.b. Fall data - Functional Groups ----

# rename selected model and use glm for "predict" function later
f.mod.fall.fg <- glm(Outcome ~ Object_Type + poly(Pioneer_Stabilizers,2) 
                    + Pioneer_Builders + Viewshed_sc, 
                    data = fall_df, family = binomial)

# Plot covariate - Pioneer Stabilizers
range <- fall_df %>%
  filter(Outcome == 1) %>% # to limit to nest site values in plot
  summarise(min_PS = min(Pioneer_Stabilizers, na.rm = TRUE), 
            max_PS = max(Pioneer_Stabilizers, na.rm = TRUE))

# Create a data frame with a sequence of values for variable of interest
# and separate rows for each level of `Object_Type`
new_data <- expand.grid(
  Pioneer_Stabilizers = seq(min(fall_df$Pioneer_Stabilizers, na.rm = TRUE), 
                            max(fall_df$Pioneer_Stabilizers, na.rm = TRUE), 
                            length.out = 100),
  Pioneer_Builders = mean(fall_df$Pioneer_Builders, na.rm = TRUE),
  Viewshed_sc = mean(fall_df$Viewshed_sc, na.rm = TRUE),
  Object_Type = c("None", "Vegetation", "Non-veg wood or shell")
)

# Get predictions with standard errors
predictions <- predict(f.mod.fall.fg, 
                       new_data, 
                       type = "link", 
                       se.fit = TRUE) 

# Calculate the predicted probabilities and confidence intervals
new_data$predicted_prob <- plogis(predictions$fit)
new_data$lower_ci <- plogis(predictions$fit - 1.96 * predictions$se.fit)
new_data$upper_ci <- plogis(predictions$fit + 1.96 * predictions$se.fit)

# Plot using ggplot2 with confidence intervals as a ribbon
PS_cover <- ggplot(new_data, aes(x = Pioneer_Stabilizers, y = predicted_prob, color = Object_Type)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = Object_Type), color = NA, alpha = 0.3) +
  labs(
    title = "% Cover Pioneer Stabilizers in 2m Radius",
    x = "Percent Cover",
    y = "Probability of nest",
    color = "Object Type",  # Change legend title for line color
    fill = "Object Type"    # Change legend title for ribbon fill
  ) +
  theme_classic() +
  xlim(range$min_PS, range$max_PS)

# Display the plot
print(PS_cover) 

# Plot covariate - Pioneer Builders
range <- fall_df %>%
  filter(Outcome == 1) %>% # to limit to nest site values in plot
  summarise(min_PB = min(Pioneer_Builders, na.rm = TRUE), 
            max_PB = max(Pioneer_Builders, na.rm = TRUE))

# Create a data frame with a sequence of values for variable of interest
# and separate rows for each level of `Object_Type`
new_data <- expand.grid(
  Pioneer_Builders = seq(min(fall_df$Pioneer_Builders, na.rm = TRUE), 
                         max(fall_df$Pioneer_Builders, na.rm = TRUE), 
                         length.out = 100),
  Pioneer_Stabilizers = mean(fall_df$Pioneer_Stabilizers, na.rm = TRUE),
  Viewshed_sc = mean(fall_df$Viewshed_sc, na.rm = TRUE),
  Object_Type = c("None", "Vegetation", "Non-veg wood or shell")
)

# Get predictions with standard errors
predictions <- predict(f.mod.fall.fg, 
                       new_data, 
                       type = "link", 
                       se.fit = TRUE) 

# Calculate the predicted probabilities and confidence intervals
new_data$predicted_prob <- plogis(predictions$fit)
new_data$lower_ci <- plogis(predictions$fit - 1.96 * predictions$se.fit)
new_data$upper_ci <- plogis(predictions$fit + 1.96 * predictions$se.fit)

# Plot using ggplot2 with confidence intervals as a ribbon
PB_cover <- ggplot(new_data, aes(x = Pioneer_Builders, y = predicted_prob, color = Object_Type)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = Object_Type), color = NA, alpha = 0.3) +
  labs(
    title = "% Cover Pioneer Builders in 2m Radius",
    x = "Percent Cover",
    y = "Probability of nest",
    color = "Object Type",  # Change legend title for line color
    fill = "Object Type"    # Change legend title for ribbon fill
  ) +
  theme_classic() +
  xlim(range$min_PB, range$max_PB)

# Display the plot
print(PB_cover) 


# Plot covariate - Viewshed
range <- fall_df %>%
  filter(Outcome == 1) %>% # to limit to nest site values in plot
  summarise(min_View = min(Viewshed_sc, na.rm = TRUE), 
            max_View = max(Viewshed_sc, na.rm = TRUE))

# Create a data frame with a sequence of values for variable of interest
# and separate rows for each level of `Object_Type`
new_data <- expand.grid(
  Viewshed_sc = seq(min(fall_df$Viewshed_sc, na.rm = TRUE), 
                         max(fall_df$Viewshed_sc, na.rm = TRUE), 
                         length.out = 100),
  Pioneer_Stabilizers = mean(fall_df$Pioneer_Stabilizers, na.rm = TRUE),
  Pioneer_Builders = mean(fall_df$Pioneer_Builders, na.rm = TRUE),
  Object_Type = c("None", "Vegetation", "Non-veg wood or shell")
)

# Get predictions with standard errors
predictions <- predict(f.mod.fall.fg, 
                       new_data, 
                       type = "link", 
                       se.fit = TRUE) 

# Calculate the predicted probabilities and confidence intervals
new_data$predicted_prob <- plogis(predictions$fit)
new_data$lower_ci <- plogis(predictions$fit - 1.96 * predictions$se.fit)
new_data$upper_ci <- plogis(predictions$fit + 1.96 * predictions$se.fit)

# Plot using ggplot2 with confidence intervals as a ribbon
View_fg <- ggplot(new_data, aes(x = Viewshed_sc, y = predicted_prob, color = Object_Type)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = Object_Type), color = NA, alpha = 0.3) +
  labs(
    title = "Viewshed",
    x = "Proportion onobscured",
    y = "Probability of nest",
    color = "Object Type",  # Change legend title for line color
    fill = "Object Type"    # Change legend title for ribbon fill
  ) +
  theme_classic() +
  xlim(range$min_View, range$max_View)

# Display the plot
print(View_fg) 










# 3.c. Fall data - Individual Species ----

# rename selected model and use glm for "predict" function later
f.mod.fall.is <- glm(Outcome ~ Object_Type + poly(BS,2) + IP + SR, 
                    data = fall_df, family = binomial)

# Plot covariate - BS
range <- fall_df %>%
  filter(Outcome == 1) %>% # to limit to nest site values in plot
  summarise(min_BS = min(BS, na.rm = TRUE), 
            max_BS = max(BS, na.rm = TRUE))

# Create a data frame with a sequence of values for variable of interest
# and separate rows for each level of `Object_Type`
new_data <- expand.grid(
  BS = seq(min(fall_df$BS, na.rm = TRUE), max(fall_df$BS, na.rm = TRUE), 
           length.out = 100),
  SR = mean(fall_df$SR, na.rm = TRUE),
  IP = mean(fall_df$IP, na.rm = TRUE),
  Object_Type = c("None", "Vegetation", "Non-veg wood or shell")
)

# Get predictions with standard errors
predictions <- predict(f.mod.fall.is, 
                       new_data, 
                       type = "link", 
                       se.fit = TRUE) 

# Calculate the predicted probabilities and confidence intervals
new_data$predicted_prob <- plogis(predictions$fit)
new_data$lower_ci <- plogis(predictions$fit - 1.96 * predictions$se.fit)
new_data$upper_ci <- plogis(predictions$fit + 1.96 * predictions$se.fit)

# Plot using ggplot2 with confidence intervals as a ribbon
BS_cover <- ggplot(new_data, aes(x = BS, y = predicted_prob, color = Object_Type)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = Object_Type), color = NA, alpha = 0.3) +
  labs(
    title = "% Cover Beach Bur Sage in 2m Radius",
    x = "Percent Cover",
    y = "Probability of nest",
    color = "Object Type",  # Change legend title for line color
    fill = "Object Type"    # Change legend title for ribbon fill
  ) +
  theme_classic() +
  xlim(range$min_BS, range$max_BS)

# Display the plot
print(BS_cover) 


# Plot covariate - SR
range <- fall_df %>%
  filter(Outcome == 1) %>% # to limit to nest site values in plot
  summarise(min_SR = min(SR, na.rm = TRUE), 
            max_SR = max(SR, na.rm = TRUE))

# Create a data frame with a sequence of values for variable of interest
# and separate rows for each level of `Object_Type`
new_data <- expand.grid(
  SR = seq(min(fall_df$SR, na.rm = TRUE), 
           max(fall_df$SR, na.rm = TRUE), 
           length.out = 100),
  BS = mean(fall_df$BS, na.rm = TRUE),
  IP = mean(fall_df$IP, na.rm = TRUE),
  Object_Type = c("None", "Vegetation", "Non-veg wood or shell")
)

# Get predictions with standard errors
predictions <- predict(f.mod.fall.is, 
                       new_data, 
                       type = "link", 
                       se.fit = TRUE) 

# Calculate the predicted probabilities and confidence intervals
new_data$predicted_prob <- plogis(predictions$fit)
new_data$lower_ci <- plogis(predictions$fit - 1.96 * predictions$se.fit)
new_data$upper_ci <- plogis(predictions$fit + 1.96 * predictions$se.fit)

# Plot using ggplot2 with confidence intervals as a ribbon
SR_cover <- ggplot(new_data, aes(x = SR, y = predicted_prob, color = Object_Type)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = Object_Type), color = NA, alpha = 0.3) +
  labs(
    title = "% Cover Sea Rocket in 2m Radius",
    x = "Percent Cover",
    y = "Probability of nest",
    color = "Object Type",  # Change legend title for line color
    fill = "Object Type"    # Change legend title for ribbon fill
  ) +
  theme_classic() +
  xlim(range$min_SR, range$max_SR)

# Display the plot
print(SR_cover) 


# Plot covariate - IP
range <- fall_df %>%
  filter(Outcome == 1) %>% # to limit to nest site values in plot
  summarise(min_IP = min(IP, na.rm = TRUE), 
            max_IP = max(IP, na.rm = TRUE))

# Create a data frame with a sequence of values for variable of interest
# and separate rows for each level of `Object_Type`
new_data <- expand.grid(
  IP = seq(min(fall_df$IP, na.rm = TRUE), 
           max(all_df$IP, na.rm = TRUE), 
           length.out = 100),
  BS = mean(fall_df$BS, na.rm = TRUE),
  SR = mean(fall_df$SR, na.rm = TRUE),
  Object_Type = c("None", "Vegetation", "Non-veg wood or shell")
)

# Get predictions with standard errors
predictions <- predict(f.mod.fall.is, 
                       new_data, 
                       type = "link", 
                       se.fit = TRUE) 

# Calculate the predicted probabilities and confidence intervals
new_data$predicted_prob <- plogis(predictions$fit)
new_data$lower_ci <- plogis(predictions$fit - 1.96 * predictions$se.fit)
new_data$upper_ci <- plogis(predictions$fit + 1.96 * predictions$se.fit)

# Plot using ggplot2 with confidence intervals as a ribbon
IP_cover <- ggplot(new_data, aes(x = IP, y = predicted_prob, color = Object_Type)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = Object_Type), color = NA, alpha = 0.3) +
  labs(
    title = "% Cover Ice Plant in 2m Radius",
    x = "Percent Cover",
    y = "Probability of nest",
    color = "Object Type",  # Change legend title for line color
    fill = "Object Type"    # Change legend title for ribbon fill
  ) +
  theme_classic() +
  xlim(range$min_IP, range$max_IP)

# Display the plot
print(IP_cover) 



# 3.d. Facet Plot for Fall Data ----

# Remove legends and apply correct titles/axis labels
TwoVegp <- TwoVeg_cover + theme(legend.position = "none",
                                plot.title = element_text(size = 11),
                                axis.title.x = element_text(size = 10), 
                                axis.title.y = element_text(size = 10)) +
  scale_color_grey() + scale_fill_grey() +
  font_theme + 
  labs(title = "(a) Vegetative cover in 2 m radius", x = "Percent cover in 2 m radius", 
       y = "Probability of nest")

Viewsmp <- View_sm + theme(legend.position = "none") +
  scale_color_grey() + scale_fill_grey() +
  font_theme + 
  labs(title = "(b) Viewshed", 
       x = "Proportion open view from nest")

PSp <- PS_cover + theme(legend.position = "none", 
                        plot.title = element_text(size = 11), 
                        axis.title.x = element_text(size = 10), 
                        axis.title.y = element_text(size = 10)) +
  scale_color_grey() + scale_fill_grey() +
  labs(title = "(a) Pioneer stabilizers", x = "Percent cover in 2 m radius", 
       y = "Probability of nest")

PBp <- PB_cover + theme(legend.position = "none") +
  scale_color_grey() + scale_fill_grey() +
  font_theme + 
  labs(title = "(b) Pioneer builders", x = "Percent cover in 2 m radius", 
       y = "Probability of nest")

Viewfgp <- View_fg + theme(legend.position = "none") +
  scale_color_grey() + scale_fill_grey() +
  font_theme + 
  labs(title = "(c) Viewshed", 
       x = "Proportion open view from nest")

BSp <- BS_cover + theme(legend.position = "none", 
                        plot.title = element_text(size = 11), 
                        axis.title.x = element_text(size = 10), 
                        axis.title.y = element_text(size = 10)) +
  scale_color_grey() + scale_fill_grey() +
  font_theme + 
  labs(title = "(a) Beach bur sage", x = "Percent cover in 2 m radius", 
       y = "Probability of nest")

SRp <- SR_cover + theme(legend.position = "none") +
  scale_color_grey() + scale_fill_grey() +
  font_theme + 
  labs(title = "(b) Sea rocket", x = "Percent cover in 2 m radius", 
       y = "Probability of nest")

IPp <- IP_cover + theme(legend.position = "none") +
  scale_color_grey() + scale_fill_grey() +
  font_theme + 
  labs(title = "(c) Dead ice plant", x = "Percent cover in 2 m radius", 
       y = "Probability of nest")


# Extract single legend before removing it from the plot
legend <- get_legend(
  BSp +
    theme(
      legend.position = "right",
      legend.direction = "horizontal",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 10)
    ) +
    guides(color = guide_legend(nrow=1))
)

# Make an empty plot for spacing
empty_plot <- ggplot() + theme_void()

# Create label plots
label_top    <- make_row_label("Summarized Microhabitat")
label_middle <- make_row_label("Functional Groups")
label_bottom <- make_row_label("Individual Species")

# Example top row: label + plots
top_row <- plot_grid(
  label_top,
  plot_grid(TwoVegp, Halfp, Viewsmp, nrow = 1, rel_widths = c(1,1,1)),
  ncol = 2,
  rel_widths = c(0.05, 0.95)  # narrow label, wide plots
)

middle_row <- plot_grid(
  label_middle,
  plot_grid(PSp, PBp, Viewfgp, nrow = 1, rel_widths = c(1,1,1)),
  ncol = 2,
  rel_widths = c(0.05, 0.95)
)

bottom_row <- plot_grid(
  label_bottom,
  plot_grid(BSp, SRp, IPp, nrow = 1, rel_widths = c(1,1,1)),
  ncol = 2,
  rel_widths = c(0.05, 0.95)
)

# Combine all rows with legend
legend_row <- plot_grid(legend, nrow = 1)
facet_falldata_3groups <- plot_grid(top_row, middle_row, bottom_row, legend_row, ncol = 1, rel_heights = c(1,1,1,0.2))

# Display
print(facet_falldata_3groups)

ggsave("fig/facet_falldata_3groups.jpg", plot = facet_falldata_3groups, width = 170, height = 170, units = "mm", dpi = 300) #adjust size as needed