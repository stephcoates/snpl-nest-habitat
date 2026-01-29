# RUN NMDS ANALYSIS
# Compare nests and random points in breeding season and post-breeding fall

# load packages and functions
source('R/packages.R')
source('R/functions_NMDS.R')

#alldata <- read_csv('data/SNPL_nest_habitat.csv')
alldata <- read_csv('data/SNPL_nest_habitat_clean.csv')

# NMDS data preparation ----
# add group type prefixes based on location type and season
alldata <- alldata %>% 
  mutate(ID = case_when(
    Season == "Nesting" & Type == "Nest" ~ paste0("SN", ID), #summer nests
    Season == "Nesting" & Type == "Random Point" ~ paste0("SP", ID), #summer random points
    Season == "Fall" & Type == "Nest" ~ paste0("FN", ID), #fall nests
    Season == "Fall" & Type == "Random Point" ~ paste0("FP", ID) #fall random points
  )) %>% 
  filter(Unspec!="In River") #removed points that fell in river in the fall

alldata <- alldata %>% mutate(Viewshed=Viewshed/360) #scaling viewshed to match other variables

## NMDS analaysis ----
# Three covariate sets; each groups individual plant species differently
# 1) summarized microhabitat
# 2) dune plant functional groups
# 3) individual plant species


## 1. Summarized Microhabitat Ordination 
# comparing point location type and season by general percent cover 
# categories such as total cover in a 2m radius, cover of dead vegetation

# Remove outliers and unnecessary columns, then set row names
nmdsdata_nv <- alldata %>%
  select(-TwoVeg) %>%  # Remove unwanted columns, TwoVeg is essentially a duplicate of Two and LiveVeg
  rename(DeadVeg = DeadVegManual, LiveVeg = Live) %>% # Rename to match manuscript table
  column_to_rownames(var = "ID") %>%  # Set 'ID' as row names
  select(Half:DeadVeg, Viewshed,NonVeg)
# Create a grouping variable
groups_nv <- case_when(
  grepl("^SN", rownames(nmdsdata_nv)) ~ "Breeding Season Nest Site",
  grepl("^FN", rownames(nmdsdata_nv)) ~ "Fall Nest Site",
  grepl("^SP", rownames(nmdsdata_nv)) ~ "Breeding Season Random Point",
  grepl("^FP", rownames(nmdsdata_nv)) ~ "Fall Random Point",
  TRUE ~ NA_character_
)

# Convert groups to a factor
groups_nv <- factor(groups_nv)
summary(groups_nv)

# run NMDS on data with scaled viewshed
set.seed(12)
myresult_nv <- nmds(nmdsdata_nv,"NMDS: All Groups")
myresult_nv
stressplot(myresult_nv$nmds_result)

# overall adonis results with scaled data
myadonis_result_nv <- adonis2(nmdsdata_nv ~ groups_nv, data = data.frame(groups = groups_nv), permutations = 999)
myadonis_result_nv
# F=4.7653 and P=0.001
# STC: F = 3.3868 and p = 0.006 w/ Two and TwoVeg
# STC: F = 3.732 and p = 0.002 w/ just Two

# pairwise permanova result
mypairwise_results2_nv <- pairwise.adonis(nmdsdata_nv,groups_nv)
mypairwise_results2_nv
# sig adjusted results:
# summer nest and summer point p=0.024 #GCU would use this to compare nests and random points
# summer nest and fall point p=0.006 
# fall nest and summer point p=0.006
# fall nest and fall point p=0.012 #GCU would use this to compare nests and random points

# STC: sig adjusted results after adding Two and TwoVeg back in:
# summer nest and fall point p = 0.036 w/ Two and TwoVeg
# STC: sig (and close) adjusted results after adding just Two back in:
# summer nests and fall point p = 0.018
# fall nest and fall point p = 0.054
# fall nest and summer point p = 0.066

# beta dispersion test
d_nv <- vegdist(nmdsdata_nv, "bray")
betadis_nv <- betadisper(d_nv,group = groups_nv) #negative eigenvalues are generated, I think we can ignore this
#betadis_nv <- betadisper(d_nv,group = groups_nv, add = TRUE) #corrects for negative values are results are the same
disp_result_nv <- permutest(betadis_nv) #tests if there are differences (non-homogeneous dispersion between groups)
disp_result_nv
# F=3.6813 and P=0.014 with 999 permutations including negative eigenvalues
# F=4.0987 and P=0.009 with 999 permutations, using add=TRUE to correct negative eigenvalues

# STC: results after adding Two and TwoVeg back in:
# F = 4.5748 and p = 0.007 with 999 permutations including negative eigenvalues
# F = 5.3598 and p = 0.003 with 999 permutations, using add= TRUE to correct negative eigenvalues

# STC: results after adding just Two back in:
# F = 4.4644 and p = 0.008 with 999 permutations including negative eigenvalues
# F = 4.9661 and p = 0.003 with 999 permutations, using add= TRUE to correct negative eigenvalues

TukeyHSD(betadis_nv)
# Significant or close to pairwise results for dispersion:
# summer nest and fall points: P=0.015
# summer nest and summer points: 0.05 #GCU would use this to compare nests and random points
# Significant or close to pairwise results for dispersion WITH Corrected Eigenvalues:
# summer nest and fall points: P=0.065
# summer nest and summer points: 0.017

# STC: Significant or close to pairwise results for dispersion after adding Two and TwoVeg back in:
# summer nest and fall points: P=0.00506
# summer nest and summer points: 0.03859

# STC: Significant or close to pairwise results for dispersion after adding just Two back in:
# summer nest and fall points: P=0.005901
# summer nest and summer points: 0.036112

# STC: Significant or close to pairwise results for dispersion WITH Corrected Eigenvalues after adding Two and TwoVeg back in:
# summer nest and fall points: P=0.019
# summer nest and summer points: p=0.009 #GCU would use this to compare nests and random points
# fall nest and summer random point: p=0.04356
# fall random point and fall nest: p = 0.0782

# STC: Significant or close to pairwise results for dispersion WITH Corrected Eigenvalues after adding just Two back in:
# summer nest and summer point: p = 0.0104458
# summer nest and fall point: p = 0.0294353
# fall nest and summper point: p = 0.0545488

# Print the NMDS plot with permanova p value
mynmds_plot_with_pvalue_nv <- myresult_nv$plot +
  labs(
    title = "NMDS: Landscape Variability",  # Updated title
    subtitle = paste0(
      #"PERMANOVA p = ", myadonis_result_nv$`Pr(>F)`[1], "; ",
      "PERMUTEST p = ", disp_result_nv$tab$`Pr(>F)`[1]
    )
  ) +  
  theme(plot.title = element_text(hjust = 0.5),  # Center title
        plot.subtitle = element_text(hjust = 0.5))  # Center subtitle
mynmds_plot_with_pvalue_nv

## 2. Functional Groups Ordination
# comparing point location type and season based on species sorted into 
# dune plant functional groups

# Remove outliers and unnecessary columns, then set row names
nmdsdata_fg <- alldata %>%
  column_to_rownames(var = "ID") %>%  # Set 'ID' as row names
  select(Pioneer_Builders:NonVeg) # Select relevant columns (just functional groups)

#Remove rows that are all 0 - don't have any cover at all so can't be in ordination
rowSums(nmdsdata_fg) %>% sort()
nmdsdata_fg$sum <- rowSums(nmdsdata_fg)
nmdsdata_fg <- nmdsdata_fg %>% filter(sum!=0) %>% select(-sum) #gets rid of 9 points

# Create a grouping variable
groups_fg <- case_when(
  grepl("^SN", rownames(nmdsdata_fg)) ~ "Breeding Season Nest Site",
  grepl("^FN", rownames(nmdsdata_fg)) ~ "Fall Nest Site",
  grepl("^SP", rownames(nmdsdata_fg)) ~ "Breeding Season Random Point",
  grepl("^FP", rownames(nmdsdata_fg)) ~ "Fall Random Point",
  TRUE ~ NA_character_
)

# Convert groups to a factor
groups_fg <- factor(groups_fg)
summary(groups_fg) #Count up number of locations left within each group

# run NMDS on data with scaled viewshed
set.seed(12)
myresult_fg <- nmds(nmdsdata_fg,"NMDS: Functional Groups")
myresult_fg
stressplot(myresult_fg$nmds_result)

# overall adonis results with my scaled data
myadonis_result_fg <- adonis2(nmdsdata_fg ~ groups_fg, data = data.frame(groups = groups_fg), permutations = 999)
myadonis_result_fg
# F=3.0706 and P=0.001
# pairwise permanova result
mypairwise_results2_fg <- pairwise.adonis(nmdsdata_fg,groups_fg)
mypairwise_results2_fg
# summer nest and fall point are p=0.006 adjusted
#  Fall Nest vs Fall Point are p =0.006

# beta dispersion test
d_fg <- vegdist(nmdsdata_fg, "bray")
betadis_fg <- betadisper(d_fg,group = groups_fg) 
disp_result_fg <- permutest(betadis_fg) #tests if there are differences (non-homogeneous dispersion between groups)
disp_result_fg
# F=4.0546 and P=0.008 with 999 permutations including negative eigenvalues
TukeyHSD(betadis_fg)
# Significant or close to pairwise results for dispersion:
# summer point and fall nest = 0.016
# summer point and summer nest = 0.0159

# Print the NMDS plot with permanova p value
mynmds_plot_with_pvalue_fg <- myresult_fg$plot +
  labs(subtitle = paste0("PERMUTEST p = ",disp_result_fg$tab$`Pr(>F)`[1])) +
  #"PERMANOVA p = ",myadonis_result_fg$`Pr(>F)`[1],"; ",) +  # Add p-value as subtitle
  theme(plot.subtitle = element_text(hjust = 0.5))  # Center the subtitle
mynmds_plot_with_pvalue_fg


## 3. Individual Species Ordination
# comparing point location type and season based on individual plant species, 
# not sorted into group

# Remove outliers and unnecessary columns, then set row names
nmdsdata_sp0 <- alldata %>%
  column_to_rownames(var = "ID") %>%  # Set 'ID' as row names
  select(SR:MI) # Select relevant columns, just individual species plus kelp, velella, etc.
# If you want to run this with ONLY vegetation species cover, use this instead:
# nmdsdata_sp0 <- alldata %>%
#   select(-Two, -Live, -Dead, -DeadVeg) %>%  # Remove unwanted columns
#   column_to_rownames(var = "ID") %>%  # Set 'ID' as row names
#   select(Sand:PW) %>% select(-Sand) #Don't include kelp, velella, woody, shell, or misc
# If Sand is not used, then you will need to exclude areas with no other cover:
rowSums(nmdsdata_sp0) %>% sort()
nmdsdata_sp0$sum <- rowSums(nmdsdata_sp0)
nmdsdata_sp0 <- nmdsdata_sp0 %>% filter(sum!=0) %>% select(-sum) #gets rid of 9 points

# Create a grouping variable
groups_sp0 <- case_when(
  grepl("^SN", rownames(nmdsdata_sp0)) ~ "Breeding Season Nest Site",
  grepl("^FN", rownames(nmdsdata_sp0)) ~ "Fall Nest Site",
  grepl("^SP", rownames(nmdsdata_sp0)) ~ "Breeding Season Random Point",
  grepl("^FP", rownames(nmdsdata_sp0)) ~ "Fall Random Point",
  TRUE ~ NA_character_
)

# Convert groups to a factor
groups_sp0 <- factor(groups_sp0)
summary(groups_sp0) #count how many locations in each group

# run NMDS
set.seed(12)
myresult_sp0 <- nmds(nmdsdata_sp0,"NMDS: All Groups", labelpoints = FALSE)
myresult_sp0
stressplot(myresult_sp0$nmds_result)

# overall adonis results
myadonis_result_sp0 <- adonis2(nmdsdata_sp0 ~ groups_sp0, data = data.frame(groups = groups_sp0), permutations = 999)
myadonis_result_sp0
# F=2.6898 and P=0.001
# pairwise permanova result
mypairwise_results2_sp0 <- pairwise.adonis(nmdsdata_sp0,groups_sp0)
mypairwise_results2_sp0
# summer nest and fall point are p=0.006
# fall nest and summer point p = 0.012
# fall nest and fall point p = 0.018


# beta dispersion test
d_sp0 <- vegdist(nmdsdata_sp0, "bray")
betadis_sp0 <- betadisper(d_sp0,group = groups_sp0) #negative eigenvalues NOT generated
disp_result_sp0 <- permutest(betadis_sp0) #tests if there are differences (non-homogeneous dispersion between groups)
disp_result_sp0
# F=5.8907 and P=0.001 with 999 permutations
TukeyHSD(betadis_sp0)
# Significant or close to pairwise results for dispersion:
# summer nest and summer points: 0.022
# summer point and fall nest --> p = 0.0003

# Print the NMDS plot with permanova p value
mynmds_plot_with_pvalue_sp0 <- myresult_sp0$plot +
  labs(title = "NMDS: Individual Cover Type",
       subtitle = paste0(
         #"PERMANOVA p = ",myadonis_result_sp0$`Pr(>F)`[1],"; ",
         "PERMUTEST p = ",disp_result_sp0$tab$`Pr(>F)`[1])) +  # Add p-value as subtitle
  theme(plot.subtitle = element_text(hjust = 0.5))  # Center the subtitle
mynmds_plot_with_pvalue_sp0

