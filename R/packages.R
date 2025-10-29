# load packages

library(tidyverse)
library(vegan)
library(grid)
library(cowplot)
library(openxlsx) #check if ok to remove since data will be included in project
library(buildmer)
library(glmmTMB)
library(performance)
library(DHARMa)
library(gt) # this was for building a table, not for analysis, so probably can remove too
library(purrr)
library(car)
library(broom.mixed)
library(stringr)

conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")