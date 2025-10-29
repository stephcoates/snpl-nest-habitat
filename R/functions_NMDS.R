# NMDS function without p-value/adonis analysis ----
nmds <- function(data, title, k=2, labelpoints=FALSE) {
  
  # Determine subset based on the title argument
  subset_data <- switch(
    title,
    "NMDS: Summer Nest vs Fall Nest" = data[grepl("^SN|^FN", rownames(data)), ],
    "NMDS: Summer Point vs Fall Point" = data[grepl("^SP|^FP", rownames(data)), ],
    "NMDS: Summer Nest vs Summer Point" = data[grepl("^SN|^SP", rownames(data)), ],
    "NMDS: Fall Nest vs Fall Point" = data[grepl("^FN|^FP", rownames(data)), ],
    "NMDS: All Groups" = data, # use the full dataset for All Groups
    "NMDS: Functional Groups" = data  # Use the full dataset for Functional Groups
  )
  
  # Add Group labels for plotting
  groups <- case_when(
    grepl("^SN", rownames(subset_data)) ~ "Breeding Season Nest Site",
    grepl("^FN", rownames(subset_data)) ~ "Fall Nest Site",
    grepl("^SP", rownames(subset_data)) ~ "Breeding Season Random Point",
    grepl("^FP", rownames(subset_data)) ~ "Fall Random Point",
    TRUE ~ NA_character_
  )
  
  # Run NMDS
  nmds_result <- metaMDS(subset_data, distance = "bray", k = k, trymax = 1000)
  
  # Check if NMDS was successful
  if (is.null(nmds_result$stress)) {
    stop("NMDS failed to converge.")
  }
  
  # Extract NMDS site scores
  nmds_scores <- as.data.frame(scores(nmds_result, display = "sites"))
  rownames(nmds_scores) <- rownames(subset_data)
  nmds_scores$Group <- groups  # Add group labels to NMDS scores
  #nmds_scores$NMDS2 <- nmds_scores$NMDS2 * scale_factor
  
  # Check for any missing group labels
  if (sum(is.na(nmds_scores$Group)) > 0) {
    print("There are missing groups in NMDS scores!")
  }
  
  # Extract species scores for plotting
  species_scores <- as.data.frame(scores(nmds_result, display = "species"))
  
  # Generate convex hulls for each group
  find_hull <- function(df) {
    df[chull(df$NMDS1, df$NMDS2), ]
  }
  
  hulls <- nmds_scores %>%
    group_by(Group) %>%
    do(find_hull(.))
  
  #Generate grob to add Stress value to plot
  grob <- grobTree(textGrob(paste0("Stress = ",round(nmds_result$stress,3)),     
                            x = 0.9, y = 0.9, hjust = 1,
                            gp=gpar(col="black", fontsize = 9)))
  
  # Plot NMDS results with convex hulls and species scores
  nmds_plot <- ggplot() +
    geom_polygon(data = hulls, aes(x = NMDS1, y = NMDS2, fill = Group), alpha = 0.2, color = NA) +
    geom_point(data = nmds_scores, aes(x = NMDS1, y = NMDS2, color = Group, shape = Group), size = 3) +
    geom_text(data = species_scores, aes(x = NMDS1/2, y = NMDS2/2, label = rownames(species_scores)), 
              color = "black", size = 3, vjust = -0.5) +
    annotation_custom(grob) +
    theme_minimal() +
    labs(title = title, x = "NMDS1", y = "NMDS2", fontsize = 10) +
    # Shapes
    scale_shape_manual(
      values = c("Breeding Season Nest Site" = 16, "Fall Nest Site" = 10,
                 "Breeding Season Random Point" = 15, "Fall Random Point" = 12)
    ) +
    # Unified color + fill scale
    # Colorblind friendly: c("#0072B2", "#009E73", "#E69F00", "#D55E00")
    scale_color_manual(
      name = "Group",  # same name for legend
      values = c(
        "Breeding Season Nest" = "#0072B2",
        "Fall Nest" = "#009E73",
        "Breeding Season Random Point" = "#D55E00",
        "Fall Random Point" = "#E69F00"
      )
    ) +
    scale_fill_manual(
      name = "Group",  # same name to merge legend
      values = c(
        "Breeding Season Nest Site" = "#0072B2",
        "Fall Nest Site" = "#009E73",
        "Breeding Season Random Point" = "#D55E00",
        "Fall Random Point" = "#E69F00"
      )
    )+
    
    
    theme(plot.title = element_text(hjust = 0.5))
  
  # Add point labels if requested
  if (labelpoints == TRUE) {
    nmds_plot <- nmds_plot +
      geom_text(data = nmds_scores, aes(x = NMDS1, y = NMDS2, color = Group, label = rownames(nmds_scores)))
  }
  
  
  # Return the NMDS plot and results
  return(list(plot = nmds_plot, nmds_result = nmds_result))
}

# Pairwise adjusted permanova function ----
# adjusted pairwise results using external function with Bonferonni adjusted p-value
# we need to adjust the p-value because we are doing multiple tests across the 
# data and have increased our likelihood of finding a false positive (sig) result
# The pairwise.adonis function from the pairwiseAdonis package
# https://github.com/pmartinezarbizu/pairwiseAdonis/tree/master
# Martinez Arbizu, P. (2020). pairwiseAdonis: Pairwise multilevel comparison using adonis. R package version 0.4

pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='bonferroni',reduce=NULL,perm=999)
{
  
  co <- combn(unique(as.character(factors)),2)
  pairs <- c()
  Df <- c()
  SumsOfSqs <- c()
  F.Model <- c()
  R2 <- c()
  p.value <- c()
  
  
  for(elem in 1:ncol(co)){
    if(inherits(x, 'dist')){
      x1=as.matrix(x)[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem])),
                      factors %in% c(as.character(co[1,elem]),as.character(co[2,elem]))]
    }
    
    else  (
      if (sim.function == 'daisy'){
        x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
      }
      else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
    )
    
    x2 = data.frame(Fac = factors[factors %in% c(co[1,elem],co[2,elem])])
    
    ad <- adonis2(x1 ~ Fac, data = x2,
                  permutations = perm);
    pairs <- c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    Df <- c(Df,ad$Df[1])
    SumsOfSqs <- c(SumsOfSqs,ad$SumOfSqs[1])
    F.Model <- c(F.Model,ad$F[1]);
    R2 <- c(R2,ad$R2[1]);
    p.value <- c(p.value,ad$`Pr(>F)`[1])
  }
  p.adjusted <- p.adjust(p.value,method=p.adjust.m)
  
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'.'
  sig[p.adjusted <= 0.01] <-'*'
  sig[p.adjusted <= 0.001] <-'**'
  sig[p.adjusted <= 0.0001] <-'***'
  pairw.res <- data.frame(pairs,Df,SumsOfSqs,F.Model,R2,p.value,p.adjusted,sig)
  
  if(!is.null(reduce)){
    pairw.res <- subset (pairw.res, grepl(reduce,pairs))
    pairw.res$p.adjusted <- p.adjust(pairw.res$p.value,method=p.adjust.m)
    
    sig = c(rep('',length(pairw.res$p.adjusted)))
    sig[pairw.res$p.adjusted <= 0.1] <-'.'
    sig[pairw.res$p.adjusted <= 0.05] <-'*'
    sig[pairw.res$p.adjusted <= 0.01] <-'**'
    sig[pairw.res$p.adjusted <= 0.001] <-'***'
    pairw.res <- data.frame(pairw.res[,1:7],sig)
  }
  class(pairw.res) <- c("pwadonis", "data.frame")
  return(pairw.res)
}


# Method summary
summary.pwadonis = function(object, ...) {
  cat("Result of pairwise.adonis:\n")
  cat("\n")
  print(object, ...)
  cat("\n")
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
}
