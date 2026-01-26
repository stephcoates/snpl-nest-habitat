# NMDS function without p-value/adonis analysis ----
nmds <- function(data, title, k = 2, labelpoints = FALSE) {
  
  # Determine subset based on the title argument
  subset_data <- switch(
    title,
    "NMDS: Summer Nest vs Fall Nest"    = data[grepl("^SN|^FN", rownames(data)), ],
    "NMDS: Summer Point vs Fall Point" = data[grepl("^SP|^FP", rownames(data)), ],
    "NMDS: Summer Nest vs Summer Point"= data[grepl("^SN|^SP", rownames(data)), ],
    "NMDS: Fall Nest vs Fall Point"    = data[grepl("^FN|^FP", rownames(data)), ],
    "NMDS: All Groups"                 = data,
    "NMDS: Functional Groups"          = data,
    stop("title not recognized in switch()")
  )
  
  # Define the canonical group levels ONCE (order matters)
  group_levels <- c(
    "Breeding Season Nest Site",
    "Fall Nest Site",
    "Breeding Season Random Point",
    "Fall Random Point"
  )
  
  # Palette and shapes keyed to those exact labels
  pal <- c(
    "Breeding Season Nest Site"     = "#0072B2",
    "Fall Nest Site"               = "#009E73",
    "Breeding Season Random Point" = "#D55E00",
    "Fall Random Point"            = "#E69F00"
  )
  
  shp <- c(
    "Breeding Season Nest Site"     = 16,
    "Fall Nest Site"               = 10,
    "Breeding Season Random Point" = 15,
    "Fall Random Point"            = 12
  )
  
  # Add Group labels for plotting
  rn <- rownames(subset_data)
  
  # defensive: trim whitespace just in case
  rn <- trimws(rn)
  
  groups <- dplyr::case_when(
    grepl("^SN", rn) ~ "Breeding Season Nest Site",
    grepl("^FN", rn) ~ "Fall Nest Site",
    grepl("^SP", rn) ~ "Breeding Season Random Point",
    grepl("^FP", rn) ~ "Fall Random Point",
    TRUE ~ NA_character_
  )
  
  # Run NMDS
  nmds_result <- vegan::metaMDS(subset_data, distance = "bray", k = k, trymax = 1000)
  
  if (is.null(nmds_result$stress)) {
    stop("NMDS failed to converge.")
  }
  
  # Extract NMDS site scores
  nmds_scores <- as.data.frame(vegan::scores(nmds_result, display = "sites"))
  rownames(nmds_scores) <- rownames(subset_data)
  
  # Force Group to a factor with fixed levels
  nmds_scores$Group <- factor(groups, levels = group_levels)
  
  # If you have NA groups, show them clearly (and you can choose to drop them)
  if (any(is.na(nmds_scores$Group))) {
    message("Missing Group labels detected. These rows will NOT be plotted:\n",
            paste(rownames(nmds_scores)[is.na(nmds_scores$Group)], collapse = ", "))
    nmds_scores <- dplyr::filter(nmds_scores, !is.na(Group))
  }
  
  # Extract species scores for plotting
  species_scores <- as.data.frame(vegan::scores(nmds_result, display = "species"))
  
  # Generate convex hulls for each group
  find_hull <- function(df) df[chull(df$NMDS1, df$NMDS2), ]
  
  hulls <- nmds_scores %>%
    dplyr::group_by(Group) %>%
    dplyr::do(find_hull(.))
  
  # Stress grob
  grob <- grid::grobTree(
    grid::textGrob(
      paste0("Stress = ", round(nmds_result$stress, 3)),
      x = 0.9, y = 0.9, hjust = 1,
      gp = grid::gpar(col = "black", fontsize = 9)
    )
  )
  
  # Plot NMDS results
  nmds_plot <- ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = hulls,
      ggplot2::aes(x = NMDS1, y = NMDS2, fill = Group),
      alpha = 0.2, color = NA
    ) +
    ggplot2::geom_point(
      data = nmds_scores,
      ggplot2::aes(x = NMDS1, y = NMDS2, color = Group, shape = Group),
      size = 3
    ) +
    ggplot2::geom_text(
      data = species_scores,
      ggplot2::aes(x = NMDS1 / 2, y = NMDS2 / 2, label = rownames(species_scores)),
      color = "black", size = 3, vjust = -0.5
    ) +
    ggplot2::annotation_custom(grob) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = title, x = "NMDS1", y = "NMDS2") +
    # LOCK mapping across plots + legend extraction
    ggplot2::scale_shape_manual(
      name = "Group",
      values = shp,
      limits = group_levels,
      breaks = group_levels,
      drop = FALSE
    ) +
    ggplot2::scale_color_manual(
      name = "Group",
      values = pal,
      limits = group_levels,
      breaks = group_levels,
      drop = FALSE,
      na.translate = FALSE
    ) +
    ggplot2::scale_fill_manual(
      name = "Group",
      values = pal,
      limits = group_levels,
      breaks = group_levels,
      drop = FALSE,
      na.translate = FALSE
    ) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  
  if (isTRUE(labelpoints)) {
    nmds_plot <- nmds_plot +
      ggplot2::geom_text(
        data = nmds_scores,
        ggplot2::aes(x = NMDS1, y = NMDS2, label = rownames(nmds_scores), color = Group),
        show.legend = FALSE
      )
  }
  
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
