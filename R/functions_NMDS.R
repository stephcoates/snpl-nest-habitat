# NMDS function without p-value/adonis analysis ----
nmds <- function(data, title, k = 2, labelpoints = FALSE, max_covariate_labels = 6) {
  
  subset_data <- switch(
    title,
    "NMDS: Summer Nest vs Fall Nest"     = data[grepl("^SN|^FN", rownames(data)), ],
    "NMDS: Summer Point vs Fall Point"  = data[grepl("^SP|^FP", rownames(data)), ],
    "NMDS: Summer Nest vs Summer Point" = data[grepl("^SN|^SP", rownames(data)), ],
    "NMDS: Fall Nest vs Fall Point"     = data[grepl("^FN|^FP", rownames(data)), ],
    "NMDS: All Groups"                  = data,
    "NMDS: Functional Groups"           = data,
    "NMDS: Landscape Variability"       = data,
    stop("title not recognized in switch()")
  )
  
  group_levels <- c(
    "Breeding Season Nest Site",
    "Fall Nest Site",
    "Breeding Season Random Point",
    "Fall Random Point"
  )
  
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
  
  rn <- trimws(rownames(subset_data))
  
  groups <- dplyr::case_when(
    grepl("^SN", rn) ~ "Breeding Season Nest Site",
    grepl("^FN", rn) ~ "Fall Nest Site",
    grepl("^SP", rn) ~ "Breeding Season Random Point",
    grepl("^FP", rn) ~ "Fall Random Point",
    TRUE ~ NA_character_
  )
  
  nmds_result <- vegan::metaMDS(
    subset_data,
    distance = "bray",
    k = k,
    trymax = 1000
  )
  
  nmds_scores <- as.data.frame(vegan::scores(nmds_result, display = "sites"))
  rownames(nmds_scores) <- rownames(subset_data)
  nmds_scores$Group <- factor(groups, levels = group_levels)
  
  nmds_scores <- nmds_scores %>%
    dplyr::filter(!is.na(Group))
  
  species_scores <- as.data.frame(vegan::scores(nmds_result, display = "species"))
  
  species_scores <- species_scores %>%
    dplyr::mutate(
      label = rownames(species_scores),
      NMDS1_lab = NMDS1 / 2,
      NMDS2_lab = NMDS2 / 2,
      dist_from_origin = sqrt(NMDS1_lab^2 + NMDS2_lab^2)
    ) %>%
    dplyr::arrange(dplyr::desc(dist_from_origin))
  
  species_scores_labeled <- species_scores %>%
    dplyr::slice_head(n = max_covariate_labels)
  
  message("Covariate labels plotted for ", title, ": ", nrow(species_scores_labeled))
  
  find_hull <- function(df) df[chull(df$NMDS1, df$NMDS2), ]
  
  hulls <- nmds_scores %>%
    dplyr::group_by(Group) %>%
    dplyr::do(find_hull(.))
  
  grob <- grid::grobTree(
    grid::textGrob(
      paste0("Stress = ", round(nmds_result$stress, 3)),
      x = 0.9,
      y = 0.9,
      hjust = 1,
      gp = grid::gpar(col = "black", fontsize = 9)
    )
  )
  
  nmds_plot <- ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = hulls,
      ggplot2::aes(x = NMDS1, y = NMDS2, fill = Group),
      alpha = 0.2,
      color = NA
    ) +
    ggplot2::geom_point(
      data = nmds_scores,
      ggplot2::aes(x = NMDS1, y = NMDS2, color = Group, shape = Group),
      size = 3
    ) +
    ggplot2::geom_text(
      data = species_scores_labeled,
      ggplot2::aes(
        x = NMDS1_lab,
        y = NMDS2_lab,
        label = label
      ),
      color = "black",
      size = 2.6,
      check_overlap = TRUE
    ) +
    ggplot2::annotation_custom(grob) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = title, x = "NMDS1", y = "NMDS2") +
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
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      plot.margin = ggplot2::margin(10, 20, 10, 10)
    )
  
  if (isTRUE(labelpoints)) {
    nmds_plot <- nmds_plot +
      ggplot2::geom_text(
        data = nmds_scores,
        ggplot2::aes(
          x = NMDS1,
          y = NMDS2,
          label = rownames(nmds_scores),
          color = Group
        ),
        show.legend = FALSE,
        check_overlap = TRUE
      )
  }
  
  return(list(
    plot = nmds_plot,
    nmds_result = nmds_result,
    species_scores = species_scores,
    species_scores_labeled = species_scores_labeled
  ))
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
