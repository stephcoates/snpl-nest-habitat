# NMDS facet plot 

# load packages and functions
source('R/packages.R')
source('R/functions_NMDS.R')

# Fig 2. NMDS plot of nests/random and breeding/fall points.

# Prepare plots with left-aligned subtitles and no legends
p_a <- mynmds_plot_with_pvalue_nv +
  labs(title = NULL, subtitle = "(a) Summarized microhabitat features") +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.subtitle = element_text(size = 11, hjust = 0, margin = margin(t=0, r=0, b=5, l=0)),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
  )

p_b <- mynmds_plot_with_pvalue_fg +
  labs(title = NULL, subtitle = "(b) Functional groups") +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.subtitle = element_text(size = 11, hjust = 0, margin = margin(t=0, r=0, b=5, l=0)),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
  )

p_c <- mynmds_plot_with_pvalue_sp0 +
  labs(title = NULL, subtitle = "(c) Individual cover type") +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.subtitle = element_text(size = 11, hjust = 0, margin = margin(t=0, r=0, b=5, l=0)),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
  )


# Extract legend (small text)
legend <- get_legend(
  mynmds_plot_with_pvalue_fg + 
    theme(legend.position = "right", legend.text = element_text(size = 10))
)

# Arrange in 2x2 grid with labels a, b, c
combined_2x2 <- plot_grid(
  p_a, p_b,
  p_c, legend,
  ncol = 2, nrow = 2,
  rel_widths = c(1, 1),
  rel_heights = c(1, 1),
  align = "hv",
  axis = "tblr"
)


# Draw final combined plot
final_plot <- ggdraw() +
  draw_plot(combined_2x2, y = 0, height = 0.95)

print(final_plot)

# Save to size formatted for Wader Study
ggsave("fig/NMDS_facet.png", final_plot, width = 170, height = 170, units = "mm", dpi = 300) 
