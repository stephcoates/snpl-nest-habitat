# NMDS facet plot 

# load packages and functions
source('R/packages.R')
source('R/functions_NMDS.R')

# Fig 2. NMDS plot of nests/random and breeding/fall points.

# Shared theme: no individual axis titles or legends
common_theme <- theme(
  axis.title = element_blank(),
  legend.position = "none",
  plot.subtitle = element_text(
    size = 11,
    hjust = 0,
    margin = margin(t = 0, r = 0, b = 5, l = 0)
  ),
  plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
)

p_a <- mynmds_plot_with_pvalue_nv +
  labs(title = NULL, subtitle = "(a) Summarized microhabitat features") +
  theme_classic() +
  common_theme

p_b <- mynmds_plot_with_pvalue_fg +
  labs(title = NULL, subtitle = "(b) Functional groups") +
  theme_classic() +
  common_theme

p_c <- mynmds_plot_with_pvalue_sp0 +
  labs(title = NULL, subtitle = "(c) Individual cover type") +
  theme_classic() +
  common_theme


# Extract legend
legend <- get_legend(
  mynmds_plot_with_pvalue_fg +
    theme_classic() +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10)
    )
)

# Put legend in a panel that is the same width as plot C
legend_panel <- plot_grid(
  NULL,
  legend,
  NULL,
  ncol = 1,
  rel_heights = c(1, 2, 1)
)

# Top row: plot A + legend panel
top_row <- plot_grid(
  p_a,
  legend_panel,
  ncol = 2,
  rel_widths = c(1, 1),
  align = "h",
  axis = "tb"
)

# Bottom row: plot B + plot C
bottom_row <- plot_grid(
  p_b,
  p_c,
  ncol = 2,
  rel_widths = c(1, 1),
  align = "h",
  axis = "tb"
)

# Combine rows
combined <- plot_grid(
  top_row,
  bottom_row,
  ncol = 1,
  rel_heights = c(1, 1),
  align = "v"
)

# Add shared axis labels
final_plot <- ggdraw() +
  draw_plot(
    combined,
    x = 0.08,
    y = 0.08,
    width = 0.90,
    height = 0.88
  ) +
  draw_label(
    "NMDS1",
    x = 0.53,
    y = 0.025,
    size = 12
  ) +
  draw_label(
    "NMDS2",
    x = 0.025,
    y = 0.52,
    angle = 90,
    size = 12
  )

print(final_plot)

ggsave(
  "fig/NMDS_facet_20July2026.jpg",
  final_plot,
  width = 170,
  height = 170,
  units = "mm",
  dpi = 300
)
