# =============================================================================
# Shared figure style for R panels
# -----------------------------------------------------------------------------
# Palette anchored on sex (female = red family, male = blue family); condition
# encoded by intensity; female in vivo phases in distinct warm tones;
# methylation direction hypo = blue / hyper = red. Compact theme.
#
# Usage:
#   source("code/_helpers/publication_style.R")
#   ... + scale_fill_manual(values = PUB$medip4) + theme_pub()
# =============================================================================

PUB <- list(
  # Four in vitro MeDIP groups (labels with or without "1")
  medip4 = c(
    "Female control" = "#F4A582", "Female CSF1" = "#B2182B", "Female CSF" = "#B2182B",
    "Male control"   = "#92C5DE", "Male CSF1"   = "#2166AC", "Male CSF"   = "#2166AC"
  ),
  # In vivo groups (oestrous phases and males)
  phase = c(
    "Proestrus" = "#E7298A", "Estrus" = "#F4C95D", "Diestrus" = "#E8826B", "Male" = "#2166AC"
  ),
  # Sex
  sex = c("Female" = "#B2182B", "Male" = "#2166AC"),
  # In vitro treatment
  treatment = c("Control" = "#9E9E9E", "CSF1" = "#B2182B", "CSF1+E2" = "#1B9E77",
                "CSF" = "#B2182B", "E2" = "#1B9E77"),
  # Methylation direction
  direction = c("Hypo" = "#2471A3", "Hyper" = "#C0392B",
                "Hypomethylated" = "#2471A3", "Hypermethylated" = "#C0392B"),
  # Generic condition
  condition = c("Control" = "#BDBDBD", "CCI" = "#D1495B", "CCI-IoN" = "#D1495B")
)

# Lighten a colour (mix with white). amount in [0, 1].
pub_lighten <- function(col, amount = 0.45) {
  m <- grDevices::col2rgb(col) / 255
  m <- m + (1 - m) * amount
  grDevices::rgb(m[1, ], m[2, ], m[3, ])
}
# Darken a colour (mix with black).
pub_darken <- function(col, amount = 0.25) {
  m <- grDevices::col2rgb(col) / 255
  m <- m * (1 - amount)
  grDevices::rgb(m[1, ], m[2, ], m[3, ])
}

# Compact publication theme.
theme_pub <- function(base_size = 16, legend = "right", grid = FALSE) {
  th <- ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.title   = ggplot2::element_text(size = base_size + 2, face = "bold", colour = "black"),
      axis.text    = ggplot2::element_text(size = base_size, colour = "black"),
      axis.line    = ggplot2::element_line(linewidth = 0.8, colour = "black"),
      axis.ticks   = ggplot2::element_line(linewidth = 0.8, colour = "black"),
      axis.ticks.length = ggplot2::unit(3, "pt"),
      plot.title   = ggplot2::element_text(size = base_size + 2, face = "bold", hjust = 0.5,
                                           margin = ggplot2::margin(b = 2)),
      legend.title = ggplot2::element_text(size = base_size, face = "bold"),
      legend.text  = ggplot2::element_text(size = base_size - 1),
      legend.key.size = ggplot2::unit(1.0, "lines"),
      legend.position = legend,
      legend.margin   = ggplot2::margin(1, 1, 1, 1),
      legend.box.spacing = ggplot2::unit(2, "pt"),
      strip.text   = ggplot2::element_text(size = base_size, face = "bold"),
      strip.background = ggplot2::element_blank(),
      panel.spacing = ggplot2::unit(0.35, "lines"),
      plot.margin  = ggplot2::margin(4, 6, 3, 3)
    )
  if (grid) th <- th + ggplot2::theme(panel.grid.major.y =
                                        ggplot2::element_line(colour = "grey90", linewidth = 0.4))
  th
}

# Sizes for significance annotations (ggplot size units).
PUB_SIG_SIZE   <- 7      # asterisks (geom_text size)
PUB_SIG_LINEW  <- 0.8    # significance bracket line width
PUB_POINT_SIZE <- 3.2

invisible(NULL)
