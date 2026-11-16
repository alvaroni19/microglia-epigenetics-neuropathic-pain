# Final-figure save controls for the CAM reproducibility bundle.
# By default, scripts write only figures listed in metadata/final_figure_allowlist.txt.

.cam_find_bundle_root <- function(start = getwd()) {
  start <- normalizePath(start, winslash = "/", mustWork = FALSE)
  current <- start
  repeat {
    if (dir.exists(file.path(current, "data")) &&
        dir.exists(file.path(current, "scripts")) &&
        dir.exists(file.path(current, "metadata"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) break
    current <- parent
  }
  stop("Could not locate bundle root from: ", start)
}

if (!exists("BUNDLE_ROOT", inherits = TRUE)) {
  BUNDLE_ROOT <- .cam_find_bundle_root(getwd())
}
if (!exists("OUTPUT_DIR", inherits = TRUE)) {
  OUTPUT_DIR <- file.path(BUNDLE_ROOT, "generated_outputs", basename(getwd()))
}
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

FINAL_FIGURES_ONLY <- !tolower(Sys.getenv("CAM_FINAL_FIGURES_ONLY", "true")) %in%
  c("0", "false", "no", "all")

FINAL_FIGURE_ALLOWLIST_FILE <- file.path(BUNDLE_ROOT, "metadata", "final_figure_allowlist.txt")
FINAL_FIGURE_ALLOWLIST <- if (file.exists(FINAL_FIGURE_ALLOWLIST_FILE)) {
  unique(trimws(readLines(FINAL_FIGURE_ALLOWLIST_FILE, warn = FALSE)))
} else {
  character()
}
FINAL_FIGURE_ALLOWLIST <- FINAL_FIGURE_ALLOWLIST[nzchar(FINAL_FIGURE_ALLOWLIST)]

cam_should_save_figure <- function(filename) {
  if (!FINAL_FIGURES_ONLY) return(TRUE)
  basename(as.character(filename)) %in% FINAL_FIGURE_ALLOWLIST
}

cam_output_figure_path <- function(filename, path = NULL) {
  filename <- as.character(filename)
  if (grepl("^[A-Za-z]:[/\\\\]", filename) || startsWith(filename, "/")) {
    return(filename)
  }
  file.path(OUTPUT_DIR, basename(filename))
}

ggsave <- function(filename, plot = ggplot2::last_plot(), device = NULL, path = NULL, ...) {
  if (!cam_should_save_figure(filename)) {
    message("[skip non-final figure] ", basename(as.character(filename)))
    return(invisible(NULL))
  }
  out <- cam_output_figure_path(filename, path)
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename = out, plot = plot, device = device, path = NULL, ...)
}
