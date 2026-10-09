#!/usr/bin/env Rscript
# Run from the repository root: Rscript run_analysis.R
# Modules exchange explicit named lists; no intermediate objects enter .GlobalEnv.
main <- function() {
  if (!file.exists("R/setup.R")) stop("Run from the kriging-incidence repository root.")
  modules <- c("setup", "variogram_helpers", "zone_helpers", "facility_helpers",
               "sensitivity_helpers", "data", "variograms", "kriging", "green_space",
               "facilities", "sensitivity", "export")
  for (module in modules) source(file.path("R", paste0(module, ".R")), local = TRUE)
  load_dependencies()
  config <- pipeline_config()
  set.seed(config$seed)
  message("Preparing data...")
  data <- prepare_data(config)
  message("Fitting and validating variograms...")
  variograms <- fit_variograms(data, config)
  message("Kriging and conditional simulation...")
  surface <- interpolate_counts(data, variograms, config)
  message("Green-space association and spatial diagnostics...")
  green <- analyze_green_space(data, variograms, surface, config)
  message("Optimizing center placement...")
  facilities <- place_centers(data, surface, config)
  message("Checking sensitivity...")
  sensitivity <- analyze_sensitivity(data, variograms, surface, facilities, config)
  export_results(data, variograms, surface, green, facilities, sensitivity, config)
}
main()
