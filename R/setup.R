# Package checks, paths, and output writers shared by all stages.
analysis_packages <- c("sf", "gstat", "spdep", "dplyr", "tidyr", "purrr", "ggplot2",
                       "ggspatial", "ggrepel", "patchwork", "units")

load_dependencies <- function() {
  missing <- analysis_packages[!vapply(analysis_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "),
                           ". Run Rscript scripts/install_dependencies.R.", call. = FALSE)
  suppressPackageStartupMessages({
    for (package in analysis_packages) library(package, character.only = TRUE)
  })
}

title_case <- function(x) gsub("\\bMrt\\b", "MRT", tools::toTitleCase(tolower(x)))

pipeline_config <- function(root = ".") {
  root <- normalizePath(root, mustWork = TRUE)
  config <- list(root = root, crs = 3414L, seed = 5105L,
                 sites_file = file.path(root, "data/raw/assignment2.csv"),
                 land_file = file.path(root, "data/cache/singapore_land.geojson"),
                 clubs_file = file.path(root, "data/cache/community_clubs.csv"),
                 exits_file = file.path(root, "data/cache/mrt_exits.geojson"),
                 verification_file = file.path(root, "data/cache/facility_verification.csv"),
                 output_dir = file.path(root, "outputs"))
  inputs <- unlist(config[c("sites_file", "land_file", "clubs_file", "exits_file", "verification_file")])
  if (any(!file.exists(inputs))) stop("Missing input files: ", paste(inputs[!file.exists(inputs)], collapse = ", "))
  for (folder in c("figures", "tables")) {
    dir.create(file.path(config$output_dir, folder), recursive = TRUE, showWarnings = FALSE)
  }
  output_dir <- config$output_dir
  config$save_fig <- function(plot, file, width, height) {
    ggplot2::ggsave(file.path(output_dir, "figures", file), plot, width = width,
                   height = height, dpi = 200, bg = "white")
  }
  config$save_tab <- function(table, file) {
    write.csv(sf::st_drop_geometry(table), file.path(output_dir, "tables", file), row.names = FALSE)
  }
  config
}
