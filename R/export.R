# export_results: explicit inputs, with figures and tables written through config.
export_results <- function(data, variograms, surface, green, facilities, sensitivity, config) {
  svy21 <- config$crs
  save_fig <- config$save_fig
  save_tab <- config$save_tab
  obs <- data$obs
  summary_stats <- data$summary_stats
  spacing <- data$spacing
  coord_check_m <- data$coord_check_m
  land_km2 <- data$land_km2
  sample_var <- variograms$sample_var
  empirical <- variograms$empirical
  tuning <- variograms$tuning
  fits <- variograms$fits
  candidates <- variograms$candidates
  control_cv <- variograms$control_cv
  selected <- variograms$selected
  vfit <- variograms$vfit
  v_sel <- variograms$v_sel
  selected_name <- variograms$selected_name
  cells <- surface$cells
  cell_pts <- surface$cell_pts
  sanity <- surface$sanity
  exceedance <- surface$exceedance
  zone_table <- surface$zone_table
  cv_sel <- surface$cv_sel
  zone_block_total <- surface$zone_block_total
  lattice <- surface$lattice
  zones <- surface$zones
  spearman <- green$spearman
  greenspace_model <- green$greenspace_model
  glm_pois <- green$glm_pois
  moran_table <- green$moran_table
  band_m <- green$band_m
  v_resid <- green$v_resid
  resid_var <- green$resid_var
  week5_tests <- green$week5_tests
  cv_moran <- green$cv_moran
  centers <- facilities$centers
  facility_summary <- facilities$facility_summary
  placement <- facilities$placement
  sensitivity <- sensitivity$sensitivity

  # ---- 7. Save everything the report needs ------------------------------------
  results <- list(
    sites = st_drop_geometry(obs), summary_stats = summary_stats, spacing = spacing,
    coord_check_m = coord_check_m, land_km2 = land_km2, sample_var = sample_var,
    empirical = select(empirical, -v), tuning = select(tuning, -model_fit),
    fits = select(fits, -model_fit), candidates = select(candidates, -model_fit),
    control_cv = control_cv, selected = select(selected, -model_fit), vfit = vfit,
    v_sel = as.data.frame(v_sel), sanity = sanity, exceedance = exceedance,
    zone_table = zone_table, spearman = spearman, greenspace_model = greenspace_model,
    glm_coef = summary(glm_pois)$coefficients, moran = moran_table, band_km = band_m / 1000,
    v_resid = as.data.frame(v_resid), resid_var = resid_var, week5_tests = week5_tests,
    cv_moran = c(I = unname(cv_moran$statistic), p = cv_moran$p.value),
    zone_block_total = zone_block_total, lattice_points = nrow(lattice),
    centers = st_drop_geometry(centers), facility_summary = facility_summary,
    pair_ranking = head(placement$ranking, 10), sensitivity = sensitivity,
    versions = sapply(c("sf", "gstat", "spdep", "dplyr", "ggplot2"),
                      \(p) as.character(packageVersion(p))),
    r_version = R.version.string, gdal = sf_extSoftVersion()[["GDAL"]],
    run_time = format(Sys.time(), "%Y-%m-%d %H:%M %Z"))
  saveRDS(results, file.path(config$output_dir, "results.rds"))
  save_tab(bind_cols(st_drop_geometry(cells), as.data.frame(st_coordinates(cell_pts))) |>
             select(cell_id, px = X, py = Y, land_km2, in_hull, nearest_site_km, pred, var,
                    risk_class, p_exceed3), "kriging_grid.csv")
  writeLines(capture.output(sessionInfo()), file.path(config$output_dir, "sessionInfo.txt"))
  message("Done. Selected: ", selected_name, " model; ", nrow(zones), " top-25% zone(s); centers: ",
          paste(centers$name, collapse = " and "))

  invisible(results)
}
