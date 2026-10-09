# interpolate_counts: explicit inputs, with figures and tables written through config.
interpolate_counts <- function(data, variograms, config) {
  svy21 <- config$crs
  save_fig <- config$save_fig
  save_tab <- config$save_tab
  sites <- data$sites
  obs <- data$obs
  land <- data$land
  subzones <- data$subzones
  map_base <- data$map_base
  vfit <- variograms$vfit
  selected_name <- variograms$selected_name
  control_model <- variograms$control_model

  # ---- 3. Ordinary kriging on a 1 km grid -------------------------------------

  # 1 km x 1 km squares aligned to whole kilometers, clipped to land.
  origin <- floor(st_bbox(land)[c("xmin", "ymin")] / 1000) * 1000
  squares <- st_sf(cell_id = seq_along(st_make_grid(land, cellsize = 1000, offset = origin)),
                   geometry = st_make_grid(land, cellsize = 1000, offset = origin))
  st_agr(squares) <- "constant"
  cells <- st_intersection(squares, land) |>
    st_collection_extract("POLYGON") |>
    group_by(cell_id) |>
    summarise(.groups = "drop")
  cells$land_km2 <- as.numeric(st_area(cells)) / 1e6

  # Each square is predicted at a point inside its land part (the center for
  # inland squares); this is point support, not a block average.
  cell_pts <- st_sf(geometry = st_point_on_surface(st_geometry(cells)))
  hull <- st_convex_hull(st_union(obs))
  # The zone-analysis domain uses squares overlapping the sampled-site hull;
  # this is a geometric convention, not a validated reliability boundary.
  cells$in_hull <- lengths(st_intersects(cells, hull)) > 0
  cells$nearest_site_km <- apply(units::drop_units(st_distance(cell_pts, obs)), 1, min) / 1000

  ok <- krige(cases ~ 1, obs, newdata = cell_pts, model = vfit, debug.level = 0)
  cells$pred <- ok$var1.pred
  cells$var <- ok$var1.var

  # High-risk zones: top 25% (and top 10%) of predicted values among squares in
  # the hull-overlap domain. Boundary prediction points may lie outside the hull.
  thresholds <- quantile(cells$pred[cells$in_hull], c(0.75, 0.90), names = FALSE)
  cells$risk_class <- case_when(cells$in_hull & cells$pred >= thresholds[2] ~ "Top 10%",
                                cells$in_hull & cells$pred >= thresholds[1] ~ "Top 25%",
                                TRUE ~ "Other")
  high <- filter(cells, risk_class != "Other")

  # Outbreak size. Estimated cases read each predicted value as the expected count
  # per km² around a point, i.e. they assume a site count reflects about 1 km²
  # around the site (the report also rescales this to a larger catchment).
  # Recorded cases are the observed counts at sampled sites inside the zone.

  top25 <- profile_zones(high, "Top-25% zone", obs, subzones, cells, cell_pts)
  top10 <- profile_zones(filter(high, risk_class == "Top 10%"), "Top-10% core", obs, subzones, cells, cell_pts)
  zones <- top25$zones
  zone_table <- bind_rows(top25$table, top10$table)
  save_tab(zone_table, "high_risk_zones.csv")

  # Check: rectangular 1 km x 1 km block predictions at the high-zone points.
  # A similar zone total is a model-based support sensitivity check, not proof of
  # an unbiased total. The rectangles are not clipped to irregular land portions.
  blocks <- krige(cases ~ 1, obs, newdata = cell_pts[match(high$cell_id, cells$cell_id), ],
                  model = vfit, block = c(1000, 1000), debug.level = 0)
  zone_block_total <- sum(blocks$var1.pred * high$land_km2)

  # Sanity checks on the surface and the variance.
  perm <- obs
  perm$cases <- sample(perm$cases)
  ok_perm <- krige(cases ~ 1, perm, newdata = cell_pts, model = vfit, debug.level = 0)
  in_hull_pred <- cells$pred[cells$in_hull]
  sanity <- list(
    pred_range = range(cells$pred), negative_preds = sum(cells$pred < 0),
    mean_pred_hull = mean(in_hull_pred), sd_pred_hull = sd(in_hull_pred),
    sd_obs = sd(sites$cases), var_range = range(cells$var),
    var_in_hull = mean(cells$var[cells$in_hull]), var_out_hull = mean(cells$var[!cells$in_hull]),
    var_distance_rho = cor(cells$var, cells$nearest_site_km, method = "spearman"),
    permutation_max_var_diff = max(abs(ok_perm$var1.var - cells$var)),
    permutation_pred_cor = cor(ok_perm$var1.pred, cells$pred),
    share_outside_hull = mean(!cells$in_hull), n_cells = nrow(cells),
    share_outside_area = sum(cells$land_km2[!cells$in_hull]) / sum(cells$land_km2),
    hull_land_km2_per_site = as.numeric(st_area(st_intersection(hull, land))) / 1e6 / nrow(obs),
    n_in_hull = sum(cells$in_hull), thresholds = thresholds,
    high_km2 = sum(high$land_km2), top10_km2 = sum(high$land_km2[high$risk_class == "Top 10%"]))

  # Exceedance from 200 conditional simulations against the kriged map (Exercise 59).
  # Each square in the sampled area is simulated at its prediction point and at a
  # 4 x 4 lattice of land points. Their combined mean approximates the average over
  # the land portion, while the main kriged map retains point support.
  # Point values are kept for comparison.
  hull_ids <- cells$cell_id[cells$in_hull]
  square_center <- floor(st_coordinates(cell_pts[cells$in_hull, ]) / 1000) * 1000 + 500
  offsets <- expand.grid(dx = c(-375, -125, 125, 375), dy = c(-375, -125, 125, 375))
  lattice <- st_as_sf(data.frame(cell_id = rep(hull_ids, each = nrow(offsets)),
                                 x = rep(square_center[, 1], each = nrow(offsets)) + offsets$dx,
                                 y = rep(square_center[, 2], each = nrow(offsets)) + offsets$dy),
                      coords = c("x", "y"), crs = svy21)
  lattice <- lattice[lengths(st_intersects(lattice, land)) > 0, ]
  sim_pts <- rbind(st_sf(cell_id = hull_ids, point = TRUE,
                         geometry = st_geometry(cell_pts)[cells$in_hull]),
                   st_sf(cell_id = lattice$cell_id, point = FALSE, geometry = st_geometry(lattice)))
  sims <- krige(cases ~ 1, obs, newdata = sim_pts, model = vfit, nsim = 200, nmax = 30,
                debug.level = 0)
  sim_values <- as.matrix(st_drop_geometry(sims))
  point_sims <- sim_values[sim_pts$point, ]
  block_sims <- rowsum(sim_values, sim_pts$cell_id) / as.vector(table(sim_pts$cell_id))
  block_sims <- block_sims[as.character(hull_ids), ]
  cells$p_exceed3 <- NA_real_
  cells$p_exceed3[cells$in_hull] <- rowMeans(block_sims > 3)
  exceedance <- tibble(threshold = c(2, 3, 4),
                       from_kriged_map = map_dbl(threshold, \(t) mean(in_hull_pred > t)),
                       block_simulations = map_dbl(threshold, \(t) mean(block_sims > t)),
                       point_simulations = map_dbl(threshold, \(t) mean(point_sims > t)))
  save_tab(exceedance, "exceedance_map_vs_simulation.csv")

  # LOOCV predictions for the figure.
  cv_sel <- krige.cv(cases ~ 1, obs, model = vfit, nfold = nrow(obs), verbose = FALSE)
  cv_ctl <- krige.cv(cases ~ 1, obs, model = control_model, nfold = nrow(obs), verbose = FALSE)
  cv_points <- bind_rows(
    tibble(model = paste(selected_name, "model (selected)"), observed = cv_sel$observed,
           predicted = cv_sel$var1.pred),
    tibble(model = "Pure nugget control", observed = cv_ctl$observed,
           predicted = cv_ctl$var1.pred)) |>
    mutate(model = factor(model, unique(model)))

  # Maps. Fill and outline share one scale so squares show no seams.
  outside <- filter(cells, !in_hull)
  hull_line <- geom_sf(data = st_sf(geometry = hull), fill = NA, color = "black",
                       linetype = "dashed", linewidth = 0.35)
  fig_pred <- ggplot() + map_base() +
    geom_sf(data = cells, aes(fill = pred, color = after_scale(fill)), linewidth = 0.05) +
    geom_sf(data = outside, fill = "white", color = NA, alpha = 0.45) +
    hull_line +
    geom_sf(data = obs, aes(size = cases), shape = 21, fill = NA, color = "black", stroke = 0.6) +
    scale_fill_distiller(palette = "YlOrRd", direction = 1, name = "Predicted\ncases") +
    scale_size_area(max_size = 4.5, breaks = c(0, 2, 4, 6), name = "Observed\ncases")
  fig_var <- ggplot() + map_base() +
    geom_sf(data = cells, aes(fill = var, color = after_scale(fill)), linewidth = 0.05) +
    hull_line +
    geom_sf(data = obs, shape = 21, size = 1.6, fill = "white", color = "black") +
    scale_fill_viridis_c(option = "viridis", name = "Kriging\nvariance\n(cases²)")
  save_fig(fig_pred / fig_var, "fig06_kriging_maps.png", 8, 9)

  fig_zones <- ggplot() + map_base() +
    geom_sf(data = cells, aes(fill = risk_class, color = after_scale(fill)), linewidth = 0.05) +
    hull_line +
    geom_sf(data = zones, fill = NA, color = "black", linewidth = 0.5) +
    geom_sf(data = obs, shape = 21, size = 1.5, fill = "white", color = "black") +
    scale_fill_manual(values = c("Top 10%" = "#a50f15", "Top 25%" = "#fb6a4a", Other = "grey88"),
                      breaks = c("Top 10%", "Top 25%", "Other"),
                      labels = c("Top 10%" = "Top 10% (90th percentile\nand above)",
                                 "Top 25%" = "75th to 90th percentile",
                                 Other = "Lower 75%, or outside\nsampled area"),
                      name = "Predicted value\nwithin sampled area")
  save_fig(fig_zones, "fig07_high_risk_zones.png", 8, 4.8)

  fig_checks <- (ggplot(st_drop_geometry(cells), aes(nearest_site_km, var, color = in_hull)) +
                   geom_point(size = 0.8, alpha = 0.6) +
                   scale_color_manual(values = c(`TRUE` = "#08306b", `FALSE` = "#bdbdbd"),
                                      labels = c(`TRUE` = "Inside sampled area",
                                                 `FALSE` = "Outside (extrapolated)"),
                                      name = NULL) +
                   labs(x = "Distance to nearest site (km)", y = "Kriging variance (cases²)",
                        subtitle = "A. Variance rises with distance from data") +
                   theme_minimal(base_size = 10) + theme(legend.position = "bottom")) |
    (ggplot(cv_points, aes(observed, predicted)) +
       geom_abline(linetype = "dashed", color = "grey50") +
       geom_point(color = "#08306b", alpha = 0.75) +
       facet_wrap(~ model) +
       coord_equal(xlim = c(0, 6), ylim = c(0, 6)) +
       labs(x = "Observed cases", y = "Leave-one-out prediction",
            subtitle = "B. Cross-validation") +
       theme_minimal(base_size = 10))
  save_fig(fig_checks + plot_layout(widths = c(1, 1.6)), "fig08_sanity_checks.png", 11, 4.4)

  fig_sim <- ggplot() + map_base() +
    geom_sf(data = filter(cells, in_hull), aes(fill = p_exceed3, color = after_scale(fill)),
            linewidth = 0.05) +
    geom_sf(data = zones, fill = NA, color = "black", linewidth = 0.5) +
    geom_sf(data = obs, shape = 21, size = 1.5, fill = "white", color = "black") +
    scale_fill_viridis_c(option = "magma", direction = -1, limits = c(0, 1),
                         name = "P(1 km² mean > 3)") +
    coord_sf(xlim = st_bbox(hull)[c(1, 3)] + c(-1500, 1500),
             ylim = st_bbox(hull)[c(2, 4)] + c(-1500, 1500))
  save_fig(fig_sim, "fig09_exceedance_probability.png", 7.5, 5)

  list(
    cells = cells,
    cell_pts = cell_pts,
    hull = hull,
    high = high,
    zones = zones,
    zone_table = zone_table,
    zone_block_total = zone_block_total,
    sanity = sanity,
    exceedance = exceedance,
    cv_sel = cv_sel,
    lattice = lattice,
    hull_line = hull_line
  )
}
