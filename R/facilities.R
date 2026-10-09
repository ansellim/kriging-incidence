# place_centers: explicit inputs, with figures and tables written through config.
place_centers <- function(data, surface, config) {
  svy21 <- config$crs
  save_fig <- config$save_fig
  save_tab <- config$save_tab
  land <- data$land
  subzones <- data$subzones
  obs <- data$obs
  map_base <- data$map_base
  high <- surface$high
  zones <- surface$zones
  cell_pts <- surface$cell_pts
  cells <- surface$cells
  hull <- surface$hull

  # ---- 5. Intervention centers ------------------------------------------------
  # Candidates are community clubs: public, staffed buildings within residential
  # areas. A club qualifies if it is within 500 m (straight line) of an MRT exit
  # and within 1 km of the top-25% zone. Among all pairs of qualifying clubs, the
  # chosen pair minimizes the average distance from high-risk squares to the
  # nearer club, weighting each square by its estimated cases (a two-facility
  # p-median problem).
  # The 2024 dataset flags some clubs "(U/C)", "(Pending U/C)" or "(Pending
  # relocation)". The current official directory (checked 8 October 2026) lists
  # the flagged clubs near the zone as operating, so the flags are treated as out
  # of date and removed. The pre-2022 Thomson CC location and venues that are not
  # community clubs (sea sports clubs, imPAct@Hong Lim Green) are dropped.
  clubs <- read.csv(config$clubs_file) |>
    filter(!grepl("Till", name)) |>
    mutate(name = trimws(sub("/.*$", "", sub("\\s*\\(.*\\)$", "", name)))) |>
    filter(grepl(" CC$", name)) |>
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) |>
    st_transform(svy21)
  clubs <- clubs[lengths(st_intersects(clubs, land)) > 0, ]
  exits <- st_read(config$exits_file, quiet = TRUE) |>
    st_transform(svy21)
  nearest_exit <- st_nearest_feature(clubs, exits)
  clubs <- clubs |>
    mutate(mrt_station = title_case(exits$STATION_NA[nearest_exit]),
           mrt_exit = exits$EXIT_CODE[nearest_exit],
           mrt_m = as.numeric(st_distance(clubs, exits[nearest_exit, ], by_element = TRUE))) |>
    st_join(subzones, join = st_within)


  placement <- choose_centers(high, st_union(zones), clubs, cell_pts, cells)
  centers <- placement$centers |>
    mutate(est_cases_served = as.numeric(tapply(placement$weight, factor(placement$nearer, 1:2), sum)),
           squares_served = tabulate(placement$nearer, 2),
           mean_km_served = as.numeric(tapply(placement$weight * placement$dist_km,
                                              factor(placement$nearer, 1:2), sum)) /
             est_cases_served)
  verification <- read.csv(config$verification_file)
  centers <- left_join(centers, select(verification, name, verified_source = source,
                                       access_notes), by = "name")
  facility_summary <- list(
    clubs = nrow(clubs), shortlist = nrow(placement$shortlist), pairs = choose(nrow(placement$shortlist), 2),
    score_km = placement$score_km,
    within_2km = sum(placement$weight[placement$dist_km <= 2]) / sum(placement$weight),
    max_km = max(placement$dist_km), runner_up = placement$ranking[2, ])
  save_tab(select(centers, name, address, planning_area, mrt_station, mrt_exit, mrt_m,
                  squares_served, est_cases_served, mean_km_served), "intervention_centers.csv")
  save_tab(head(placement$ranking, 10), "center_pair_ranking.csv")
  save_tab(bind_cols(st_drop_geometry(placement$shortlist)[, c("name", "mrt_m")],
                     as.data.frame(st_coordinates(placement$shortlist))), "center_shortlist.csv")

  near_exits <- st_filter(exits, st_as_sfc(st_bbox(st_buffer(hull, 2500))))
  center_xy <- bind_cols(st_drop_geometry(centers), as.data.frame(st_coordinates(centers))) |>
    mutate(label = sprintf("%s\n(%s, %s, %.0f m)", name, mrt_station, mrt_exit, mrt_m))
  fig_centers <- ggplot() + map_base() +
    geom_sf(data = filter(cells, in_hull), aes(fill = risk_class, color = after_scale(fill)),
            linewidth = 0.05) +
    geom_sf(data = zones, fill = NA, color = "black", linewidth = 0.5) +
    geom_sf(data = st_buffer(centers, 2000), fill = NA, color = "#08306b", linetype = "dashed",
            linewidth = 0.45) +
    geom_sf(data = near_exits, aes(shape = "MRT exit"), size = 0.6, color = "grey35") +
    geom_sf(data = placement$shortlist, aes(shape = "Shortlisted community club"), size = 1.8,
            color = "#08306b") +
    geom_sf(data = obs, aes(shape = "Sampled site"), size = 1.4, fill = "white", color = "black") +
    geom_sf(data = centers, aes(shape = "Proposed center"), size = 4, fill = "#2171b5",
            color = "black") +
    scale_shape_manual(values = c("Proposed center" = 23, "Shortlisted community club" = 2,
                                  "Sampled site" = 21, "MRT exit" = 16), name = NULL,
                       breaks = c("Proposed center", "Shortlisted community club", "Sampled site",
                                  "MRT exit"),
                       guide = guide_legend(override.aes = list(
                         size = c(3.5, 1.8, 1.6, 1.2), fill = c("#2171b5", NA, "white", NA),
                         color = c("black", "#08306b", "black", "grey35")))) +
    geom_label_repel(data = center_xy, aes(X, Y, label = label), size = 2.7, seed = 2,
                     nudge_x = -7500, nudge_y = ifelse(center_xy$Y == min(center_xy$Y), -800, 1200),
                     min.segment.length = 0, label.size = 0.2, fill = "white") +
    scale_fill_manual(values = c("Top 10%" = "#a50f15", "Top 25%" = "#fb6a4a", Other = "grey90"),
                      breaks = c("Top 10%", "Top 25%", "Other"),
                      labels = c("Top 10%" = "Top 10% (90th percentile\nand above)",
                                 "Top 25%" = "75th to 90th percentile", Other = "Lower 75%"),
                      name = "Predicted value\nwithin sampled area") +
    coord_sf(xlim = st_bbox(hull)[c(1, 3)] + c(-2500, 2500),
             ylim = st_bbox(hull)[c(2, 4)] + c(-2500, 2500))
  save_fig(fig_centers, "fig12_intervention_centers.png", 8, 5.6)

  list(clubs = clubs, centers = centers, placement = placement, facility_summary = facility_summary)
}
