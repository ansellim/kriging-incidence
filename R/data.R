# prepare_data: explicit inputs, with figures and tables written through config.
prepare_data <- function(config) {
  svy21 <- config$crs
  save_fig <- config$save_fig
  save_tab <- config$save_tab

  # ---- 1. Data preparation and exploration ------------------------------------

  sites <- read.csv(config$sites_file)
  required <- c("id", "latitude", "longitude", "x", "y", "green_space", "cases")
  if (!identical(names(sites), required) || !all(vapply(sites, is.numeric, logical(1))))
    stop("Site data must contain numeric columns: ", paste(required, collapse = ", "))
  if (nrow(sites) != 30 || anyNA(sites) || any(!is.finite(as.matrix(sites))) ||
      anyDuplicated(sites$id) || any(sites$cases < 0 | sites$cases != floor(sites$cases)) ||
      any(sites$green_space < 0 | sites$green_space > 100))
    stop("Expected 30 unique sites, finite coordinates, integer counts, and green space in [0, 100].")

  # The spatial object is built from the projected x/y columns (meters).
  obs <- st_as_sf(sites, coords = c("x", "y"), crs = svy21, remove = FALSE)

  # Check: latitude/longitude (WGS84, EPSG:4326) projected to SVY21 should give x/y.
  lonlat <- st_as_sf(sites, coords = c("longitude", "latitude"), crs = 4326)
  coord_check_m <- max(as.numeric(
    st_distance(obs, st_transform(lonlat, svy21), by_element = TRUE)))

  # Singapore land: URA 2019 subzones without sea, reprojected and dissolved.
  subzones <- st_read(config$land_file, quiet = TRUE) |>
    st_make_valid() |>
    st_transform(svy21) |>
    transmute(planning_area = title_case(PLN_AREA_N))
  land <- st_make_valid(st_union(subzones))
  land_km2 <- as.numeric(st_area(land)) / 1e6

  # Planning area of each site, used only to describe locations.
  obs <- st_join(obs, subzones, join = st_within)
  stopifnot(nrow(obs) == 30, !anyNA(obs$planning_area))

  describe <- function(x) {
    tibble(n = length(x), mean = mean(x), sd = sd(x), median = median(x),
           q1 = quantile(x, 0.25, names = FALSE), q3 = quantile(x, 0.75, names = FALSE),
           min = min(x), max = max(x))
  }
  summary_stats <- bind_rows(green_space = describe(sites$green_space),
                             cases = describe(sites$cases), .id = "variable")
  save_tab(summary_stats, "summary_statistics.csv")

  # Site spacing. Straight-line distance is valid in meters because SVY21 is projected.
  dist_m <- as.matrix(dist(sites[, c("x", "y")]))
  nearest_m <- apply(dist_m + diag(Inf, nrow(sites)), 1, min)
  closest <- which(dist_m == min(nearest_m), arr.ind = TRUE)[1, ]
  spacing <- list(pairs = choose(nrow(sites), 2),
                  nn_median_km = median(nearest_m) / 1000,
                  nn_max_km = max(nearest_m) / 1000,
                  max_sep_km = max(dist_m) / 1000,
                  extent_km = c(x = diff(range(sites$x)), y = diff(range(sites$y))) / 1000,
                  closest_ids = sort(sites$id[closest]), closest_m = min(nearest_m))

  # Shared map layers: land outline, scale bar, north arrow.
  map_base <- function() {
    list(geom_sf(data = land, fill = "grey96", color = "grey65", linewidth = 0.2),
         annotation_scale(location = "bl", width_hint = 0.2, text_cex = 0.7),
         annotation_north_arrow(location = "tr", style = north_arrow_orienteering(text_size = 7),
                                height = unit(0.9, "cm"), width = unit(0.6, "cm")),
         theme_void(base_size = 10),
         theme(legend.position = "right", plot.margin = margin(4, 4, 4, 4)))
  }
  site_xy <- bind_cols(st_drop_geometry(obs), as.data.frame(st_coordinates(obs)))

  fig_locations <- ggplot() + map_base() +
    geom_sf(data = obs, size = 1.8, color = "black") +
    geom_text_repel(data = site_xy, aes(X, Y, label = id), size = 2.6, seed = 1,
                    box.padding = 0.25, segment.color = "grey50")
  save_fig(fig_locations, "fig01_locations.png", 7.5, 4.6)

  fig_site_cases <- ggplot() + map_base() +
    geom_sf(data = filter(obs, cases > 0), aes(size = cases), shape = 21, fill = "#d7301f",
            color = "black", alpha = 0.8) +
    geom_sf(data = filter(obs, cases == 0), aes(shape = "No cases"), size = 2, stroke = 0.8) +
    scale_size_area(max_size = 6, breaks = c(1, 2, 4, 6), name = "Cases") +
    scale_shape_manual(values = c("No cases" = 4), name = NULL) +
    labs(subtitle = "A. Case count at each site")
  fig_site_green <- ggplot() + map_base() +
    geom_sf(data = obs, aes(fill = green_space), shape = 21, size = 3, color = "black") +
    scale_fill_distiller(palette = "Greens", direction = 1, name = "Green\nspace (%)") +
    labs(subtitle = "B. Green space around each site")
  fig_hist_green <- ggplot(sites, aes(green_space)) +
    geom_histogram(binwidth = 5, boundary = 25, fill = "#41ab5d", color = "white") +
    geom_vline(xintercept = mean(sites$green_space), linetype = "dashed") +
    labs(x = "Green space (%)", y = "Sites", subtitle = "C. Green space (dashed line: mean)") +
    theme_minimal(base_size = 10)
  fig_bar_cases <- ggplot(sites, aes(cases)) +
    geom_bar(fill = "#d7301f", width = 0.7) +
    scale_x_continuous(breaks = 0:6) +
    scale_y_continuous(breaks = seq(0, 10, 2)) +
    labs(x = "Cases at site", y = "Sites", subtitle = "D. Case counts") +
    theme_minimal(base_size = 10)
  save_fig((fig_site_cases | fig_site_green) / (fig_hist_green | fig_bar_cases) +
             plot_layout(heights = c(1.25, 1)), "fig02_site_attributes.png", 10, 7)

  list(
    sites = sites,
    obs = obs,
    coord_check_m = coord_check_m,
    subzones = subzones,
    land = land,
    land_km2 = land_km2,
    summary_stats = summary_stats,
    spacing = spacing,
    dist_m = dist_m,
    map_base = map_base
  )
}
