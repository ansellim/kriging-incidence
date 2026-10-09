# Reusable numerical helpers. Distances are in EPSG:3414 meters.

profile_zones <- function(squares, label, obs, subzones, cells, cell_pts) {
  zones <- st_sf(geometry = st_cast(st_union(squares), "POLYGON")) |>
    mutate(area_km2 = as.numeric(st_area(geometry)) / 1e6) |>
    filter(area_km2 >= 0.5) |>
    mutate(zone = paste0(label, " ", row_number()))
  in_zone <- st_join(st_sf(st_drop_geometry(squares), geometry = st_geometry(cell_pts)[
    match(squares$cell_id, cells$cell_id)]), zones, join = st_within)
  sites_in <- st_join(obs, zones, join = st_within) |> st_drop_geometry() |> filter(!is.na(zone))
  areas <- suppressWarnings(st_intersection(zones, subzones)) |>
    mutate(a = as.numeric(st_area(geometry)) / 1e6) |>
    st_drop_geometry() |>
    group_by(zone, planning_area) |>
    summarise(a = sum(a), .groups = "drop") |>
    filter(a >= 1) |>
    arrange(zone, desc(a)) |>
    group_by(zone) |>
    summarise(planning_areas = paste(planning_area, collapse = ", "))
  st_drop_geometry(in_zone) |>
    filter(!is.na(zone)) |>
    group_by(zone, area_km2) |>
    summarise(squares = n(), mean_pred = mean(pred), max_pred = max(pred),
              est_cases = sum(pred * land_km2), .groups = "drop") |>
    mutate(est_per_km2 = est_cases / area_km2) |>
    left_join(sites_in |> group_by(zone) |>
                summarise(sites = n(), site_ids = paste(sort(id), collapse = ", "),
                          recorded_cases = sum(cases)), by = "zone") |>
    mutate(sites = coalesce(sites, 0L), recorded_cases = coalesce(recorded_cases, 0L),
           recorded_per_km2 = recorded_cases / area_km2) |>
    left_join(areas, by = "zone") |>
    list(table = _, zones = zones)
}
