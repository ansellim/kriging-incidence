# Reusable numerical helpers. Distances are in EPSG:3414 meters.

choose_centers <- function(high_squares, zone_union, clubs, cell_pts, cells) {
  shortlist <- filter(clubs, mrt_m <= 500,
                      as.numeric(st_distance(clubs, zone_union))[seq_len(nrow(clubs))] <= 1000)
  demand <- cell_pts[match(high_squares$cell_id, cells$cell_id), ]
  weight <- high_squares$pred * high_squares$land_km2   # estimated cases per square
  d_km <- units::drop_units(st_distance(demand, shortlist)) / 1000
  pairs <- combn(nrow(shortlist), 2)
  score <- apply(pairs, 2, \(p) sum(weight * pmin(d_km[, p[1]], d_km[, p[2]])) / sum(weight))
  best <- pairs[, which.min(score)]
  nearer <- apply(d_km[, best], 1, which.min)
  list(shortlist = shortlist, centers = shortlist[best, ], score_km = min(score),
       ranking = tibble(club_a = shortlist$name[pairs[1, ]], club_b = shortlist$name[pairs[2, ]],
                        mean_km = score) |> arrange(mean_km),
       nearer = nearer, dist_km = apply(d_km[, best], 1, min), weight = weight)
}
