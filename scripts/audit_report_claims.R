# Recompute evidence for the numerical interpretation audit without changing
# analytical outputs. Run from the repository root after run_analysis.R.
required <- c("sf", "gstat", "dplyr", "knitr", "jsonlite")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({ library(sf); library(gstat) })
out <- file.path("outputs", "audit")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
report <- readLines("report.qmd", warn = FALSE, encoding = "UTF-8")
first <- match("```{r}", report)
if (is.na(first)) stop("Report setup chunk not found.")
last <- which(seq_along(report) > first & report == "```")[1]
if (is.na(last)) stop("Report setup chunk is not terminated.")
env <- new.env(parent = globalenv())
eval(parse(text = report[(first + 1):(last - 1)]), envir = env)

# Evaluate every inline R expression, retaining source lines and rendered values.
# This checks executable report claims, not the scientific validity of their
# prose. The separate reviewed_claims.json records that interpretive review.
inline <- lapply(seq_along(report), function(i) {
  hits <- regmatches(report[i], gregexpr("`r [^`]+`", report[i]))[[1]]
  if (!length(hits)) return(NULL)
  rows <- lapply(hits, function(hit) {
    expression <- substr(hit, 4, nchar(hit) - 1)
    value <- tryCatch(eval(parse(text = expression), envir = env),
                      error = function(e) stop("Inline expression at line ", i,
                                               " failed: ", conditionMessage(e)))
    if (is.numeric(value) && any(!is.finite(value)))
      stop("Nonfinite inline result at line ", i)
    data.frame(line = i, expression = expression,
               rendered_value = paste(as.character(value), collapse = "; "))
  })
  do.call(rbind, rows)
})
inline <- do.call(rbind, inline)
if (is.null(inline) || nrow(inline) == 0) stop("No inline report expressions found.")
write.csv(inline, file.path(out, "report_inline_values.csv"), row.names = FALSE)

sites <- read.csv("data/raw/assignment2.csv")
r <- readRDS("outputs/results.rds")
grid <- read.csv("outputs/tables/kriging_grid.csv")
xy <- as.matrix(sites[, c("x", "y")])
D <- as.matrix(dist(xy))
nn <- D
diag(nn) <- Inf
range_m <- r$vfit$range[2]
nugget <- r$vfit$psill[1]
partial <- r$vfit$psill[2]
sill <- nugget + partial
stopifnot(as.character(r$vfit$model[2]) == "Sph")
C <- partial * ifelse(D < range_m,
                       1 - 1.5 * D / range_m + 0.5 * (D / range_m)^3, 0)
diag(C) <- sill
inverse_one <- solve(C, rep(1, nrow(C)))
gls_mean <- sum(inverse_one * sites$cases) / sum(inverse_one)
mean_variance <- 1 / sum(inverse_one)
remote <- grid$nearest_site_km > range_m / 1000
stopifnot(any(remote),
          max(abs(grid$pred[remote] - gls_mean)) < 1e-8,
          max(abs(grid$var[remote] - sill - mean_variance)) < 1e-8)

obs <- st_as_sf(sites, coords = c("x", "y"), crs = 3414)
hull <- st_convex_hull(st_union(obs))
land <- st_read("data/cache/singapore_land.geojson", quiet = TRUE) |>
  st_make_valid() |> st_transform(3414) |> st_union() |> st_make_valid()
hull_land <- as.numeric(st_area(st_intersection(hull, land))) / 1e6
stopifnot(abs(hull_land / nrow(sites) - r$sanity$hull_land_km2_per_site) < 1e-8)
points <- st_as_sf(grid, coords = c("px", "py"), crs = 3414, remove = FALSE)
point_in_hull <- lengths(st_intersects(points, hull)) > 0

# Reconstruct the actual grid-cell land portions to verify the facility-location
# claims and identify any boundary prediction points outside the strict hull.
origin <- floor(st_bbox(land)[c("xmin", "ymin")] / 1000) * 1000
squares <- st_make_grid(land, cellsize = 1000, offset = origin)
cells <- st_sf(cell_id = seq_along(squares), geometry = squares)
st_agr(cells) <- "constant"
cells <- st_intersection(cells, land) |> st_collection_extract("POLYGON")
high <- grid$risk_class != "Other"
high_shape <- st_union(cells[cells$cell_id %in% grid$cell_id[high], ])
core_shape <- st_union(cells[cells$cell_id %in% grid$cell_id[grid$risk_class == "Top 10%"], ])
shortlist <- read.csv("outputs/tables/center_shortlist.csv")
centers <- shortlist[match(r$centers$name, shortlist$name), ]
stopifnot(!anyNA(centers))
center_points <- st_as_sf(centers, coords = c("X", "Y"), crs = 3414)
centers_in_zone <- lengths(st_intersects(center_points, high_shape)) > 0
centers_in_core <- lengths(st_intersects(center_points, core_shape)) > 0

chosen <- paste(sort(r$centers$name), collapse = " + ")
same_pair <- r$sensitivity$centers == chosen
ordinary <- !grepl("Universal", r$sensitivity$scenario)
lowest_rmse <- r$candidates[which.min(r$candidates$rmse), ]
lowest_mae <- r$candidates[which.min(r$candidates$mae), ]
stopifnot(lowest_rmse$model == r$selected$model,
          lowest_mae$model == "Exp", sum(same_pair) == 5,
          all(same_pair[ordinary]), !same_pair[!ordinary])
v_first <- lapply(c(1000, 1500, 2000, 2500, 3000), function(width) {
  v <- variogram(cases ~ 1, obs, width = width, cutoff = 10000)
  data.frame(width_m = width, mean_lag_m = v$dist[1],
             pairs = v$np[1], semivariance = v$gamma[1])
})
v_first <- do.call(rbind, v_first)
v_tail <- variogram(cases ~ 1, obs, width = 1500, cutoff = 11000)
tail_bin <- tail(v_tail[, c("np", "dist", "gamma")], 1)
core <- grid$risk_class == "Top 10%"
stopifnot(max(grid$p_exceed3, na.rm = TRUE) == max(grid$p_exceed3[core]))

evidence <- list(
  site_count = nrow(sites), total_cases = sum(sites$cases), zero_count_sites = sum(sites$cases == 0),
  count_range = range(sites$cases), sample_mean = mean(sites$cases), sample_variance = var(sites$cases),
  coordinate_reprojection_max_m = r$coord_check_m,
  sampled_hull_land_km2 = hull_land, land_proxy_km2 = as.numeric(st_area(land)) / 1e6,
  illustrative_catchment_km2 = hull_land / nrow(sites),
  conditional_zone_density_1km2 = r$zone_table$est_per_km2,
  conditional_zone_density_hull_per_site = r$zone_table$est_per_km2 / (hull_land / nrow(sites)),
  catchment_ratio_is_a_scenario_not_an_uncertainty_estimate = TRUE,
  selected_rmse = r$selected$rmse, control_rmse = r$control_cv$rmse,
  relative_rmse_difference = 1 - r$selected$rmse / r$control_cv$rmse,
  selected_mae = r$selected$mae, lowest_mae = lowest_mae$mae,
  lowest_mae_model = lowest_mae$model,
  first_empirical_bins = v_first, final_1_5km_bin_at_11km_cutoff = tail_bin,
  grid_cells = nrow(grid), hull_overlapping_cells = sum(grid$in_hull),
  prediction_points_in_hull = sum(point_in_hull),
  hull_overlapping_prediction_points_outside_hull = sum(grid$in_hull & !point_in_hull),
  high_zone_prediction_points_outside_hull = sum(high & !point_in_hull),
  beyond_range_prediction_points = sum(remote),
  beyond_range_cell_land_km2 = sum(grid$land_km2[remote]),
  beyond_range_gls_mean = gls_mean, beyond_range_mean_estimation_variance = mean_variance,
  beyond_range_kriging_sd = sqrt(sill + mean_variance),
  loo_max_nearest_training_km = max(apply(nn, 1, min)) / 1000,
  loo_sites_beyond_selected_range = sum(apply(nn, 1, min) > range_m),
  poisson_count_ratio = r$greenspace_model$ratio,
  poisson_profile_interval = unname(r$greenspace_model$ratio_ci),
  poisson_deviance_reduction = r$greenspace_model$deviance_explained,
  raw_moran = r$moran$cases_I, pearson_residual_moran = r$moran$residual_I,
  attributable_fraction_of_spatial_clustering_estimated = FALSE,
  exceedance_frequencies = r$exceedance,
  single_cell_simulation_frequency_nominal_se_upper_bound = sqrt(0.25 / 200),
  maximum_exceedance_frequency = max(grid$p_exceed3, na.rm = TRUE),
  maximum_frequency_in_top10 = max(grid$p_exceed3[core]),
  top10_fraction_frequency_over_half = mean(grid$p_exceed3[core] > 0.5),
  center_names = r$centers$name, centers_in_zone = centers_in_zone, centers_in_core = centers_in_core,
  centers_mrt_straight_line_m = r$centers$mrt_m,
  sensitivity_scenarios = nrow(r$sensitivity), same_center_pair_scenarios = sum(same_pair),
  sensitivity_center_pairs = r$sensitivity[, c("scenario", "centers")],
  distinction = "Numerical consistency and conditional model calculations do not establish causal, predictive, or operational validity.")
jsonlite::write_json(evidence, file.path(out, "claim_numerics.json"), pretty = TRUE,
                     auto_unbox = TRUE, digits = 15, na = "null")
reviews <- jsonlite::fromJSON(file.path("scripts", "reviewed_claims.json"), simplifyVector = FALSE)
reviews$inline_expressions_evaluated <- nrow(inline)
reviews$numerical_evidence_file <- "claim_numerics.json"
reviews$reviewed_report_md5 <- unname(tools::md5sum("report.qmd"))
reviews$generated_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
jsonlite::write_json(reviews, file.path(out, "claim_review.json"), pretty = TRUE,
                     auto_unbox = TRUE, digits = 15, na = "null")
message("Claim evidence verified; ", nrow(inline), " inline expressions evaluated; ",
        length(reviews$reviewed_changes), " reviewed edits recorded.")
