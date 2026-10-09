# analyze_sensitivity: explicit inputs, with figures and tables written through config.
analyze_sensitivity <- function(data, variograms, surface, facilities, config) {
  svy21 <- config$crs
  save_fig <- config$save_fig
  save_tab <- config$save_tab
  obs <- data$obs
  candidates <- variograms$candidates
  selected <- variograms$selected
  selected_name <- variograms$selected_name
  vfit <- variograms$vfit
  cell_pts <- surface$cell_pts
  cells <- surface$cells
  high <- surface$high
  clubs <- facilities$clubs

  # ---- 6. Sensitivity analyses ------------------------------------------------
  # Re-krige under alternative defensible choices and compare with the selected
  # model: cross-validation, overlap of the top-25% squares (Jaccard index), and
  # whether the same two centers are chosen.

  # Refit with the sensible start; if no structured model converges, fit a pure
  # nugget instead and say so (the structure was absorbed elsewhere).

  cell_xy <- bind_cols(cell_pts, as.data.frame(st_coordinates(cell_pts))) |> rename(x = X, y = Y)
  no25 <- filter(obs, id != 25)
  v_uk <- variogram(cases ~ x + y, obs, width = selected$width, cutoff = selected$cutoff)
  scenarios <- list(
    list(name = paste(selected_name, "(selected)"), model = vfit, formula = cases ~ 1, data = obs),
    list(name = "Exponential, same bins", model = pick_fit(candidates, "Exp", 3000, 10000, 7),
         formula = cases ~ 1, data = obs),
    list(name = "Gaussian, same bins, pair weights", model = pick_fit(candidates, "Gau", 3000, 10000, 1),
         formula = cases ~ 1, data = obs),
    list(name = "Spherical, 2.5 km bins, pair weights", model = pick_fit(candidates, "Sph", 2500, 10000, 1),
         formula = cases ~ 1, data = obs),
    list(name = "Spherical, site 25 removed",
         model = fit_or_nugget(variogram(cases ~ 1, no25, width = 3000, cutoff = 10000), "Sph"),
         formula = cases ~ 1, data = no25),
    list(name = "Universal kriging, linear trend in x and y", model = fit_or_nugget(v_uk, "Sph"),
         formula = cases ~ x + y, data = obs))
  jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))
  sensitivity <- map_dfr(scenarios, \(s) {
    k <- krige(s$formula, s$data, newdata = cell_xy, model = s$model, debug.level = 0)
    in_hull_k <- k$var1.pred[cells$in_hull]
    top <- cells$in_hull & k$var1.pred >= quantile(in_hull_k, 0.75)
    top_squares <- mutate(cells, pred = k$var1.pred)[top, ]
    pick <- choose_centers(top_squares, st_union(top_squares), clubs, cell_pts, cells)
    sp <- s$model[s$model$model != "Nug", ]
    bind_cols(tibble(scenario = s$name, nugget = sum(s$model$psill[s$model$model == "Nug"]),
                     psill = if (nrow(sp)) sp$psill else 0,
                     range_km = if (nrow(sp)) practical_range(as.character(sp$model), sp$range) / 1000
                                else NA_real_),
              loocv(s$model, s$formula, s$data),
              tibble(pred_min = min(k$var1.pred), pred_max = max(k$var1.pred),
                     top25_jaccard = jaccard(cells$cell_id[top], high$cell_id),
                     centers = paste(sort(pick$centers$name), collapse = " + ")))
  })
  save_tab(sensitivity, "sensitivity.csv")

  list(sensitivity = sensitivity)
}
