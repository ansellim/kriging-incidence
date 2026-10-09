# fit_variograms: explicit inputs, with figures and tables written through config.
fit_variograms <- function(data, config) {
  svy21 <- config$crs
  save_fig <- config$save_fig
  save_tab <- config$save_tab
  sites <- data$sites
  obs <- data$obs
  dist_m <- data$dist_m
  model_names <- c(Sph = "Spherical", Exp = "Exponential", Gau = "Gaussian")

  # ---- 2. Variogram construction and fitting ----------------------------------
  # Semivariance gamma(h) = half the mean squared difference in cases between pairs
  # of sites about h meters apart. Units: cases squared.

  sample_var <- var(sites$cases)
  max_sep_m <- max(dist_m)

  # 2a. Empirical variograms over a grid of bin widths and cutoffs (Exercise 52).
  widths <- c(1000, 1500, 2000, 2500, 3000)
  cutoffs <- c(8000, 9000, 10000, 11000)  # all at most half the maximum separation
  empirical <- expand_grid(width = widths, cutoff = cutoffs) |>
    mutate(v = map2(width, cutoff, \(w, ct) variogram(cases ~ 1, obs, width = w, cutoff = ct)),
           bins = map_int(v, nrow),
           min_pairs = map_int(v, \(x) min(x$np)),
           first_bin_pairs = map_int(v, \(x) x$np[1]))
  save_tab(select(empirical, -v), "variogram_bin_settings.csv")

  # 2b. Fit spherical, exponential and Gaussian models under three weighting
  # schemes and four sets of starting values (Exercise 53). The practical range is
  # reported for every model: Sph = range, Exp = 3 x range, Gau = sqrt(3) x range.




  starts <- tribble(
    ~start,               ~start_nugget, ~start_psill, ~start_range_m,
    "sensible",           1.5,           1.5,          6000,   # read off the plot
    "range too short",    1.5,           1.5,          600,    # poor starts
    "range too long",     1.5,           1.5,          60000,
    "nugget equals sill", 2.9,           0.1,          6000)
  methods <- c("1" = "Weights: pairs (method 1)",
               "2" = "Weights: Cressie (method 2)",
               "7" = "Weights: pairs/distance² (default, method 7)")

  fits <- empirical |>
    select(width, cutoff, bins, min_pairs, v) |>
    expand_grid(model = c("Sph", "Exp", "Gau"), method = c(1, 2, 7), starts) |>
    mutate(result = pmap(list(v, model, start_nugget, start_psill, start_range_m, method),
                         fit_one)) |>
    select(-v) |>
    unnest(result)

  # A setting is "stable" when the sensible start converges cleanly, at least two
  # starts converge cleanly and agree within 1%, the parameters are plausible
  # (practical range within the cutoff, since no bin informs longer distances;
  # sill at most twice the sample variance) and there are more bins than fitted
  # parameters.
  agree <- function(x) (max(x) - min(x)) <= 0.01 * max(abs(x), 1e-9)
  tuning <- fits |>
    group_by(width, cutoff, bins, min_pairs, model, method) |>
    summarise(clean_starts = sum(clean),
              sensible_clean = clean[start == "sensible"],
              starts_agree = sum(clean) >= 2 && agree(nugget[clean]) &&
                agree(psill[clean]) && agree(range_m[clean]),
              nugget = nugget[start == "sensible"],
              psill = psill[start == "sensible"],
              range_m = range_m[start == "sensible"],
              model_fit = model_fit[start == "sensible"],
              .groups = "drop") |>
    mutate(sill = nugget + psill,
           plausible = nugget >= 0 & psill > 0 & range_m <= cutoff & sill <= 2 * sample_var,
           stable = sensible_clean & starts_agree & plausible & bins >= 4,
           status = case_when(stable ~ "Stable fit",
                              sensible_clean ~ "Converged, but unstable or implausible",
                              TRUE ~ "No convergence or singular"))

  # 2c. Leave-one-out cross-validation (Exercise 56) of every stable fit and of a
  # pure-nugget control, which asserts no spatial correlation (Week 6 slide 56).

  control_model <- vgm(psill = sample_var, model = "Nug", range = 0)
  control_cv <- loocv(control_model, data = obs)
  candidates <- tuning |>
    filter(stable) |>
    mutate(cv = map(model_fit, \(model) loocv(model, data = obs))) |>
    unnest(cv)

  # 2d. Selection rule, fixed before looking at the maps: prefer settings in which
  # every bin has at least 30 pairs (Exercise 52 rule), then the lowest LOOCV RMSE
  # among fits whose z-score SD is between 0.8 and 1.2 and that beat the control.
  pool <- filter(candidates, min_pairs >= 30)
  if (nrow(pool) == 0) pool <- candidates
  selected <- pool |>
    filter(z_sd > 0.8, z_sd < 1.2, rmse < control_cv$rmse) |>
    slice_min(rmse, n = 1, with_ties = FALSE)
  stopifnot(nrow(selected) == 1)
  vfit <- selected$model_fit[[1]]
  v_sel <- variogram(cases ~ 1, obs, width = selected$width, cutoff = selected$cutoff)
  selected_name <- unname(model_names[selected$model])

  save_tab(select(fits, -model_fit), "variogram_fit_attempts.csv")
  save_tab(select(tuning, -model_fit), "variogram_tuning.csv")
  save_tab(select(candidates, -model_fit), "variogram_candidates_cv.csv")

  # Figure: tuning map of all 180 settings.
  fig_tuning <- tuning |>
    mutate(model = factor(model_names[model], model_names),
           method = factor(methods[as.character(method)], methods),
           label = ifelse(sensible_clean, sprintf("%.0f", range_m / 1000), "")) |>
    ggplot(aes(factor(width / 1000), factor(cutoff / 1000), fill = status)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = label), size = 2.6) +
    facet_grid(method ~ model) +
    scale_fill_manual(values = c("Stable fit" = "#1a9850",
                                 "Converged, but unstable or implausible" = "#fee08b",
                                 "No convergence or singular" = "#d9d9d9"), name = NULL) +
    labs(x = "Bin width (km)", y = "Cutoff (km)") +
    theme_minimal(base_size = 9) +
    theme(legend.position = "bottom", panel.grid = element_blank())
  save_fig(fig_tuning, "fig03_variogram_tuning.png", 9, 7.5)

  # Figure: selected fit, and two settings where gstat returned parameters
  # despite a convergence warning.

  same_bins <- filter(candidates, width == selected$width, cutoff == selected$cutoff,
                      method == selected$method)
  failed <- fits |>
    filter(model == selected$model, method == 7, start == "sensible", !clean) |>
    filter((width == 2000 & cutoff == 10000) | (width == 1500 & cutoff == 11000))
  panels <- c(sprintf("A. %g km bins, %g km cutoff: converged", selected$width / 1000,
                      selected$cutoff / 1000),
              "B. 2 km bins, 10 km cutoff: no convergence",
              "C. 1.5 km bins, 11 km cutoff: no convergence")
  curves <- bind_rows(
    pmap(list(same_bins$model_fit, same_bins$model),
         \(m, nm) curve_of(m, model_names[nm], panels[1])),
    curve_of(failed$model_fit[[which(failed$width == 2000)]], "Returned despite warning",
             panels[2]),
    curve_of(failed$model_fit[[which(failed$width == 1500)]], "Returned despite warning",
             panels[3]))
  points <- bind_rows(
    mutate(as_tibble(v_sel)[, c("np", "dist", "gamma")], panel = panels[1]),
    mutate(as_tibble(variogram(cases ~ 1, obs, width = 2000, cutoff = 10000))[
      , c("np", "dist", "gamma")], panel = panels[2]),
    mutate(as_tibble(variogram(cases ~ 1, obs, width = 1500, cutoff = 11000))[
      , c("np", "dist", "gamma")], panel = panels[3])) |>
    mutate(dist = dist / 1000)
  fig_variogram <- ggplot() +
    geom_hline(yintercept = sample_var, linetype = "dashed", color = "grey55") +
    geom_line(data = curves, aes(dist, gamma, color = curve, linetype = curve), linewidth = 0.9) +
    geom_point(data = points, aes(dist, gamma, size = np), color = "#08306b") +
    geom_text(data = points, aes(dist, gamma, label = np), vjust = -1.2, size = 2.7) +
    facet_wrap(~ factor(panel, panels), ncol = 3) +
    scale_color_manual(values = c(Spherical = "#cb181d", Exponential = "#fd8d3c",
                                  Gaussian = "#6a51a3", "Returned despite warning" = "grey30"),
                       breaks = c("Spherical", "Exponential", "Gaussian",
                                  "Returned despite warning"), name = NULL) +
    scale_linetype_manual(values = c(Spherical = "solid", Exponential = "solid",
                                     Gaussian = "solid", "Returned despite warning" = "dashed"),
                          breaks = c("Spherical", "Exponential", "Gaussian",
                                     "Returned despite warning"), name = NULL) +
    scale_size_area(max_size = 5, guide = "none") +
    coord_cartesian(xlim = c(0, 11), ylim = c(0, 5)) +
    labs(x = "Distance between sites (km)", y = "Semivariance (cases²)") +
    theme_minimal(base_size = 10) + theme(legend.position = "bottom")
  save_fig(fig_variogram, "fig04_variogram_fits.png", 11, 4.4)

  # Figure: the selected curve over the empirical variogram at other bin widths.
  sel_curve <- curve_of(vfit, "Selected model", "")
  fig_bins <- empirical |>
    filter(cutoff == selected$cutoff) |>
    mutate(bins_label = sprintf("%g km bins", width / 1000),
           pts = map(v, \(x) as_tibble(x)[, c("np", "dist", "gamma")])) |>
    select(bins_label, pts) |>
    unnest(pts) |>
    ggplot(aes(dist / 1000, gamma)) +
    geom_hline(yintercept = sample_var, linetype = "dashed", color = "grey55") +
    geom_line(data = select(sel_curve, -panel), aes(dist, gamma), color = "#cb181d") +
    geom_point(aes(fill = np >= 30, size = np), shape = 21) +
    geom_text(aes(label = np), vjust = -1.1, size = 2.5) +
    facet_wrap(~ bins_label, nrow = 1) +
    scale_fill_manual(values = c(`TRUE` = "#08306b", `FALSE` = "#fd8d3c"),
                      labels = c(`TRUE` = "30 or more pairs", `FALSE` = "fewer than 30 pairs"),
                      name = NULL) +
    scale_size_area(max_size = 4, guide = "none") +
    coord_cartesian(ylim = c(0, 5)) +
    labs(x = "Distance between sites (km)", y = "Semivariance (cases²)") +
    theme_minimal(base_size = 9) + theme(legend.position = "bottom")
  save_fig(fig_bins, "fig05_bin_sensitivity.png", 11, 3.6)

  list(
    sample_var = sample_var,
    empirical = empirical,
    fits = fits,
    tuning = tuning,
    control_model = control_model,
    control_cv = control_cv,
    candidates = candidates,
    selected = selected,
    vfit = vfit,
    v_sel = v_sel,
    selected_name = selected_name
  )
}
