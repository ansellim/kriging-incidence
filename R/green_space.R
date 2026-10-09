# analyze_green_space: explicit inputs, with figures and tables written through config.
analyze_green_space <- function(data, variograms, surface, config) {
  svy21 <- config$crs
  save_fig <- config$save_fig
  save_tab <- config$save_tab
  sites <- data$sites
  obs <- data$obs
  map_base <- data$map_base
  selected <- variograms$selected
  v_sel <- variograms$v_sel
  sample_var <- variograms$sample_var
  cv_sel <- surface$cv_sel
  cells <- surface$cells
  hull <- surface$hull
  hull_line <- surface$hull_line

  # ---- 4. Green space and case counts -----------------------------------------

  spearman <- cor.test(sites$cases, sites$green_space, method = "spearman", exact = FALSE)

  # Poisson regression for counts; green space per 10 percentage points.
  # No offset: population and observation time are not supplied.
  glm_pois <- glm(cases ~ I(green_space / 10), family = poisson, data = sites)
  dispersion <- sum(residuals(glm_pois, type = "pearson")^2) / df.residual(glm_pois)
  ratio_ci <- exp(confint(glm_pois))[2, ]   # profile-likelihood 95% CI
  lr_test <- anova(glm_pois, test = "Chisq")
  glm_quasi <- glm(cases ~ I(green_space / 10), family = quasipoisson, data = sites)
  quasi_ci <- exp(confint.default(glm_quasi))[2, ]
  greenspace_model <- list(
    ratio = exp(coef(glm_pois))[2], ratio_ci = ratio_ci, lr_p = lr_test$`Pr(>Chi)`[2],
    wald_p = summary(glm_pois)$coefficients[2, 4], dispersion = dispersion,
    quasi_ci = quasi_ci, deviance_explained = 1 - glm_pois$deviance / glm_pois$null.deviance)

  # Moran's I on the counts and on the model residuals (Week 5), under two
  # neighbor definitions, both row-standardized.
  xy <- st_coordinates(obs)
  nb_knn <- make.sym.nb(knn2nb(knearneigh(xy, k = 4)))
  # Smallest band that leaves no site without a neighbor (as in the Week 5 code).
  nn1 <- suppressWarnings(knn2nb(knearneigh(xy, k = 1)))  # warns about sub-graphs; harmless
  band_m <- max(unlist(nbdists(nn1, xy))) * 1.001
  nb_band <- dnearneigh(xy, 0, band_m)
  neighbors <- list("4 nearest neighbors (symmetric)" = nb_knn,
                    "Distance band (smallest with no isolates)" = nb_band)
  pearson_resid <- residuals(glm_pois, type = "pearson")
  moran_table <- imap_dfr(neighbors, \(nb, name) {
    lw <- nb2listw(nb, style = "W")
    raw <- moran.mc(sites$cases, lw, nsim = 9999)
    res <- moran.mc(pearson_resid, lw, nsim = 9999)
    tibble(weights = name, mean_links = mean(card(nb)), isolates = sum(card(nb) == 0),
           cases_I = raw$statistic, cases_p = raw$p.value,
           residual_I = res$statistic, residual_p = res$p.value)
  })
  save_tab(moran_table, "moran_tests.csv")

  # The Week 5 procedure on a linear model of the same covariate: Moran's I for
  # regression residuals and the robust Lagrange multiplier tests (adjRSerr and
  # adjRSlag, the robust RLMerr and RLMlag of the course code, renamed in spdep).
  ols <- lm(cases ~ green_space, data = sites)
  week5_tests <- imap_dfr(neighbors, \(nb, name) {
    lw <- nb2listw(nb, style = "W")
    mt <- lm.morantest(ols, lw)
    rs <- lm.RStests(ols, lw, test = c("adjRSerr", "adjRSlag"))
    tibble(weights = name, ols_I = mt$estimate[[1]], ols_p = mt$p.value,
           robust_error_p = rs$adjRSerr$p.value, robust_lag_p = rs$adjRSlag$p.value)
  })
  save_tab(week5_tests, "week5_ols_tests.csv")

  # Does green space account for the spatial structure? Compare the variogram of
  # the counts with the variogram of residuals after a linear green-space trend.
  v_resid <- variogram(cases ~ green_space, obs, width = selected$width, cutoff = selected$cutoff)
  cv_moran <- moran.mc(cv_sel$residual, nb2listw(nb_knn, style = "W"), nsim = 9999)
  resid_var <- var(residuals(lm(cases ~ green_space, data = sites)))

  curve_glm <- tibble(green_space = seq(min(sites$green_space), max(sites$green_space),
                                        length.out = 100))
  link <- predict(glm_pois, curve_glm, type = "link", se.fit = TRUE)
  curve_glm <- mutate(curve_glm, fit = exp(link$fit), lo = exp(link$fit - 1.96 * link$se.fit),
                      hi = exp(link$fit + 1.96 * link$se.fit))
  fig_green_scatter <- ggplot(sites, aes(green_space, cases)) +
    geom_ribbon(data = curve_glm, aes(green_space, ymin = lo, ymax = hi), inherit.aes = FALSE,
                fill = "#9ecae1", alpha = 0.5) +
    geom_line(data = curve_glm, aes(green_space, fit), color = "#08519c", linewidth = 0.9) +
    geom_jitter(width = 0, height = 0.08, color = "#d7301f", size = 2, alpha = 0.85) +
    labs(x = "Green space around site (%)", y = "Cases",
         subtitle = "A. Poisson regression with 95% confidence band") +
    theme_minimal(base_size = 10)
  fig_green_vario <- bind_rows(
    mutate(as_tibble(v_sel)[, c("dist", "gamma", "np")], series = "Case counts"),
    mutate(as_tibble(v_resid)[, c("dist", "gamma", "np")],
           series = "Residuals after green space")) |>
    ggplot(aes(dist / 1000, gamma, color = series)) +
    geom_hline(yintercept = c(sample_var, resid_var), color = c("#cb181d", "#2171b5"),
               linetype = "dashed") +
    geom_line() + geom_point(aes(size = np)) +
    scale_color_manual(values = c("Case counts" = "#cb181d",
                                  "Residuals after green space" = "#2171b5"), name = NULL) +
    scale_size_area(max_size = 4, guide = "none") +
    coord_cartesian(ylim = c(0, 3.5)) +
    labs(x = "Distance between sites (km)", y = "Semivariance (cases²)",
         subtitle = "B. Variogram before and after green space") +
    theme_minimal(base_size = 10) + theme(legend.position = "bottom")
  save_fig(fig_green_scatter | fig_green_vario, "fig10_greenspace_association.png", 11, 4.4)

  fig_green_overlay <- ggplot() + map_base() +
    geom_sf(data = filter(cells, in_hull), aes(fill = pred, color = after_scale(fill)),
            linewidth = 0.05) +
    hull_line +
    geom_sf(data = obs, size = 4.2, color = "black") +
    geom_sf(data = obs, aes(color = green_space), size = 3.2) +
    scale_fill_distiller(palette = "YlOrRd", direction = 1, name = "Predicted\ncases") +
    scale_color_distiller(palette = "Greens", direction = 1, name = "Green\nspace (%)") +
    coord_sf(xlim = st_bbox(hull)[c(1, 3)] + c(-1500, 1500),
             ylim = st_bbox(hull)[c(2, 4)] + c(-1500, 1500))
  save_fig(fig_green_overlay, "fig11_greenspace_overlay.png", 7.5, 5)

  list(
    spearman = spearman,
    glm_pois = glm_pois,
    greenspace_model = greenspace_model,
    moran_table = moran_table,
    band_m = band_m,
    v_resid = v_resid,
    resid_var = resid_var,
    week5_tests = week5_tests,
    cv_moran = cv_moran
  )
}
