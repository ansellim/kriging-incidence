# Reusable numerical helpers. Distances are in EPSG:3414 meters.

practical_range <- function(model, range) {
  switch(model, Sph = range, Exp = 3 * range, Gau = sqrt(3) * range)
}

fit_one <- function(v, model, nugget, psill, start_range_m, method) {
  start <- vgm(psill = psill, model = model, nugget = nugget,
               range = start_range_m / practical_range(model, 1))
  messages <- character()
  fit <- withCallingHandlers(
    fit.variogram(v, start, fit.method = method),
    warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  tibble(nugget = fit$psill[1], psill = fit$psill[2],
         range_m = practical_range(model, fit$range[2]),
         singular = isTRUE(attr(fit, "singular")),
         warning = paste(unique(messages), collapse = "; "),
         model_fit = list(fit)) |>
    mutate(clean = !singular & warning == "")
}

loocv <- function(model, formula = cases ~ 1, data) {
  cv <- krige.cv(formula, data, model = model, nfold = nrow(data), verbose = FALSE)
  tibble(rmse = sqrt(mean(cv$residual^2)), mae = mean(abs(cv$residual)),
         mean_error = mean(cv$residual), z_mean = mean(cv$zscore), z_sd = sd(cv$zscore))
}

curve_of <- function(model, label, panel) {
  line <- variogramLine(model, maxdist = 11000, n = 300)
  tibble(dist = line$dist / 1000, gamma = line$gamma, curve = label, panel = panel)
}
