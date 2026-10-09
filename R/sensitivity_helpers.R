# Reusable numerical helpers. Distances are in EPSG:3414 meters.

pick_fit <- function(candidates, model_code, width, cutoff, method) {
  row <- filter(candidates, model == model_code, .data$width == .env$width,
                .data$cutoff == .env$cutoff, .data$method == .env$method)
  stopifnot(nrow(row) == 1)
  row$model_fit[[1]]
}

fit_or_nugget <- function(v, model) {
  out <- fit_one(v, model, 1.5, 1.5, 6000, 7)
  if (out$clean) return(out$model_fit[[1]])
  message("No clean ", model, " fit (", out$warning, "); using a fitted pure nugget.")
  fit.variogram(v, vgm(psill = 1.5, model = "Nug", range = 0))
}
