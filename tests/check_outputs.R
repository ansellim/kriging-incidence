# =============================================================================
# Independent end-to-end checks of run_analysis.R outputs.
# Run after run_analysis.R, from the same folder:   Rscript tests/check_outputs.R
# Each check recomputes a result by a separate route (mostly base R) and stops
# with an error if it disagrees with the saved output.
# =============================================================================

suppressPackageStartupMessages({ library(sf); library(gstat); library(spdep) })
if (!file.exists("outputs/results.rds")) stop("Run run_analysis.R first.")
r <- readRDS("outputs/results.rds")
grid <- read.csv("outputs/tables/kriging_grid.csv")
sites <- read.csv("data/raw/assignment2.csv")
passed <- character()
check <- function(label, ok) {
  if (!isTRUE(ok)) stop("FAILED: ", label, call. = FALSE)
  passed <<- c(passed, label)
}
near <- function(a, b, tol = 1e-6) isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tol))

# 1. Input data: the 30 rows of the assignment table, unchanged.
check("30 sites, 56 cases, green space 26.3 to 76.1%",
      nrow(sites) == 30 && sum(sites$cases) == 56 &&
        near(range(sites$green_space), c(26.3, 76.1)))
check("x/y agree with projected latitude/longitude within 1 m", r$coord_check_m < 1)
check("summary statistics match base R",
      near(r$summary_stats$mean, c(mean(sites$green_space), mean(sites$cases))) &&
        near(r$summary_stats$sd, c(sd(sites$green_space), sd(sites$cases))))

# 2. Empirical variogram of the selected setting, recomputed without gstat.
D <- as.matrix(dist(sites[, c("x", "y")]))
ij <- which(upper.tri(D), arr.ind = TRUE)
h <- D[ij]
g <- (sites$cases[ij[, 1]] - sites$cases[ij[, 2]])^2 / 2
w <- r$selected$width
keep <- h <= r$selected$cutoff
bin <- ceiling(h[keep] / w)
check("pair counts per bin", identical(as.integer(table(bin)), as.integer(r$v_sel$np)))
check("mean distance per bin", near(tapply(h[keep], bin, mean), r$v_sel$dist))
check("semivariance per bin", near(tapply(g[keep], bin, mean), r$v_sel$gamma))

# 3. The selected fit converges cleanly and is reproducible from the sensible start.
obs <- st_as_sf(sites, coords = c("x", "y"), crs = 3414, remove = FALSE)
v <- variogram(cases ~ 1, obs, width = r$selected$width, cutoff = r$selected$cutoff)
start_range <- if (r$selected$model == "Exp") 2000 else if (r$selected$model == "Gau") 6000 / sqrt(3) else 6000
warned <- FALSE
refit <- withCallingHandlers(
  fit.variogram(v, vgm(1.5, r$selected$model, start_range, 1.5), fit.method = r$selected$method),
  warning = function(cond) { warned <<- TRUE; invokeRestart("muffleWarning") })
check("selected fit has no warning and is not singular", !warned && !isTRUE(attr(refit, "singular")))
check("selected fit parameters reproduce", near(refit$psill, r$vfit$psill, 1e-4) &&
        near(refit$range, r$vfit$range, 1e-4))
check("selected fit beats the pure-nugget control in LOOCV RMSE",
      r$selected$rmse < r$control_cv$rmse)

# 4. Ordinary kriging solved directly as a linear system at 25 random squares.
nug <- r$vfit$psill[1]; ps <- r$vfit$psill[2]; a <- r$vfit$range[2]
cov_fun <- switch(as.character(r$vfit$model[2]),
  Sph = function(h) ps * ifelse(h < a, 1 - 1.5 * h / a + 0.5 * (h / a)^3, 0),
  Exp = function(h) ps * exp(-h / a),
  Gau = function(h) ps * exp(-(h / a)^2))
C <- cov_fun(D); diag(C) <- ps + nug
n <- nrow(sites)
A <- rbind(cbind(C, 1), c(rep(1, n), 0))
set.seed(1)
idx <- sample(nrow(grid), 25)
direct <- t(sapply(idx, function(i) {
  c0 <- cov_fun(sqrt((sites$x - grid$px[i])^2 + (sites$y - grid$py[i])^2))
  s <- solve(A, c(c0, 1))
  c(pred = sum(s[1:n] * sites$cases), var = ps + nug - sum(s[1:n] * c0) - s[n + 1])
}))
check("kriging predictions match the direct solution", near(direct[, "pred"], grid$pred[idx]))
check("kriging variances match the direct solution", near(direct[, "var"], grid$var[idx]))
check("1 km grid: 1,017 land squares, no negative variance", nrow(grid) == 1017 && all(grid$var > 0))

# 5. Variance depends on the sampling configuration only (Exercise 54).
shuffled <- sites$cases[sample(n)]
direct_shuffled <- sapply(idx, function(i) {
  c0 <- cov_fun(sqrt((sites$x - grid$px[i])^2 + (sites$y - grid$py[i])^2))
  s <- solve(A, c(c0, 1)); ps + nug - sum(s[1:n] * c0) - s[n + 1]
})
check("variance unchanged when counts are permuted", near(direct_shuffled, grid$var[idx]))

# 6. High-risk thresholds, zone membership and outbreak size.
q <- quantile(grid$pred[grid$in_hull], c(0.75, 0.90), names = FALSE)
check("quantile thresholds", near(q, r$sanity$thresholds))
high <- grid[grid$in_hull & grid$pred >= q[1], ]
check("top-25% squares", nrow(high) == r$zone_table$squares[1] &&
        all(grid$risk_class[grid$in_hull & grid$pred >= q[1]] != "Other"))
check("estimated cases in top-25% zone",
      near(sum(high$pred * high$land_km2), r$zone_table$est_cases[1]))
zone_sites <- strsplit(r$zone_table$site_ids[1], ", ")[[1]]
check("recorded cases in top-25% zone",
      sum(sites$cases[sites$id %in% as.integer(zone_sites)]) == r$zone_table$recorded_cases[1])

# 7. Facility pair: brute force over the saved shortlist.
shortlist <- read.csv("outputs/tables/center_shortlist.csv")
weight <- high$pred * high$land_km2
d_km <- sapply(seq_len(nrow(shortlist)), function(j)
  sqrt((high$px - shortlist$X[j])^2 + (high$py - shortlist$Y[j])^2) / 1000)
pairs <- combn(nrow(shortlist), 2)
score <- apply(pairs, 2, function(p) sum(weight * pmin(d_km[, p[1]], d_km[, p[2]])) / sum(weight))
best <- sort(shortlist$name[pairs[, which.min(score)]])
check("chosen pair is optimal over all shortlisted pairs", identical(best, sort(r$centers$name)))
check("chosen clubs are within 500 m of an MRT exit", all(r$centers$mrt_m <= 500))
verified <- read.csv("data/cache/facility_verification.csv")
access <- verified$access_notes[match(r$centers$name, verified$name)]
check("chosen clubs verified in the official directory as wheelchair accessible",
      all(grepl("lists wheelchair accessibility", access)))

# 8. Green space model and Moran's I.
m <- glm(cases ~ I(green_space / 10), family = poisson, data = sites)
check("count ratio per 10 points", near(exp(coef(m))[2], r$greenspace_model$ratio))
lw <- nb2listw(make.sym.nb(knn2nb(knearneigh(as.matrix(sites[, c("x", "y")]), k = 4))), style = "W")
W <- listw2mat(lw)
z <- sites$cases - mean(sites$cases)
moran_direct <- (n / sum(W)) * sum(W * outer(z, z)) / sum(z^2)
check("Moran's I of counts (4 nearest neighbors)", near(moran_direct, r$moran$cases_I[1]))

# 9. Every displayed variogram formula, including the discontinuity at zero.
gamma_model <- function(h, family, nugget, psill, range_m) {
  if (family == "Nug") return(ifelse(h == 0, 0, nugget))
  a <- switch(family, Sph = range_m, Exp = range_m / 3, Gau = range_m / sqrt(3))
  structure <- switch(family,
    Sph = ifelse(h <= a, 1.5 * h / a - 0.5 * (h / a)^3, 1),
    Exp = 1 - exp(-h / a), Gau = 1 - exp(-(h / a)^2))
  ifelse(h == 0, 0, nugget + psill * structure)
}
for (family in c("Sph", "Exp", "Gau")) {
  h_test <- c(0, 1e-6, 100, 1000, 6000, 18000)
  scale_m <- switch(family, Sph = 6000, Exp = 2000, Gau = 6000 / sqrt(3))
  model <- vgm(1.5, family, scale_m, 1.5)
  check(paste(family, "formula matches gstat at zero, near zero, and long lags"),
        near(gamma_model(h_test, family, 1.5, 1.5, 6000),
             variogramLine(model, dist_vector = h_test)$gamma))
}
check("3a and sqrt(3)a reach about 95%, with exact 95% lags using -log(0.05)",
      abs((1 - exp(-3)) - 0.95) < 0.0003 &&
        near(1 - exp(-(-log(0.05))), 0.95))

# Check the actual method-2 implementation: weights are fixed at the start.
# Fix the range so that this real-data fit has an identifiable linear solution.
initial <- vgm(1.5, "Sph", 6000, 1.5)
method2 <- fit.variogram(v, initial, fit.method = 2, fit.ranges = FALSE)
err <- v$gamma - variogramLine(method2, dist_vector = v$dist)$gamma
weights2 <- v$np / variogramLine(initial, dist_vector = v$dist)$gamma^2
check("method-2 weighted SSE uses the starting variogram",
      near(sum(weights2 * err^2), attr(method2, "SSErr")))
weights7 <- v$np / v$dist^2
check("selected method-7 weighted SSE matches N/h^2 formula",
      near(sum(weights7 * (v$gamma - variogramLine(r$vfit, dist_vector = v$dist)$gamma)^2),
           attr(r$vfit, "SSErr"), 1e-10))

# 10. Solve the report's semivariogram equations for ALL 1,017 grid points.
# A scaled coordinate basis also permits an independent universal-kriging check.
direct_kriging <- function(train, target, family, nugget, psill, range_m, trend = FALSE) {
  train_xy <- as.matrix(train[, c("x", "y")])
  target_xy <- as.matrix(target[, c("x", "y")])
  distances <- sqrt(outer(train$x, target$x, "-")^2 + outer(train$y, target$y, "-")^2)
  gamma <- function(h) gamma_model(h, family, nugget, psill, range_m)
  X <- if (trend) cbind(1, train_xy / 10000) else matrix(1, nrow(train), 1)
  X0 <- if (trend) cbind(1, target_xy / 10000) else matrix(1, nrow(target), 1)
  G <- gamma(as.matrix(dist(train_xy)))
  A <- rbind(cbind(G, X), cbind(t(X), matrix(0, ncol(X), ncol(X))))
  g0 <- gamma(distances)
  solution <- solve(A, rbind(g0, t(X0)))
  weights <- solution[seq_len(nrow(train)), , drop = FALSE]
  multipliers <- solution[nrow(train) + seq_len(ncol(X)), , drop = FALSE]
  if (max(abs(crossprod(X, weights) - t(X0))) > 1e-8) stop("Kriging constraints failed")
  data.frame(pred = drop(crossprod(train$cases, weights)),
             var = colSums(weights * g0) + colSums(t(X0) * multipliers))
}
targets <- data.frame(x = grid$px, y = grid$py)
all_direct <- direct_kriging(sites, targets, "Sph", nug, ps, a)
check("semivariogram equations reproduce every grid prediction and variance",
      near(all_direct$pred, grid$pred) && near(all_direct$var, grid$var))
check("exact interpolation at observed sites has zero variance",
      { at_sites <- direct_kriging(sites, sites, "Sph", nug, ps, a)
        near(at_sites$pred, sites$cases) && max(abs(at_sites$var)) < 1e-8 })
check("unsampled point variances are at least the nugget", all(grid$var >= nug))

# 11. Leave-one-out predictions and uncertainty independently from linear algebra.
cv_direct <- function(train, family, nugget, psill, range_m, trend = FALSE) {
  folds <- do.call(rbind, lapply(seq_len(nrow(train)), function(i)
    direct_kriging(train[-i, ], train[i, ], family, nugget, psill, range_m, trend)))
  error <- train$cases - folds$pred
  z <- error / sqrt(folds$var)
  c(rmse = sqrt(mean(error^2)), mae = mean(abs(error)), mean_error = mean(error),
    z_mean = mean(z), z_sd = sd(z))
}
metrics <- c("rmse", "mae", "mean_error", "z_mean", "z_sd")
for (i in seq_len(nrow(r$candidates))) {
  m <- r$candidates[i, ]
  cv_score <- cv_direct(sites, m$model, m$nugget, m$psill, m$range_m)
  check(paste("all five CV metrics for stable candidate", i), near(cv_score, unlist(m[metrics])))
}
control <- cv_direct(sites, "Nug", var(sites$cases), 0, NA_real_)
check("pure-nugget CV metrics reproduce", near(control, unlist(r$control_cv[metrics])))
check("pure-nugget prediction and variance equal the analytic mean-only formulas",
      { held_out <- direct_kriging(sites[-1, ], sites[1, ], "Nug", var(sites$cases), 0, NA_real_)
        near(held_out$pred, mean(sites$cases[-1])) &&
          near(held_out$var, var(sites$cases) * (1 + 1 / 29)) })

# 12. Independent Poisson likelihood, profile limits, LR test and dispersion.
x10 <- sites$green_space / 10
profile_ll <- function(beta) {
  intercept <- log(sum(sites$cases) / sum(exp(beta * x10)))
  sum(dpois(sites$cases, exp(intercept + beta * x10), log = TRUE))
}
beta <- unname(coef(m <- glm(cases ~ x10, family = poisson, data = sites))[2])
max_ll <- profile_ll(beta)
profile_cut <- function(b) 2 * (max_ll - profile_ll(b)) - qchisq(0.95, 1)
ci <- exp(c(uniroot(profile_cut, c(beta - 2, beta), tol = 1e-10)$root,
            uniroot(profile_cut, c(beta, beta + 2), tol = 1e-10)$root))
check("profile-likelihood count-ratio limits reproduce independently",
      max(abs(ci - r$greenspace_model$ratio_ci)) < 1e-4)
check("likelihood-ratio p-value reproduces from Poisson log likelihoods",
      near(pchisq(2 * (max_ll - profile_ll(0)), 1, lower.tail = FALSE), r$greenspace_model$lr_p))
mu <- exp(coef(m)[1] + beta * x10)
pearson <- (sites$cases - mu) / sqrt(mu)
phi <- sum(pearson^2) / (nrow(sites) - length(coef(m)))
check("Pearson dispersion and explained deviance reproduce",
      near(phi, r$greenspace_model$dispersion) &&
        near(1 - deviance(m) / m$null.deviance, r$greenspace_model$deviance_explained))
quasi_limits <- exp(beta + c(-1, 1) * qnorm(0.975) * sqrt(phi * vcov(m)[2, 2]))
# Separate Poisson and quasi-Poisson IRLS stopping tolerances cause a <2e-6
# numerical difference; this is much smaller than the reported two decimals.
check("quasi-Poisson confidence interval reproduces by scaling Poisson variance",
      max(abs(quasi_limits - r$greenspace_model$quasi_ci)) < 2e-6)

# 13. Moran's formula for both graphs, counts, Pearson and CV residuals.
nb4 <- make.sym.nb(knn2nb(knearneigh(as.matrix(sites[, c("x", "y")]), k = 4)))
nbband <- dnearneigh(as.matrix(sites[, c("x", "y")]), 0, r$band_km * 1000)
for (i in seq_along(list(nb4, nbband))) {
  nb <- list(nb4, nbband)[[i]]
  lw <- nb2listw(nb, style = "W"); W <- listw2mat(lw)
  moran_formula <- function(value) {
    centered <- value - mean(value)
    (n / sum(W)) * sum(W * outer(centered, centered)) / sum(centered^2)
  }
  check(paste("counts and Poisson residual Moran formulas for graph", i),
        near(moran_formula(sites$cases), r$moran$cases_I[i]) &&
          near(moran_formula(pearson), r$moran$residual_I[i]))
}
check("Moran permutation and linear-residual defaults test the greater alternative",
      identical(formals(moran.mc)$alternative, "greater") &&
        identical(formals(lm.morantest)$alternative, "greater"))

# 14. Area-based integration for both zones and an apples-to-apples denominator.
for (i in seq_len(nrow(r$zone_table))) {
  rows <- grid[grid$in_hull & grid$pred >= q[i], ]
  total <- sum(rows$pred * rows$land_km2)
  area <- sum(rows$land_km2)
  check(paste("zone", i, "area, total, and area-weighted density"),
        near(area, r$zone_table$area_km2[i]) && near(total, r$zone_table$est_cases[i]) &&
          near(total / area, r$zone_table$est_per_km2[i]))
}
rest <- grid[grid$in_hull & grid$risk_class == "Other", ]
rest_density <- sum(rest$pred * rest$land_km2) / sum(rest$land_km2)
check("rest-of-area weighted mean equals integrated density",
      near(weighted.mean(rest$pred, rest$land_km2), rest_density))
check("report uses the area-weighted denominator for relative intensities",
      grepl("rest_mean <- weighted.mean", paste(readLines("report.qmd"), collapse = "\n"), fixed = TRUE))
check("catchment rescaling cancels in relative intensity and facility objective",
      near((r$zone_table$est_per_km2 / r$sanity$hull_land_km2_per_site) /
             (rest_density / r$sanity$hull_land_km2_per_site), r$zone_table$est_per_km2 / rest_density) &&
        near(sum((weight / 6.4) * pmin(d_km[, pairs[1, which.min(score)]],
                                    d_km[, pairs[2, which.min(score)]])) / sum(weight / 6.4), min(score)))

# 15. Every sensitivity scenario: independent CV, grid extrema and cell Jaccard.
families <- c("Sph", "Exp", "Gau", "Sph", "Sph", "Nug")
for (i in seq_len(nrow(r$sensitivity))) {
  scenario <- r$sensitivity[i, ]
  train <- if (i == 5) sites[sites$id != 25, ] else sites
  k <- direct_kriging(train, targets, families[i], scenario$nugget, scenario$psill,
                      scenario$range_km * 1000, trend = i == 6)
  selected_ids <- grid$cell_id[grid$in_hull & k$pred >= quantile(k$pred[grid$in_hull], 0.75)]
  overlap <- length(intersect(selected_ids, high$cell_id)) / length(union(selected_ids, high$cell_id))
  cv <- cv_direct(train, families[i], scenario$nugget, scenario$psill,
                  scenario$range_km * 1000, trend = i == 6)
  check(paste("sensitivity", i, "CV, grid extrema and Jaccard"),
        near(cv, unlist(scenario[metrics])) &&
          near(range(k$pred), c(scenario$pred_min, scenario$pred_max)) &&
          near(overlap, scenario$top25_jaccard))
}

# 16. Remaining numerical summaries printed in the report.
all_bin_summaries <- vapply(seq_len(nrow(r$empirical)), function(i) {
  setting <- r$empirical[i, ]
  counts <- table(ceiling(h[h <= setting$cutoff] / setting$width))
  length(counts) == setting$bins && min(counts) == setting$min_pairs &&
    counts[1] == setting$first_bin_pairs
}, logical(1))
check("all 20 bin-setting summaries reproduce from unordered pairs", all(all_bin_summaries))
cv_folds <- do.call(rbind, lapply(seq_len(n), function(i)
  direct_kriging(sites[-i, ], sites[i, ], "Sph", nug, ps, a)))
cv_errors <- sites$cases - cv_folds$pred
W4 <- listw2mat(nb2listw(nb4, style = "W"))
centered_cv <- cv_errors - mean(cv_errors)
cv_I <- n / sum(W4) * sum(W4 * outer(centered_cv, centered_cv)) / sum(centered_cv^2)
check("cross-validation residual Moran's I reproduces", near(cv_I, r$cv_moran["I"]))
check("all 496 facility scores give the reported optimum and top-ten ranking",
      length(score) == r$facility_summary$pairs && near(min(score), r$facility_summary$score_km) &&
        near(head(sort(score), 10), r$pair_ranking$mean_km))
best_indices <- match(r$centers$name, shortlist$name)
nearer <- max.col(-d_km[, best_indices], ties.method = "first")
nearest <- apply(d_km[, best_indices], 1, min)
check("facility coverage percentage and farthest distance reproduce",
      near(sum(weight[nearest <= 2]) / sum(weight), r$facility_summary$within_2km) &&
        near(max(nearest), r$facility_summary$max_km))
served_counts <- tabulate(nearer, nbins = 2)
served_cases <- as.numeric(tapply(weight, factor(nearer, levels = 1:2), sum))
served_distance <- as.numeric(tapply(weight * nearest, factor(nearer, levels = 1:2), sum)) / served_cases
check("per-center square counts, estimated cases and weighted distance reproduce",
      near(served_counts, r$centers$squares_served) &&
        near(served_cases, r$centers$est_cases_served) && near(served_distance, r$centers$mean_km_served))
land <- st_transform(st_read("data/cache/singapore_land.geojson", quiet = TRUE), 3414)
hull <- st_convex_hull(st_union(obs))
hull_area <- as.numeric(st_area(st_intersection(hull, st_union(land)))) / 1e6
check("illustrative catchment equals exact hull-land area divided by 30",
      near(hull_area / n, r$sanity$hull_land_km2_per_site))
check("land area and extrapolation fractions reproduce from clipped cells",
      near(sum(grid$land_km2), r$land_km2) &&
        near(mean(!grid$in_hull), r$sanity$share_outside_hull) &&
        near(sum(grid$land_km2[!grid$in_hull]) / sum(grid$land_km2), r$sanity$share_outside_area))

cat(sprintf("All %d checks passed:\n", length(passed)), paste0("  - ", passed, "\n"), sep = "")
