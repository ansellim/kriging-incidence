#!/usr/bin/env Rscript
# Restore the versions recorded in renv.lock, including transitive dependencies.
if (!file.exists("renv.lock")) stop("Run from the kriging-incidence repository root.")
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = "https://cloud.r-project.org")
}
renv::restore(project = ".", prompt = FALSE)
renv::status(project = ".")
