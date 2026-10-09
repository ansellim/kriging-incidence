# kriging-incidence

Reproducible spatial analysis of disease case counts at 30 sites in Singapore for SPH5105 Individual Assessment 2. The analysis selects a variogram, predicts counts on a 1 km grid, profiles high-risk zones, examines the green-space association, and proposes two intervention centers.

**These are case counts, not incidence rates.** Population denominators, observation time, and site catchment areas were not supplied. The report treats conversions to cases per square kilometer as illustrative support assumptions.

Read the [PDF report](reports/report.pdf), edit the [Quarto source](report.qmd), or inspect the [Word export](reports/report.docx). Cached inputs and computed outputs are included, so the report can be rendered without downloading data or rerunning the simulations.

## Repository layout

| Path | Purpose |
|---|---|
| `run_analysis.R` | Entry point; calls the stages with explicit inputs and outputs. |
| `R/setup.R` | Dependency checks, paths, seed, and output writers. |
| `R/data.R` | Input validation, coordinate checks, spatial layers, and exploration. |
| `R/variograms.R`, `R/variogram_helpers.R` | Variogram tuning, convergence checks, and cross-validation. |
| `R/kriging.R`, `R/zone_helpers.R` | Ordinary kriging, zones, block sensitivity, and conditional simulation. |
| `R/green_space.R` | Poisson models, Moran's I, and spatial regression diagnostics. |
| `R/facilities.R`, `R/facility_helpers.R` | Community-club filtering and exhaustive two-center optimization. |
| `R/sensitivity.R`, `R/sensitivity_helpers.R` | Alternative variograms, influential-site removal, and universal kriging. |
| `R/export.R` | Saves report inputs, grid predictions, and session metadata. |
| `data/` | Site counts, cached government layers, and source documentation. |
| `outputs/` | Twelve analysis figures, fourteen tables, and report results. |
| `report.qmd`, `reports/`, `references/` | Report source, rendered documents, bibliography, citation style, and filter. |
| `tests/` | Independent scientific checks and exact table reproduction against the audited baseline. |
| `renv.lock`, `scripts/` | Pinned R dependencies and installation script. |

## Reproduce

Use R 4.4.1 (the audited version), Python 3, Quarto, and XeLaTeX. The spatial packages need GDAL, GEOS, PROJ, and UDUNITS system libraries. On Ubuntu:

```sh
sudo apt install gdal-bin libgdal-dev libgeos-dev libproj-dev libudunits2-dev
```

From the repository root:

```sh
Rscript scripts/install_dependencies.R
Rscript run_analysis.R
Rscript tests/check_outputs.R
python3 tests/check_reproduction.py
quarto render report.qmd
```

The pipeline takes about one minute in the original environment, uses seed 5105, and runs offline once R dependencies are installed. Quarto writes PDF and Word documents to `reports/`. To render the included results only, run the last command. For the PDF, install a LaTeX distribution with XeLaTeX, for example `quarto install tinytex`.

The R checks independently verify 63 results, including all 1,017 kriging predictions and variances, every empirical variogram setting, cross-validation metrics, Poisson confidence limits, spatial statistics, zone totals, and all 496 candidate-center scores. The Python check compares all fourteen CSV tables with SHA-256 hashes from the audited analysis before modularization. Exact equality assumes the recorded package versions and numerical environment; changes in spatial libraries can affect geometry or floating-point details.

## Main results and limitations

The selected spherical variogram has a nugget of 1.47 cases squared, partial sill of 1.09 cases squared, and range of 7.05 km. Leave-one-out RMSE is 1.40 versus 1.75 for the pure-nugget control. The top-quarter zone covers about 55.7 square kilometers. Braddell Heights and Kampong Glam Community Clubs minimize the case-weighted distance among the 496 shortlisted pairs.

The variogram was held fixed during cross-validation, and its selection reused those errors. Gaussian kriging is an approximation for sparse counts. Most land is outside the sample hull. Green-space associations are ecological and do not establish causality. Facility distances are straight-line distances; capacity, opening hours, walking routes, and operator agreement need separate verification. See the report for the full assumptions and equations.

## Data attribution

Contains information from URA's *Master Plan 2019 Subzone Boundary (No Sea)*, the People's Association's *Community Clubs (2024)*, and LTA's *MRT Station Exits*, accessed on October 7, 2026 from [data.gov.sg](https://data.gov.sg), made available under the [Singapore Open Data Licence](https://data.gov.sg/open-data-licence). See [data/README.md](data/README.md) for dataset links, fields, and limitations. This project has no official status or agency endorsement.

The 30-site table comes from the SPH5105 assessment. The assignment question, lecture slides, textbooks, and session handoffs are not included. The Chicago citation style retains its Creative Commons attribution and share-alike notice in [references/chicago-author-date.csl](references/chicago-author-date.csl). No separate license is granted here for the assessment table, code, or report.

Claude Code (Anthropic) and Codex (OpenAI) assisted with the code and report. The author takes responsibility for the work, as stated in the report's AI declaration.
