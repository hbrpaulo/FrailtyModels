# FrailtyModels

This repository contains R scripts and an R Markdown notebook for exploring frailty models in survival analysis. The code currently focuses on Gamma frailty with parametric baseline hazards and demonstrates estimation and simulation strategies.

While this repository begins with standalone scripts, the long‑term goal is to turn these components into a comprehensive R package. The planned package will support:

- **Parametric frailty models** with flexible baseline choices (Weibull, exponential, Gompertz, gamma, lognormal, loglogistic, piecewise exponential, and others).
- The **Cox proportional hazards model** with frailty terms.
- **Shape and rate parameterization** for baseline hazards to simplify prior specification and interpretation.
- A variety of frailty distributions, starting with Gamma but open to future extensions.

## Contents

- **`initialScript.R`** – A prototype script that fits parametric frailty models via maximum likelihood and returns empirical Bayes frailty estimates for a user-supplied data set.
- **`frailtyWeibullEstimation.Rmd`** – An R Markdown notebook walking through simulation, estimation, and coverage analysis for a Gamma frailty model with a Weibull baseline.
- **`scripts/R/legacy/analysis/`** – Modular R functions for simulating data,
  estimating parameters, computing coverage, and visualizing results for the
  Weibull baseline with Gamma frailty.
- **`usage_examples/vignette_weibull.Rmd`** – Reproducible walkthrough sourcing
  the analysis scripts and demonstrating the full Weibull frailty workflow.
- **`LICENSE`** – MIT License.

## Requirements

These examples rely on R packages such as `survival`, `tidyverse`, `parfm`, `kableExtra`, and `ggplot2`. Install them from CRAN before running the code:

```r
install.packages(c("survival", "tidyverse", "parfm", "kableExtra", "ggplot2"))
```

## Usage

1. Edit `initialScript.R` to point to your own CSV file with survival data.
2. Source the script or run it line by line in R to inspect the fitted parameters and frailty estimates.
3. Open `frailtyWeibullEstimation.Rmd` in RStudio (or a similar environment) and knit it to HTML to reproduce the simulation and coverage study.
4. Alternatively, source the functions under `scripts/R/legacy/analysis/` and
   run `usage_examples/vignette_weibull.Rmd` to explore the modular Weibull
   frailty workflow including simulation, estimation, coverage studies, and
   plotting.

Example data are not included, so you will need to supply your own data set to reproduce the analyses.

## Future work

Development will continue toward a full R package that exposes user-friendly functions for fitting parametric and semiparametric frailty models. The package will emphasize modularity, allowing new baseline hazards or frailty distributions to be added easily.

## License

This project is distributed under the terms of the MIT License. See the [LICENSE](LICENSE) file for details.
