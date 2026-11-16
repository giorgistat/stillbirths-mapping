# stillbirthmap

R package with the data and code for a Bayesian small-area analysis of stillbirth risk
across English local authority districts, prepared for *Paediatric and Perinatal
Epidemiology*. It reanalyses and extends Ranveer Basra's MSc thesis (University of
Birmingham, 2026).

## Installation

```r
# CmdStan is needed for model fitting (once):
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::install_cmdstan()

# the package
remotes::install_github("giorgistat/stillbirths-mapping")
```

## Data

All datasets load with `data()`:

| Dataset | Contents |
|---|---|
| `stillbirths_lad_2024` | 294 modelling units: live births, stillbirths (≥24 weeks) and covariates, 2024 registrations (main analysis) |
| `stillbirths_lad_2025` | Same for 2025 registrations |
| `lad_sf` | Boundaries of the 294 units (`sf`, EPSG:27700) |
| `lad_nb` | Queen adjacency (`spdep::nb`); Isle of Wight linked to Portsmouth, Southampton, New Forest |
| `lad_utla_lookup` | Unit → upper-tier local authority |
| `obesity_pregnancy_utla` | Obesity in early pregnancy (Fingertips 94131), upper-tier authorities, 2022/23–2023/24 |

Modelling units are LADs (May 2024 boundaries), with Cornwall + Isles of Scilly and
Hackney + City of London merged because their data are published combined.

Covariates:
- deprivation: mean LSOA deprivation percentile, IoD 2025;
- % Black and % South Asian: Census 2021, whole population;
- % deliveries to mothers aged 35+, % smoking at time of delivery and % adult obesity: Fingertips 2024/25;
- obesity in early pregnancy: upper-tier values only, modelled at LAD level (see below).

Sources are in `data-raw/sources/`, and `data-raw/make_datasets.R` rebuilds `data/`. All are public
aggregate statistics under the Open Government Licence.

## Model

```r
library(stillbirthmap)
data(stillbirths_lad_2024)
fit <- fit_stillbirth_model(stillbirths_lad_2024)   # main model
fit                                                 # rate ratios
exc <- exceedance(fit)                              # P(rate > England, +10%, +20%, 2.6)
map_exceedance(exc)
```

- **Stillbirth model.** Stillbirths are Poisson with offset log(total births) and IID LAD random
  effects, with a half-normal(0, 1) prior on their SD. BYM2 and Leroux CAR random effects are
  available as sensitivity analyses via `random_effect =`.
- **Obesity in early pregnancy.** If `"pct_obese_preg"` is among the covariates (as in the default,
  `main_covariates()`), obesity in pregnancy is a latent LAD-level variable. It is predicted from
  adult obesity, maternal age, deprivation, ethnicity and smoking, plus an upper-tier random
  effect. It is linked to the published upper-tier values (births-weighted averages with binomial
  error) and fitted jointly with the stillbirth model. Hospital booking data can be added through
  `obesity_extra`.
- **Covariate choice.** All covariates are pre-specified and kept. Cross-validation is used only
  descriptively: `cross_validate()` and `summarise_cv()` compare models, and `drop_one_cv()`
  measures each covariate's predictive contribution (log score and ranked probability score).

## Reproducing the paper

```
R CMD INSTALL .
Rscript analysis/run_analysis_2024.R     # ~45 min on 2 cores; results in output/2024/
```

`legacy/ranveer_msc/` holds the original thesis R Markdown files for reference.
