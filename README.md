# stillbirthmap

R package with the data and code to reproduce the analysis in

> Giorgi E, Basra R, Heazell AEP. *Mapping stillbirth risk across local authorities in England
> by integrating multiple national data sources.* Manuscript in preparation for *Paediatric and
> Perinatal Epidemiology*.

The analysis estimates stillbirth rates for the 294 local authority districts (LADs) of England
using 2024 birth registrations. It combines several national data sources published at
different geographical levels in a single statistical model. Exceedance probabilities are then
used to identify districts whose stillbirth rate is likely to be above national benchmarks. The
work builds on Ranveer Basra's MSc thesis (University of Birmingham, 2026).

## Installation

Model fitting uses Stan through `cmdstanr`, which needs a one-off installation of CmdStan:

```r
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::install_cmdstan()

remotes::install_github("giorgistat/stillbirths-mapping")
```

## Data

All datasets are public aggregate statistics released under the Open Government Licence. They
load with `data()`:

| Dataset | Contents |
|---|---|
| `stillbirths_lad_2024` | 294 LADs: live births, stillbirths (24+ weeks) and covariates, 2024 registrations (main analysis) |
| `stillbirths_lad_2025` | Same for 2025 registrations |
| `lad_sf` | LAD boundaries (`sf`, EPSG:27700) |
| `lad_nb` | Queen adjacency (`spdep::nb`); Isle of Wight linked to Portsmouth, Southampton and New Forest |
| `lad_utla_lookup` | LAD to upper-tier local authority lookup |
| `obesity_pregnancy_utla` | Obesity in early pregnancy (OHID Fingertips 94131), upper-tier authorities, 2022/23 and 2023/24 |

LADs follow the May 2024 boundaries. Cornwall and the Isles of Scilly, and Hackney and the City of
London, are merged because their data are published combined.

Stillbirths and births come from ONS birth registrations. The covariates are:

| Covariate | Source | Level published |
|---|---|---|
| Deprivation (mean LSOA percentile) | English Indices of Deprivation 2025 | LSOA, averaged to LAD |
| % Black, % South Asian (whole population) | Census 2021 | LAD |
| % deliveries to mothers aged 35+ | OHID Fingertips, 2024/25 | LAD |
| % smoking at time of delivery | OHID Fingertips, 2024/25 | LAD |
| % adult obesity (obesity sub-model and sensitivity analysis) | OHID Fingertips (Active Lives Survey), 2024/25 | LAD |
| % women with BMI 30+ at the booking appointment | OHID Fingertips 94131 | Upper-tier authority, about a third withheld |

The raw files are in `data-raw/sources/`, and `data-raw/make_datasets.R` rebuilds `data/`.

## Model

Stillbirth counts are analysed with a Poisson generalised linear mixed model. The log of total
births is used as an offset, and each LAD has an independent normal random effect.

Obesity in early pregnancy is published only for upper-tier authorities, and about a third of
those values are withheld. It is therefore treated as a latent LAD-level variable, predicted from
adult obesity, maternal age, deprivation, ethnicity and smoking, with an upper-tier random effect.
The published values are modelled as births-weighted averages of the latent LAD values, with
binomial sampling error, extra measurement error and a shift between the two years. This
measurement model is fitted jointly with the stillbirth model, so uncertainty about obesity is
carried through to all outputs. The Stan code is in `inst/stan/`.

All covariates are pre-specified and kept in the model. Ten-fold cross-validation, scored with
the log score and the ranked probability score, is used descriptively to compare models and to
measure each covariate's predictive contribution.

Exceedance probabilities are computed for four thresholds:

- the England rate;
- England +10%;
- England +20%;
- the national ambition of 2.6 per 1,000 births.

LADs are classified as high risk if the probability is above 0.8 and low risk if it is below 0.2.

BYM2 and Leroux spatial random effects are available through `random_effect =`. They did not
improve the fit and are reported only in the supplementary material.

```r
library(stillbirthmap)
data(stillbirths_lad_2024)
fit <- fit_stillbirth_model(stillbirths_lad_2024)   # main model (joint obesity sub-model)
fit                                                 # rate ratios
exc <- exceedance(fit)                              # exceedance probabilities and classes
map_exceedance(exc)
```

## Reproducing the paper

```
R CMD INSTALL .
Rscript analysis/run_analysis_2024.R     # about 45 minutes on 2 cores
```

Results are written to `output/2024/`: rate ratios, model parameters, exceedance probabilities
for every LAD, latent obesity estimates, cross-validation summaries, residual checks and
figures.

`legacy/ranveer_msc/` holds the original thesis R Markdown files for reference.

## Licence

MIT. The data are reused under the Open Government Licence v3.0.
