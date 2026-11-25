# stillbirthmap

R package with the data and code to reproduce the analysis in

> Giorgi E, Basra R, Heazell AEP. *Mapping stillbirth risk across local authorities in England
> by integrating multiple national data sources.* Manuscript, 2026 (*Paediatric and Perinatal
> Epidemiology*).

Stillbirth rates in England remain above the national ambition and vary substantially between
local areas. Because stillbirth is a relatively rare outcome, local rates are difficult to
estimate precisely, and the national data on maternal risk factors are released at different
geographies and are often incomplete. This analysis combines birth registration, Census,
deprivation and maternity data in a Poisson generalised linear mixed model to estimate stillbirth
risk for all 294 local authority districts (LADs) of England in 2024. Exceedance probabilities
then identify the districts with a high probability of rates above the England rate and other
policy-relevant thresholds. The work builds on Ranveer Basra's MSc thesis (University of
Birmingham, 2026).

## Main findings

- District characteristics explained most of the variation in stillbirth risk between districts
  (an 82.1% reduction in the variance of the district random effects). The South Asian population
  share was the most predictive indicator and the most clearly associated with risk (rate ratio
  1.047 per 5 percentage points, 95% credible interval 1.012, 1.084).
- Thirty-three districts had a probability above 0.8 of exceeding the England rate (3.83 per
  1,000 total births). Seven had a probability above 0.8 of exceeding it by 20%: Leicester,
  Sandwell, Birmingham, Barking and Dagenham, Luton, Newham and Wolverhampton.
- In 291 of the 294 districts, the probability of exceeding the national ambition of 2.6 per
  1,000 was above 0.8, so the gap to the target is not confined to a few high-risk areas.

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
| Black and South Asian population shares (%) | Census 2021, table TS022 | LAD |
| % deliveries to mothers aged 35+ | OHID Fingertips, 2024/25 | LAD |
| % smoking at time of delivery | OHID Fingertips, 2024/25 | LAD |
| % adult obesity (obesity sub-model and sensitivity analysis) | OHID Fingertips (Active Lives Survey), 2024/25 | LAD |
| % women with BMI 30+ at the booking appointment | OHID Fingertips 94131 | Upper-tier authority, about a third withheld |

Fingertips flags 50 LADs for data-quality concerns in the maternal age or smoking indicators. They
are kept in the main analysis and excluded in a sensitivity analysis. For Torbay, maternal age was
not published for 2022/23 to 2024/25, so its most recent published value (2021/22) is used.

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

All six covariates are pre-specified and kept in the model. Ten-fold cross-validation, scored
with the logarithmic and ranked probability scores, is used to describe each covariate's
predictive contribution, not to select covariates. The calibration of the cross-validated
predictions is checked against replicate counts simulated from the model, and Moran's I of the
residuals is used to check for residual spatial correlation.

For each LAD, the exceedance probability is the posterior probability that its stillbirth rate is
above a threshold. Thresholds are the England rate (3.83 per 1,000 total births), 10% and 20%
above it (4.21 and 4.59), and the national ambition (2.6). LADs with a probability above 0.8 are
classed as high-probability areas and those below 0.2 as low-probability areas.

BYM2 and Leroux spatial random effects are available through `random_effect =`. They did not
improve the fit, which suggests that any residual spatial variation occurs at a scale finer than
the LAD.

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
Rscript analysis/cv_calibration_2024.R   # calibration checks, Figures S3 and S4
```

Results are written to `output/2024/`: rate ratios, model parameters, exceedance probabilities
for every LAD, latent obesity estimates, cross-validation summaries, residual checks, sensitivity
analyses and figures. The Supplementary Material covers data harmonisation (S1), the main model
(S2), the joint obesity model (S3), residual spatial correlation (S4), cross-validation and
calibration (S5) and sensitivity analyses (S6).

`legacy/ranveer_msc/` holds the original thesis R Markdown files for reference.

## Licence

MIT. The data are reused under the Open Government Licence v3.0.
