#' LAD-level stillbirths and covariates, England
#'
#' One row per modelling unit (294 units: English local authority districts,
#' May 2024 geography, with Cornwall + Isles of Scilly and Hackney + City of
#' London merged). Use `data(stillbirths_lad_2024)` (main analysis) or
#' `data(stillbirths_lad_2025)`.
#'
#' @format A data frame with 294 rows and columns:
#' \describe{
#'   \item{unit}{LAD24 code of the modelling unit}
#'   \item{name}{Unit name}
#'   \item{year}{Registration year of births and stillbirths}
#'   \item{live_births, stillbirths, total_births}{Counts by mother's usual residence (ONS
#'     birth registrations, Table 2). Stillbirths are >= 24 weeks; total = live + still.}
#'   \item{imd_dep_pctile}{Deprivation: mean over LSOAs of 100 x (1 - rank/33,755), IoD 2025
#'     (higher = more deprived)}
#'   \item{pct_lsoa_most_dep_quintile}{\% of LSOAs in the most deprived 20\% nationally}
#'   \item{pct_black, pct_south_asian, pct_non_white_british}{Census 2021, whole population}
#'   \item{pct_age35plus}{\% of deliveries to women aged 35+ (Fingertips 92250)}
#'   \item{age_dq_flag, age_period}{Fingertips data-quality flag; period used (Torbay is
#'     suppressed from 2022/23 and uses 2021/22)}
#'   \item{pct_smoking, smoking_dq_flag}{Smoking at time of delivery (Fingertips 93085;
#'     resolved at NHS sub-ICB level)}
#'   \item{pct_obese_adult}{Adult obesity prevalence (Fingertips 93881, all adults)}
#'   \item{utla}{Upper-tier local authority of the unit}
#'   \item{fingertips_period}{Fingertips period used for the maternity indicators}
#' }
#' @source ONS Births in England and Wales: birth registrations (2024, 2025 editions);
#'   MHCLG English Indices of Deprivation 2025; ONS Census 2021 (TS ethnic group, detailed);
#'   OHID Fingertips. Open Government Licence v3.0.
"stillbirths_lad_2024"

#' @rdname stillbirths_lad_2024
"stillbirths_lad_2025"

#' Boundaries of the 294 modelling units
#'
#' `sf` polygons (British National Grid, EPSG:27700), derived from ONS LAD May 2024
#' BGC boundaries, merged as in [stillbirths_lad_2024] and lightly simplified for mapping.
#' Rows are in the same order as the stillbirth datasets.
#' @format An `sf` data frame with columns `unit`, `name`, `geometry`.
#' @source ONS Open Geography Portal, OGL v3.0.
"lad_sf"

#' Adjacency of the modelling units
#'
#' `spdep` `nb` object (queen contiguity on the full-resolution boundaries, 50 m snap).
#' The Isle of Wight is linked to Portsmouth, Southampton and New Forest (ferry routes).
#' Same order as [lad_sf].
#' @format An object of class `nb`.
"lad_nb"

#' Modelling unit to upper-tier local authority lookup
#' @format A data frame with columns `unit`, `utla`, `utla_name`.
#' @source OHID Fingertips area hierarchy.
"lad_utla_lookup"

#' Obesity in early pregnancy, upper-tier local authorities
#'
#' Published values of Fingertips indicator 94131 (\% of women with BMI >= 30 at booking).
#' Only upper-tier authorities are published, and many are suppressed.
#' @format A data frame with columns `utla`, `utla_name`, `period` ("2022/23", "2023/24"),
#'   `pct_obese_preg`, `n_obese`, `n_bookings`.
#' @source OHID Fingertips, OGL v3.0.
"obesity_pregnancy_utla"
