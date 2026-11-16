#' stillbirthmap: Bayesian small-area mapping of stillbirth risk in England
#'
#' Datasets: [stillbirths_lad_2024], [stillbirths_lad_2025], [lad_sf], [lad_nb],
#' [lad_utla_lookup], [obesity_pregnancy_utla] (load with `data()`).
#' Main functions: [fit_stillbirth_model()], [rate_ratios()], [exceedance()],
#' [map_exceedance()], [cross_validate()], [drop_one_cv()].
#' @keywords internal
#' @importFrom ggplot2 .data
#' @importFrom sf st_geometry
#' @importFrom spdep card
"_PACKAGE"
