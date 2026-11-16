#' Covariate labels and reporting units
#'
#' @return A data frame with the covariate name (`var`), a label, and the increment
#'   (`per`, with text `unit`) used to report rate ratios.
#' @export
covariate_info <- function() {
  data.frame(
    var   = c("imd_dep_pctile", "pct_lsoa_most_dep_quintile", "pct_black", "pct_south_asian",
              "pct_non_white_british", "pct_age35plus", "pct_smoking", "pct_obese_adult",
              "pct_obese_preg"),
    label = c("Deprivation (IMD 2025, mean LSOA percentile)",
              "LSOAs in most deprived 20% (%)",
              "Black ethnicity (% population)", "South Asian ethnicity (% population)",
              "Non-White-British (% population)", "Mothers aged 35+ (% deliveries)",
              "Smoking at time of delivery (%)", "Adult obesity (%)",
              "Obesity in early pregnancy (%)"),
    per   = c(10, 10, 5, 5, 10, 5, 2, 5, 5),
    unit  = c("10 percentile points", "10 percentage points", "5 percentage points",
              "5 percentage points", "10 percentage points", "5 percentage points",
              "2 percentage points", "5 percentage points", "5 percentage points"),
    stringsAsFactors = FALSE)
}

#' Main covariate set of the paper
#' @return Character vector: deprivation, \% Black, \% South Asian, mothers 35+, smoking and
#'   obesity in early pregnancy (`"pct_obese_preg"`, modelled jointly).
#' @export
main_covariates <- function() {
  c("imd_dep_pctile", "pct_black", "pct_south_asian", "pct_age35plus", "pct_smoking",
    "pct_obese_preg")
}

# Compile (once per session/cache) one of the package's Stan models
stan_model_of <- function(name) {
  cache <- tools::R_user_dir("stillbirthmap", "cache")
  dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  f <- system.file("stan", paste0(name, ".stan"), package = "stillbirthmap", mustWork = TRUE)
  local <- file.path(cache, basename(f))
  if (!file.exists(local) || !identical(readLines(local), readLines(f))) file.copy(f, local, overwrite = TRUE)
  cmdstanr::cmdstan_model(local, quiet = TRUE)
}

#' Graph quantities for areal models
#'
#' @param nb An `spdep` `nb` object (single connected component for BYM2).
#' @return A list with the edge list (`node1`, `node2`), the BYM2 scaling factor
#'   (Riebler et al. 2016) and the eigen-decomposition of `D - W` (used by the Leroux model).
#' @export
lad_graph <- function(nb) {
  N <- length(nb)
  e <- do.call(rbind, lapply(seq_len(N), function(i) {
    j <- nb[[i]]; j <- j[j > i]; if (length(j)) cbind(i, j) else NULL }))
  W <- matrix(0, N, N); W[e] <- 1; W <- W + t(W)
  R <- diag(rowSums(W)) - W
  Sig <- solve(R + diag(max(diag(R)) * sqrt(.Machine$double.eps), N))
  A <- matrix(1, 1, N)
  Sig_c <- Sig - Sig %*% t(A) %*% solve(A %*% Sig %*% t(A)) %*% A %*% Sig
  list(node1 = e[, 1], node2 = e[, 2], N_edges = nrow(e),
       scaling_factor = exp(mean(log(diag(Sig_c)))), eig = eigen(R, symmetric = TRUE))
}
