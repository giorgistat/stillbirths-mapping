#' Fit the stillbirth small-area model
#'
#' Poisson log-linear model for LAD stillbirth counts with offset log(total births)
#' and LAD random effects:
#' \deqn{Y_i \sim \mathrm{Poisson}(B_i\lambda_i),\quad \log\lambda_i = \beta_0 + x_i^\top\beta + b_i.}
#' If `"pct_obese_preg"` is among the covariates, obesity in early pregnancy is treated as a
#' latent LAD-level variable, \eqn{o_i = a_0 + z_i^\top\gamma + u_{g(i)}}, observed only as
#' births-weighted averages over upper-tier authorities (`obesity_data`) with binomial
#' sampling error, and estimated jointly with the stillbirth model.
#'
#' @param data A stillbirth dataset, e.g. [stillbirths_lad_2024].
#' @param covariates Character vector of covariate columns; include `"pct_obese_preg"`
#'   for the joint obesity model. Default [main_covariates()].
#' @param random_effect `"iid"` (default), `"bym2"` or `"leroux"`. The joint obesity model
#'   supports `"iid"` only.
#' @param nb Adjacency (`nb`) for spatial random effects; default [lad_nb].
#' @param obesity_data Upper-tier observations, default [obesity_pregnancy_utla].
#' @param obesity_extra Optional LAD/provider-level observations (e.g. hospital booking
#'   data): a data frame with `unit`, `pct_obese_preg`, `n_bookings`, `period`.
#' @param obesity_predictors LAD covariates of the obesity sub-model.
#' @param sigma_scale Scale of the half-normal prior on the random-effect SD.
#' @param holdout Optional logical/integer vector of units whose counts are held out of the
#'   likelihood (used for cross-validation).
#' @param chains,iter_warmup,iter_sampling,seed,parallel_chains Passed to cmdstanr.
#' @return An object of class `stillbirth_fit`.
#' @export
fit_stillbirth_model <- function(data, covariates = main_covariates(),
                                 random_effect = c("iid", "bym2", "leroux"),
                                 nb = NULL, obesity_data = NULL, obesity_extra = NULL,
                                 obesity_predictors = c("pct_obese_adult", "pct_age35plus",
                                                        "imd_dep_pctile", "pct_black",
                                                        "pct_south_asian", "pct_smoking"),
                                 sigma_scale = 1, holdout = NULL, chains = 4,
                                 iter_warmup = 1000, iter_sampling = 1000, seed = 2026,
                                 parallel_chains = min(chains, max(1, parallel::detectCores() - 0))) {
  random_effect <- match.arg(random_effect)
  N <- nrow(data)
  in_lik <- rep(1L, N)
  if (!is.null(holdout)) {
    if (is.logical(holdout)) holdout <- which(holdout)
    in_lik[holdout] <- 0L
  }
  joint <- "pct_obese_preg" %in% covariates
  xvars <- setdiff(covariates, "pct_obese_preg")
  X <- if (length(xvars)) scale(as.matrix(data[, xvars, drop = FALSE])) else matrix(0, N, 0)
  base <- list(N = N, y = as.integer(data$stillbirths), log_offset = log(data$total_births),
               K = ncol(X), X = matrix(X, N), in_lik = in_lik, sigma_scale = sigma_scale)
  obs_info <- NULL
  if (!joint) {
    if (is.null(nb)) nb <- get_data("lad_nb")
    re_type <- c(iid = 0L, bym2 = 1L, leroux = 2L)[[random_effect]]
    g <- lad_graph(nb)
    sd <- c(base, list(N_edges = g$N_edges, node1 = g$node1, node2 = g$node2, re_type = re_type,
                       scaling_factor = g$scaling_factor, lambda_W = pmax(g$eig$values, 0),
                       V_W = if (re_type == 2) g$eig$vectors else matrix(0, 0, 0)))
    mod <- stan_model_of("areal_poisson")
  } else {
    if (random_effect != "iid") stop("The joint obesity model supports random_effect = 'iid' only.")
    if (is.null(obesity_data)) obesity_data <- get_data("obesity_pregnancy_utla")
    obs_info <- obesity_observations(data, obesity_data, obesity_extra)
    Z <- scale(as.matrix(data[, obesity_predictors, drop = FALSE]))
    sd <- c(base, list(G = obs_info$G, grp = obs_info$grp, P = ncol(Z), Z = unclass(Z)[, , drop = FALSE],
                       M = length(obs_info$z), Wobs = obs_info$W, zobs = obs_info$z,
                       se_obs = obs_info$se, is_2223 = obs_info$is2223,
                       ob_c = mean(obs_info$z), ob_s = stats::sd(obs_info$z)))
    mod <- stan_model_of("joint_obesity")
  }
  init_f <- function() {
    l <- list(beta0 = log(sum(data$stillbirths) / sum(data$total_births)),
              beta = as.array(rep(0, sd$K)), sigma = 0.1, theta = stats::rnorm(N, 0, 0.1))
    if (!joint) { l$rho <- 0.5; l$phi <- if (sd$re_type == 1) stats::rnorm(N, 0, 0.1) else numeric(0) }
    l
  }
  fit <- mod$sample(data = sd, chains = chains, parallel_chains = parallel_chains,
                    iter_warmup = iter_warmup, iter_sampling = iter_sampling, seed = seed,
                    init = init_f, adapt_delta = 0.95, max_treedepth = 12, refresh = 0,
                    show_messages = FALSE, show_exceptions = FALSE)
  structure(list(fit = fit, data = data, covariates = covariates, xvars = xvars,
                 x_scale = if (length(xvars)) attr(X, "scaled:scale") else numeric(0),
                 ob_scale = if (joint) stats::sd(obs_info$z) else NA_real_,
                 joint = joint, random_effect = random_effect, in_lik = in_lik,
                 obesity = obs_info),
            class = "stillbirth_fit")
}

get_data <- function(name) {
  e <- new.env()
  utils::data(list = name, package = "stillbirthmap", envir = e)
  e[[name]]
}

obesity_observations <- function(data, obesity_data, obesity_extra = NULL) {
  utlas <- sort(unique(data$utla)); grp <- match(data$utla, utlas); N <- nrow(data)
  ob <- obesity_data[obesity_data$utla %in% utlas & !is.na(obesity_data$pct_obese_preg), ]
  W <- t(vapply(seq_len(nrow(ob)), function(k) {
    w <- numeric(N); m <- which(data$utla == ob$utla[k])
    w[m] <- data$total_births[m] / sum(data$total_births[m]); w }, numeric(N)))
  p <- ob$pct_obese_preg / 100
  z <- ob$pct_obese_preg; se <- 100 * sqrt(p * (1 - p) / ob$n_bookings)
  is2223 <- as.numeric(ob$period == "2022/23")
  if (!is.null(obesity_extra) && nrow(obesity_extra)) {
    for (k in seq_len(nrow(obesity_extra))) {
      w <- numeric(N); w[match(obesity_extra$unit[k], data$unit)] <- 1
      W <- rbind(W, w); pk <- obesity_extra$pct_obese_preg[k] / 100
      z <- c(z, 100 * pk); se <- c(se, 100 * sqrt(pk * (1 - pk) / obesity_extra$n_bookings[k]))
      is2223 <- c(is2223, as.numeric(obesity_extra$period[k] == "2022/23"))
    }
  }
  list(G = length(utlas), grp = grp, W = W, z = z, se = se, is2223 = is2223,
       n_utla_observed = length(unique(ob$utla)), units_covered = sum(colSums(W) > 0))
}

#' @export
print.stillbirth_fit <- function(x, ...) {
  cat("Stillbirth small-area model,", x$data$year[1], "births\n")
  cat("Random effect:", x$random_effect, if (x$joint) "| obesity in pregnancy modelled jointly", "\n")
  cat("Covariates:", paste(x$covariates, collapse = ", "), "\n")
  dg <- x$fit$diagnostic_summary(quiet = TRUE)
  cat("Divergences:", sum(dg$num_divergent), "\n")
  print(rate_ratios(x)[, c("label", "per", "RR", "lo", "hi")], digits = 3, row.names = FALSE)
  invisible(x)
}

#' Rate ratios per reporting unit
#'
#' @param object A `stillbirth_fit`.
#' @param probs Credible interval probabilities.
#' @return Data frame with posterior median rate ratio, interval and P(RR > 1) per covariate,
#'   per the increments in [covariate_info()].
#' @export
rate_ratios <- function(object, probs = c(0.025, 0.975)) {
  info <- covariate_info()
  vars <- c(object$xvars, if (object$joint) "pct_obese_preg")
  if (!length(vars)) return(data.frame())
  B <- posterior::as_draws_matrix(object$fit$draws(c(if (length(object$xvars)) "beta",
                                                      if (object$joint) "beta_ob")))
  scl <- c(object$x_scale, if (object$joint) object$ob_scale)
  do.call(rbind, lapply(seq_along(vars), function(k) {
    i <- match(vars[k], info$var)
    x <- exp(as.numeric(B[, k]) * info$per[i] / scl[k])
    data.frame(var = vars[k], label = info$label[i], per = info$unit[i], RR = stats::median(x),
               lo = stats::quantile(x, probs[1]), hi = stats::quantile(x, probs[2]),
               p_gt1 = mean(x > 1), row.names = NULL)
  }))
}

#' Posterior draws of LAD stillbirth rates
#' @param object A `stillbirth_fit`.
#' @return Matrix (draws x units) of rates per 1,000 total births.
#' @export
rate_draws <- function(object) {
  1000 * exp(unclass(posterior::as_draws_matrix(object$fit$draws("log_rate"))))
}

#' Latent obesity in early pregnancy by LAD (joint model)
#' @param object A joint `stillbirth_fit`.
#' @return Data frame of posterior median and 95\% interval per unit.
#' @export
latent_obesity <- function(object) {
  stopifnot(object$joint)
  O <- unclass(posterior::as_draws_matrix(object$fit$draws("ob")))
  data.frame(unit = object$data$unit, name = object$data$name, utla = object$data$utla,
             median = unname(apply(O, 2, stats::median)), lo = unname(apply(O, 2, stats::quantile, 0.025)),
             hi = unname(apply(O, 2, stats::quantile, 0.975)))
}
