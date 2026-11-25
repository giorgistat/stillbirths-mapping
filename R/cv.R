#' Cross-validation folds stratified by births
#' @param data A stillbirth dataset.
#' @param K Number of folds.
#' @param seed Random seed.
#' @return Integer vector of fold labels (one per unit).
#' @export
make_folds <- function(data, K = 10, seed = 20251116) {
  set.seed(seed)
  ord <- order(data$total_births); fold <- integer(nrow(data))
  for (blk in split(ord, ceiling(seq_along(ord) / K))) fold[blk] <- sample(K, length(blk))
  fold
}

#' Proper scoring rules for held-out counts
#'
#' @param y Observed counts.
#' @param log_offset log(total births).
#' @param log_rate_draws Matrix (draws x units) of posterior log rates.
#' @return Data frame with log predictive density (`lpd`), ranked probability score (`rps`),
#'   predictive mean, squared and absolute error, and 90\% interval coverage.
#' @export
score_counts <- function(y, log_offset, log_rate_draws) {
  mu <- exp(sweep(log_rate_draws, 2, log_offset, "+")); S <- nrow(mu)
  lpd <- vapply(seq_along(y), function(i) {
    lp <- stats::dpois(y[i], mu[, i], log = TRUE); m <- max(lp); m + log(mean(exp(lp - m))) }, 0)
  rps <- vapply(seq_along(y), function(i) {
    k <- 0:(max(y[i], stats::qpois(0.9999, max(mu[, i]))) + 1)
    Fk <- colMeans(outer(mu[, i], k, function(m, kk) stats::ppois(kk, m)))
    sum((Fk - as.numeric(k >= y[i]))^2) }, 0)
  pm <- colMeans(mu)
  cover <- vapply(seq_along(y), function(i) {
    q <- stats::quantile(stats::rpois(S, mu[, i]), c(0.05, 0.95)); y[i] >= q[1] & y[i] <= q[2] }, TRUE)
  data.frame(lpd = lpd, rps = rps, pred_mean = pm, sq_err = (y - pm)^2, abs_err = abs(y - pm),
             cover90 = cover)
}

#' K-fold cross-validation of a model specification
#'
#' Each fold's units keep their covariates (and, in the joint model, their obesity data)
#' but their stillbirth counts are removed from the likelihood and predicted.
#' @inheritParams fit_stillbirth_model
#' @param folds Fold labels from [make_folds()].
#' @param iter_warmup,iter_sampling MCMC iterations per fold.
#' @param ... Further arguments to [fit_stillbirth_model()].
#' @param keep_draws If `TRUE`, the posterior draws of the held-out rates (per 1,000 total
#'   births; draws x units, columns in the order of `data`) are attached as attribute
#'   `"rate_draws"`, so that replicate counts can be simulated from the predictive distribution.
#' @return Pointwise scores with unit, fold and observed count.
#' @export
cross_validate <- function(data, covariates, folds = make_folds(data), random_effect = "iid",
                           iter_warmup = 500, iter_sampling = 500, keep_draws = FALSE, ...) {
  draws <- NULL
  out <- lapply(sort(unique(folds)), function(k) {
    test <- which(folds == k)
    f <- fit_stillbirth_model(data, covariates, random_effect = random_effect, holdout = test,
                              iter_warmup = iter_warmup, iter_sampling = iter_sampling,
                              seed = 100 + k, ...)
    rd <- rate_draws(f)[, test, drop = FALSE]
    if (keep_draws) draws[[k]] <<- rd
    lr <- log(rd / 1000)
    cbind(data.frame(unit = data$unit[test], fold = k, y = data$stillbirths[test],
                     births = data$total_births[test]),
          score_counts(data$stillbirths[test], log(data$total_births[test]), lr))
  })
  res <- do.call(rbind, out)
  if (keep_draws) {
    m <- matrix(NA_real_, nrow(draws[[1]]), nrow(data), dimnames = list(NULL, data$unit))
    for (k in seq_along(draws)) m[, folds == sort(unique(folds))[k]] <- draws[[k]]
    attr(res, "rate_draws") <- m
  }
  res
}

#' Summarise cross-validation results against a reference model
#' @param cv_list Named list of outputs of [cross_validate()].
#' @param reference Name of the reference element (e.g. the null model).
#' @return Data frame with ELPD, difference from the reference with paired SE, RPS change,
#'   RMSE, MAE and 90\% coverage.
#' @export
summarise_cv <- function(cv_list, reference) {
  ref <- cv_list[[reference]]
  do.call(rbind, lapply(names(cv_list), function(m) {
    x <- cv_list[[m]]; r <- ref[match(x$unit, ref$unit), ]
    d <- x$lpd - r$lpd
    data.frame(model = m, ELPD = sum(x$lpd), dELPD = sum(d), se_dELPD = sqrt(length(d) * stats::var(d)),
               dRPS_pct = 100 * (mean(x$rps) / mean(r$rps) - 1), RMSE = sqrt(mean(x$sq_err)),
               MAE = mean(x$abs_err), coverage90 = mean(x$cover90))
  }))
}

#' Predictive contribution of each covariate (drop-one cross-validation)
#'
#' Refits the full model leaving out one covariate at a time. The contribution of a covariate is
#' ELPD(full) - ELPD(full without it), with a paired standard error. Used descriptively, to rank
#' covariates by predictive value, not to select them.
#' @inheritParams cross_validate
#' @param cv_full Output of [cross_validate()] for the full model (computed if `NULL`).
#' @return List with `summary` (one row per covariate) and `pointwise` results.
#' @export
drop_one_cv <- function(data, covariates = main_covariates(), folds = make_folds(data),
                        cv_full = NULL, ...) {
  if (is.null(cv_full)) cv_full <- cross_validate(data, covariates, folds, ...)
  pw <- lapply(covariates, function(v) cross_validate(data, setdiff(covariates, v), folds, ...))
  names(pw) <- covariates
  info <- covariate_info()
  summ <- do.call(rbind, lapply(covariates, function(v) {
    x <- pw[[v]]; f <- cv_full[match(x$unit, cv_full$unit), ]
    d <- f$lpd - x$lpd
    data.frame(dropped = v, label = info$label[match(v, info$var)], dELPD = sum(d),
               se = sqrt(length(d) * stats::var(d)),
               dRPS_pct = 100 * (mean(x$rps) / mean(f$rps) - 1))
  }))
  list(summary = summ[order(-summ$dELPD), ], pointwise = pw, full = cv_full)
}
