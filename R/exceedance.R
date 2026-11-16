#' Exceedance probabilities for policy thresholds
#'
#' @param object A `stillbirth_fit`.
#' @param thresholds Named numeric vector of rates per 1,000 total births. Defaults to the
#'   England rate of the data year, England +10\% and +20\%, and the 2.6 national ambition.
#' @return Data frame per unit with observed and posterior rates (median, 95\% interval),
#'   exceedance probabilities `p_exceed_<name>` and classes (P > 0.8 high, P < 0.2 low).
#' @export
exceedance <- function(object, thresholds = NULL) {
  d <- object$data; R <- rate_draws(object)
  eng <- 1000 * sum(d$stillbirths) / sum(d$total_births)
  if (is.null(thresholds))
    thresholds <- c(ambition = 2.6, england = round(eng, 2), plus10 = round(1.1 * eng, 2),
                    plus20 = round(1.2 * eng, 2))
  res <- data.frame(unit = d$unit, name = d$name, total_births = d$total_births,
                    stillbirths = d$stillbirths, observed_rate = 1000 * d$stillbirths / d$total_births,
                    post_median = unname(apply(R, 2, stats::median)),
                    post_lo = unname(apply(R, 2, stats::quantile, 0.025)),
                    post_hi = unname(apply(R, 2, stats::quantile, 0.975)))
  for (t in names(thresholds)) {
    p <- unname(colMeans(R > thresholds[[t]]))
    res[[paste0("p_exceed_", t)]] <- p
    res[[paste0("class_", t)]] <- cut(p, c(-Inf, 0.2, 0.8, Inf),
                                      labels = c("Low (P < 0.2)", "Uncertain", "High (P > 0.8)"))
  }
  attr(res, "thresholds") <- thresholds
  attr(res, "england_rate") <- eng
  attr(res, "year") <- d$year[1]
  res
}

#' Map exceedance probabilities
#'
#' @param exc Output of [exceedance()].
#' @param which Threshold names to map.
#' @param sf Boundaries, default [lad_sf].
#' @param london_inset Add a London inset.
#' @return A patchwork/ggplot object.
#' @export
map_exceedance <- function(exc, which = c("england", "plus10", "plus20"), sf = NULL,
                           london_inset = TRUE) {
  if (is.null(sf)) sf <- get_data("lad_sf")
  thr <- attr(exc, "thresholds")
  m <- join_sf(sf, exc)
  brks <- c(0, 0.2, 0.5, 0.8, 0.95, 1)
  labs <- c("< 0.20", "0.20-0.50", "0.50-0.80", "0.80-0.95", "> 0.95")
  cols <- stats::setNames(c("#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"), labs)
  ttl <- c(ambition = "national ambition", england = paste("England", attr(exc, "year")),
           plus10 = "England +10%", plus20 = "England +20%")
  panels <- lapply(seq_along(which), function(k) {
    w <- which[k]
    mm <- m; mm$cls <- cut(mm[[paste0("p_exceed_", w)]], brks, labs, include.lowest = TRUE)
    mk <- function(x) ggplot2::ggplot(x) +
      ggplot2::geom_sf(ggplot2::aes(fill = .data$cls), colour = "white", linewidth = 0.08,
                       show.legend = TRUE) +
      ggplot2::scale_fill_manual(values = cols, drop = FALSE, name = "Posterior exceedance probability") +
      ggplot2::theme_void(base_size = 9) +
      ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1, title.position = "top"))
    p <- mk(mm) + ggplot2::ggtitle(sprintf("%s. P(rate > %.2f%s)", LETTERS[k], thr[[w]],
                                           if (w %in% names(ttl)) paste0(", ", ttl[[w]]) else "")) +
      ggplot2::theme(legend.position = if (k == ceiling(length(which) / 2)) "bottom" else "none",
                     plot.title = ggplot2::element_text(size = 9, face = "bold"))
    if (london_inset) {
      lon <- mm[substr(mm$unit, 1, 3) == "E09", ]
      p <- p + patchwork::inset_element(mk(lon) + ggplot2::theme(legend.position = "none",
                 panel.border = ggplot2::element_rect(fill = NA, colour = "grey60", linewidth = 0.3)),
               left = 0, bottom = 0.55, right = 0.32, top = 0.85, align_to = "panel")
    }
    p
  })
  patchwork::wrap_plots(panels, nrow = 1)
}

#' Map observed and model-smoothed rates
#' @inheritParams map_exceedance
#' @return A patchwork object (observed | posterior median), diverging scale centred on England.
#' @export
map_rates <- function(exc, sf = NULL, london_inset = TRUE) {
  if (is.null(sf)) sf <- get_data("lad_sf")
  eng <- attr(exc, "england_rate"); m <- join_sf(sf, exc)
  one <- function(var, title) {
    mk <- function(x) ggplot2::ggplot(x) +
      ggplot2::geom_sf(ggplot2::aes(fill = log(pmax(.data[[var]], 1.5))), colour = "white", linewidth = 0.08) +
      ggplot2::scale_fill_gradient2(low = "#1c5cab", mid = "#f0efec", high = "#c0302f",
        midpoint = log(eng), limits = log(c(1.5, 10)), oob = scales::squish,
        breaks = log(c(2, 3, eng, 5, 7, 10)), labels = c("2", "3", sprintf("%.1f", eng), "5", "7", "10"),
        name = "Stillbirths per\n1,000 births") + ggplot2::theme_void(base_size = 9)
    p <- mk(m) + ggplot2::ggtitle(title) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 9, face = "bold"))
    if (london_inset) p <- p + patchwork::inset_element(
      mk(m[substr(m$unit, 1, 3) == "E09", ]) + ggplot2::theme(legend.position = "none",
        panel.border = ggplot2::element_rect(fill = NA, colour = "grey60", linewidth = 0.3)),
      left = 0, bottom = 0.55, right = 0.32, top = 0.85, align_to = "panel")
    p
  }
  patchwork::wrap_plots(one("observed_rate", "A. Observed rate"),
                        one("post_median", "B. Model-smoothed rate (posterior median)"), nrow = 1)
}

join_sf <- function(sf, df) {
  requireNamespace("sf", quietly = TRUE)
  idx <- match(sf$unit, df$unit)
  for (v in setdiff(names(df), c("unit", "name"))) sf[[v]] <- df[[v]][idx]
  sf
}
