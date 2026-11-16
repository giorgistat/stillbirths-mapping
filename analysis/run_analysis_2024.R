# =============================================================================
# analysis/run_analysis_2024.R -- full analysis for the paper (2024 births)
# Run from the repository root after installing the package:
#   R CMD INSTALL .   (or devtools::install())
#   Rscript analysis/run_analysis_2024.R
# Outputs go to output/2024/. Run time on 2 cores: ~45 minutes.
# =============================================================================
library(stillbirthmap)
library(ggplot2)
data(stillbirths_lad_2024); data(lad_nb); data(lad_sf); data(obesity_pregnancy_utla)
dat <- stillbirths_lad_2024
out <- file.path("output", "2024"); dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")

# hospital/provider booking-BMI data, if supplied (unit, pct_obese_preg, n_bookings, period)
extra_file <- file.path("data-raw", "sources", "hospital_booking_obesity.csv")
obesity_extra <- if (file.exists(extra_file)) read.csv(extra_file) else NULL

# ---- 1. main model: all pre-specified covariates, obesity in pregnancy joint, IID ------
main <- fit_stillbirth_model(dat, main_covariates(), obesity_extra = obesity_extra)
print(main)
rr_main <- rate_ratios(main)
s <- main$fit$summary(c("beta0", "beta", "beta_ob", "sigma", "a0", "gamma", "tau", "delta", "omega"),
                      "median", ~quantile(.x, c(0.025, 0.975)), "rhat", "ess_bulk")
dg <- main$fit$diagnostic_summary(quiet = TRUE)
write.csv(rr_main, file.path(out, "rate_ratios_main.csv"), row.names = FALSE)
write.csv(s, file.path(out, "parameters_main.csv"), row.names = FALSE)
write.csv(latent_obesity(main), file.path(out, "latent_obesity_pregnancy_LAD.csv"), row.names = FALSE)
log_msg("main model: divergences", sum(dg$num_divergent), "max Rhat", round(max(s$rhat), 3),
        "min ESS", round(min(s$ess_bulk)))

# ---- 2. exceedance and maps (main model) ------------------------------------------
exc <- exceedance(main)
write.csv(exc, file.path(out, "exceedance_main.csv"), row.names = FALSE)
ggsave(file.path(out, "figures", "fig_exceedance_main.png"), map_exceedance(exc), width = 12, height = 6.2, dpi = 300)
ggsave(file.path(out, "figures", "fig_rates_main.png"), map_rates(exc), width = 9, height = 5.6, dpi = 300)
thr <- attr(exc, "thresholds")
exc_tab <- do.call(rbind, lapply(names(thr), function(t) data.frame(
  threshold = t, rate = thr[[t]], n_high = sum(exc[[paste0("p_exceed_", t)]] > 0.8),
  n_low = sum(exc[[paste0("p_exceed_", t)]] < 0.2),
  high_units = paste(exc$name[exc[[paste0("p_exceed_", t)]] > 0.8 & thr[[t]] > attr(exc, "england_rate")],
                     collapse = "; "))))
write.csv(exc_tab, file.path(out, "exceedance_summary.csv"), row.names = FALSE)

# ---- 3. comparison and sensitivity fits -------------------------------------------
adult <- c("imd_dep_pctile", "pct_black", "pct_south_asian", "pct_age35plus", "pct_smoking", "pct_obese_adult")
comp <- list(
  null_iid        = fit_stillbirth_model(dat, character(0)),
  deprivation_iid = fit_stillbirth_model(dat, "imd_dep_pctile"),
  adult_obesity_iid    = fit_stillbirth_model(dat, adult),
  adult_obesity_bym2   = fit_stillbirth_model(dat, adult, "bym2", nb = lad_nb),
  adult_obesity_leroux = fit_stillbirth_model(dat, adult, "leroux", nb = lad_nb),
  null_bym2       = fit_stillbirth_model(dat, character(0), "bym2", nb = lad_nb),
  main_wider_prior = fit_stillbirth_model(dat, main_covariates(), sigma_scale = 2.5,
                                          obesity_extra = obesity_extra),
  main_excl_dq    = fit_stillbirth_model(dat[!(dat$age_dq_flag | dat$smoking_dq_flag), ], main_covariates(),
                                         obesity_extra = obesity_extra))
rr_comp <- do.call(rbind, lapply(names(comp), function(m) {
  r <- rate_ratios(comp[[m]]); if (nrow(r)) cbind(model = m, r) else NULL }))
write.csv(rr_comp, file.path(out, "rate_ratios_comparison_sensitivity.csv"), row.names = FALSE)
par_comp <- do.call(rbind, lapply(names(comp), function(m) {
  v <- intersect(c("sigma", "rho"), comp[[m]]$fit$metadata()$stan_variables)
  if (comp[[m]]$random_effect == "iid") v <- "sigma"
  x <- comp[[m]]$fit$summary(v, "median", ~quantile(.x, c(0.025, 0.975)), "rhat")
  cbind(model = m, x) }))
write.csv(par_comp, file.path(out, "random_effect_parameters.csv"), row.names = FALSE)
log_msg("comparison fits done")

# overdispersion and spatial autocorrelation diagnostics
E <- dat$total_births * sum(dat$stillbirths) / sum(dat$total_births)
X2 <- sum((dat$stillbirths - E)^2 / E)
lw <- spdep::nb2listw(lad_nb, style = "W")
set.seed(1)
eb <- spdep::EBImoran.mc(dat$stillbirths, dat$total_births, lw, nsim = 9999)
g1 <- glm(stillbirths ~ imd_dep_pctile + offset(log(total_births)), family = poisson, data = dat)
mr <- spdep::moran.mc(residuals(g1, type = "pearson"), lw, nsim = 9999)
diag_tab <- data.frame(
  statistic = c("Pearson X2 (constant rate)", "df", "dispersion", "p",
                "EB Moran's I (rates)", "p", "Moran's I (Pearson residuals | deprivation)", "p"),
  value = c(X2, nrow(dat) - 1, X2 / (nrow(dat) - 1), pchisq(X2, nrow(dat) - 1, lower.tail = FALSE),
            eb$statistic, eb$p.value, mr$statistic, mr$p.value))
write.csv(diag_tab, file.path(out, "dispersion_moran.csv"), row.names = FALSE)

# ---- 4. cross-validation (descriptive: which covariates are most predictive) -----
folds <- make_folds(dat)
cv <- list(
  main            = cross_validate(dat, main_covariates(), folds, obesity_extra = obesity_extra),
  null_iid        = cross_validate(dat, character(0), folds),
  null_bym2       = cross_validate(dat, character(0), folds, random_effect = "bym2", nb = lad_nb),
  deprivation_iid = cross_validate(dat, "imd_dep_pctile", folds),
  adult_obesity_iid = cross_validate(dat, adult, folds))
log_msg("CV comparison done")
cv_summ <- summarise_cv(cv, reference = "null_iid")
write.csv(cv_summ, file.path(out, "cv_summary.csv"), row.names = FALSE)

d1 <- drop_one_cv(dat, main_covariates(), folds, cv_full = cv$main, obesity_extra = obesity_extra)
write.csv(d1$summary, file.path(out, "cv_drop_one.csv"), row.names = FALSE)
saveRDS(list(cv = cv, drop_one = d1$pointwise), file.path(out, "cv_pointwise.rds"))
log_msg("drop-one CV done")
print(cv_summ); print(d1$summary)

# ---- 5. figures: rate ratios and covariate contributions --------------------------
rr_plot <- rbind(cbind(model = "Main (obesity in pregnancy, joint)", rr_main),
                 cbind(model = "Adult obesity instead", rate_ratios(comp$adult_obesity_iid)))
rr_plot$lab <- paste0(rr_plot$label, "\n(per ", rr_plot$per, ")")
rr_plot$lab <- factor(rr_plot$lab, levels = rev(unique(rr_plot$lab)))
p <- ggplot(rr_plot, aes(RR, lab, colour = model)) +
  geom_vline(xintercept = 1, colour = "grey55", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.8, position = position_dodge(0.5)) +
  geom_point(size = 2, position = position_dodge(0.5)) +
  scale_colour_manual(values = c("Main (obesity in pregnancy, joint)" = "#1c5cab",
                                 "Adult obesity instead" = "#86b6ef"), name = NULL) +
  scale_x_log10() + labs(x = "Rate ratio (95% credible interval), log scale", y = NULL) +
  theme_minimal(base_size = 9) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())
ggsave(file.path(out, "figures", "fig_rate_ratios.png"), p, width = 6.5, height = 4.5, dpi = 300)

dd <- d1$summary; dd$label <- factor(dd$label, levels = rev(dd$label))
p2 <- ggplot(dd, aes(dELPD, label)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = dELPD - 2 * se, xmax = dELPD + 2 * se), height = 0, colour = "#86b6ef", linewidth = 0.9) +
  geom_point(colour = "#1c5cab", size = 2) +
  labs(x = "Loss in 10-fold CV log predictive density when the covariate is dropped (±2 SE)", y = NULL) +
  theme_minimal(base_size = 9) + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())
ggsave(file.path(out, "figures", "fig_cv_drop_one.png"), p2, width = 7, height = 3.2, dpi = 300)

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
log_msg("all done")
