## Cross-validated calibration checks (Supplementary Section S5, Figures S3 and S4).
## Refits the ten-fold cross-validation keeping the held-out posterior draws of the rates, then
## simulates replicate counts from Poisson(births x rate) for each draw, so that the predictive
## intervals include both parameter uncertainty and Poisson sampling variation.
## Run from the package root after analysis/run_analysis_2024.R:  Rscript analysis/cv_calibration_2024.R
suppressPackageStartupMessages({ library(stillbirthmap); library(ggplot2); library(dplyr); library(patchwork) })
out <- "output/2024"; dir.create(out, showWarnings = FALSE, recursive = TRUE)

data(stillbirths_lad_2024); dat <- stillbirths_lad_2024
folds <- make_folds(dat)
cvd <- cross_validate(dat, main_covariates(), folds, keep_draws = TRUE)
saveRDS(list(scores = cvd, draws = attr(cvd, "rate_draws")), file.path(out, "cv_draws_main.rds"))

set.seed(20251116)
sc <- cvd; R <- attr(cvd, "rate_draws")[, sc$unit, drop = FALSE] / 1000     # rate per birth
B <- sc$births; y <- sc$y
mu <- sweep(R, 2, B, "*")                                                  # draws x units
yrep <- matrix(rpois(length(mu), mu), nrow(mu))                            # Poisson layer
Fy  <- colMeans(ppois(matrix(y,     nrow(mu), ncol(mu), byrow = TRUE), mu))
Fy1 <- colMeans(ppois(matrix(y - 1, nrow(mu), ncol(mu), byrow = TRUE), mu))
pit <- Fy1 + runif(length(y)) * (Fy - Fy1)                                 # randomised PIT

A <- sc |> mutate(obs = 1000 * y / births,
                  lo = 1000 * apply(yrep, 2, quantile, 0.05) / births,
                  hi = 1000 * apply(yrep, 2, quantile, 0.95) / births,
                  pred = 1000 * colMeans(yrep) / births,
                  inside = obs >= lo & obs <= hi)
cat(sprintf("Observed rates inside 90%% predictive interval: %d of %d (%.1f%%)\n",
            sum(A$inside), nrow(A), 100 * mean(A$inside)))

grp <- cut(colMeans(R), quantile(colMeans(R), 0:10 / 10), include.lowest = TRUE, labels = FALSE)
D <- do.call(rbind, lapply(1:10, function(g) {
  j <- grp == g; pooled <- 1000 * rowSums(yrep[, j, drop = FALSE]) / sum(B[j])
  data.frame(g = g, pred = 1000 * mean(rowSums(mu[, j, drop = FALSE])) / sum(B[j]),
             obs = 1000 * sum(y[j]) / sum(B[j]),
             lo = quantile(pooled, 0.025), hi = quantile(pooled, 0.975))
}))
print(D, row.names = FALSE, digits = 4)

th <- theme_minimal(base_size = 9) + theme(panel.grid.minor = element_blank(), legend.position = "bottom",
                                           plot.title = element_text(face = "bold", size = 9))
pa <- ggplot(A, aes(pred, obs)) +
  geom_linerange(aes(ymin = lo, ymax = hi), colour = "grey85", linewidth = 0.2, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, colour = "grey40", linetype = 2, linewidth = 0.4) +
  geom_point(aes(size = births), shape = 21, fill = "#3987e5", colour = "white", alpha = 0.8, stroke = 0.3) +
  scale_size_area(max_size = 6, breaks = c(1000, 5000, 15000), labels = scales::comma, name = "Total births") +
  coord_cartesian(xlim = c(2.5, 5.8), ylim = c(0, 10)) +
  labs(x = "Predicted rate (cross-validated), per 1,000", y = "Observed rate, per 1,000",
       title = "A. Individual LADs (90% predictive interval)") + th
pb <- ggplot(D, aes(pred, obs)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey40", linetype = 2, linewidth = 0.4) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.03, colour = "grey50") +
  geom_point(size = 2, colour = "#1c5cab") +
  labs(x = "Mean predicted rate, per 1,000", y = "Observed pooled rate, per 1,000",
       title = "B. Tenths of predicted rate (95% predictive interval)") + th
ggsave(file.path(out, "figS3_cv_calibration.png"), pa | pb, width = 10, height = 5, dpi = 300, bg = "white")
pc <- ggplot(data.frame(pit), aes(pit)) +
  geom_histogram(aes(y = after_stat(density)), breaks = seq(0, 1, 0.1), fill = "#3987e5", colour = "white") +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey40") +
  labs(x = "Randomised PIT", y = "Density") + th
ggsave(file.path(out, "figS4_pit.png"), pc, width = 5, height = 3.5, dpi = 300, bg = "white")
