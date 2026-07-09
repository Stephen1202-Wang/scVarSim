#!/usr/bin/env Rscript
# Pool heterozygous allelic ratios across ALL ASARP RDD samples into one figure.
#
# RDD file format (whitespace-delimited, no header):
#   col1 chr | col2 position | col3 ref>alt | col4 (unused) | col5 ref:alt:other counts | col6 (unused)
# Allelic ratio = alt / (ref + alt).
# Heterozygous SNPs: 0 < allelic ratio < 1 (ratio of 0 or 1 = homozygous, excluded).
#
# Usage:
#   Rscript plot_allelic_ratio_all.R [out.pdf]

rdd_dir <- "/u/project/gxxiao/alisonki/projects/ASARP/RDD"

args    <- commandArgs(trailingOnly = TRUE)
outfile <- if (length(args) >= 1) args[1] else "allelic_ratio_all_samples.pdf"

files <- sort(list.files(rdd_dir, pattern = "\\.rdd$", full.names = TRUE))
cat(sprintf("Found %d samples\n", length(files)))

# Heterozygous allelic ratios for one file
get_het_ratio <- function(f) {
  dat <- read.table(f, header = FALSE, stringsAsFactors = FALSE, colClasses = "character")
  counts <- do.call(rbind, strsplit(dat[[5]], ":", fixed = TRUE))
  ref <- as.numeric(counts[, 1])
  alt <- as.numeric(counts[, 2])
  depth <- ref + alt
  ratio <- alt[depth > 0] / depth[depth > 0]
  ratio[ratio > 0 & ratio < 1]   # heterozygous only
}

# Pool ratios from every sample (cache to RDS to avoid re-reading 1000+ files)
cache <- "allelic_ratio_all_pooled.rds"
if (file.exists(cache)) {
  cat(sprintf("Loading pooled ratios from cache: %s\n", cache))
  all_ratio <- readRDS(cache)
} else {
  all_ratio <- unlist(lapply(files, function(f) {
    r <- get_het_ratio(f)
    cat(sprintf("  %-12s %d het sites\n", sub("\\.final.*$", "", basename(f)), length(r)))
    r
  }))
  saveRDS(all_ratio, cache)
}

cat(sprintf("\nTotal pooled heterozygous sites: %d (across %d samples)\n",
            length(all_ratio), length(files)))
cat("Summary of pooled allelic ratio:\n")
print(summary(all_ratio))

n <- length(all_ratio)

## ---- Fit 1: Normal (MLE) ------------------------------------------------
mu       <- mean(all_ratio)
sigma    <- sqrt(mean((all_ratio - mu)^2))   # MLE sd (denominator n)
ll_norm  <- sum(dnorm(all_ratio, mu, sigma, log = TRUE))
cat(sprintf("\nNormal fit (MLE):  mean = %.4f, sd = %.4f\n", mu, sigma))

## ---- Fit 2: Beta --------------------------------------------------------
# Bounded on (0,1): the natural model for an allelic ratio / proportion.
# Method-of-moments gives a closed-form estimate (instant). At n=45M it is
# statistically indistinguishable from MLE, and we refine it with an MLE on a
# random subsample (sub-second) -- avoids MLE-optimising over all 45M points,
# which iterates dbeta() over the full vector dozens of times (very slow).
m <- mean(all_ratio); v <- var(all_ratio)
common  <- m * (1 - m) / v - 1                 # method-of-moments
a_mom <- m * common; b_mom <- (1 - m) * common
set.seed(1)
sub   <- sample(all_ratio, min(n, 2e5))
fit_b <- MASS::fitdistr(sub, dbeta, start = list(shape1 = a_mom, shape2 = b_mom))
a <- fit_b$estimate["shape1"]; b <- fit_b$estimate["shape2"]
ll_beta <- sum(dbeta(all_ratio, a, b, log = TRUE))   # full-data logLik: one pass
cat(sprintf("Beta fit (MLE on %d subsample): shape1(alpha) = %.4f, shape2(beta) = %.4f  (mean = %.4f)\n",
            length(sub), a, b, a / (a + b)))

## ---- Goodness-of-fit ----------------------------------------------------
gof <- function(name, ll, k, cdf) {
  aic <- -2 * ll + 2 * k
  bic <- -2 * ll + k * log(n)
  # KS distance D (computed on the full data; p-values are meaningless at n=45M,
  # so we report only the distance statistic). Sample for a usable p-value.
  set.seed(1)
  s   <- sample(all_ratio, min(n, 1e5))
  ks  <- suppressWarnings(ks.test(s, cdf))
  # Binned chi-square over the 50 histogram bins
  br  <- seq(0, 1, by = 0.02)
  obs <- as.numeric(table(cut(all_ratio, br, include.lowest = TRUE)))
  exp <- n * diff(cdf(br))
  chi <- sum((obs - exp)^2 / exp)
  df  <- length(obs) - 1 - k
  cat(sprintf("  %-7s  logLik=%.0f  AIC=%.0f  BIC=%.0f  KS-D=%.4f (p=%.3g, n=1e5)  chi2=%.0f (df=%d)\n",
              name, ll, aic, bic, ks$statistic, ks$p.value, chi, df))
  c(aic = aic, bic = bic, ksD = unname(ks$statistic))
}

cat("\nGoodness-of-fit:\n")
g_norm <- gof("Normal", ll_norm, 2, function(q) pnorm(q, mu, sigma))
g_beta <- gof("Beta",   ll_beta, 2, function(q) pbeta(q, a, b))
better <- if (g_beta["aic"] < g_norm["aic"]) "Beta" else "Normal"
cat(sprintf("  -> lower AIC/BIC = better fit: %s\n", better))

## ---- Plot: histogram + both fitted curves -------------------------------
pdf(outfile, width = 8, height = 6)
hist(all_ratio,
     breaks = seq(0, 1, by = 0.02),
     freq = FALSE,
     col = "steelblue", border = "white",
     main = sprintf("Pooled heterozygous allelic ratio (%d samples)", length(files)),
     xlab = "Allelic ratio  alt/(ref+alt)",
     ylab = "Density")
xs <- seq(0.001, 0.999, length.out = 500)
lines(xs, dbeta(xs, a, b),                  col = "darkgreen", lwd = 2.5)        # Beta (bounded)
lines(xs, dnorm(xs, mean = mu, sd = sigma), col = "darkred",   lwd = 2, lty = 3) # Normal (ref)
legend("topright",
       legend = c(sprintf("n = %d sites", n),
                  sprintf("Beta(%.2f, %.2f)  AIC=%.0f", a, b, g_beta["aic"]),
                  sprintf("Normal(%.3f, %.3f)  AIC=%.0f", mu, sigma, g_norm["aic"])),
       col = c(NA, "darkgreen", "darkred"),
       lwd = c(NA, 2.5, 2), lty = c(NA, 1, 3), bty = "n")

## ---- Page 2: Normal Q-Q sanity check ------------------------------------
# Quantile-grid Q-Q (fast at n=45M): sample quantiles vs theoretical Normal
# quantiles. Points on the diagonal => data ~ Normal.
probs  <- ppoints(500)                       # 500 evenly spaced probabilities
sample_q <- quantile(all_ratio, probs, names = FALSE)
theo_q   <- qnorm(probs, mu, sigma)

# Numeric sanity check: sample vs theoretical Normal quantiles at key probs
cat("\nNormal Q-Q numeric check (sample vs theoretical quantiles):\n")
key <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)
tab <- data.frame(prob = key,
                  sample      = round(quantile(all_ratio, key, names = FALSE), 4),
                  theoretical = round(qnorm(key, mu, sigma), 4))
tab$diff <- round(tab$sample - tab$theoretical, 4)
print(tab, row.names = FALSE)

qqplot_fun <- function() {
plot(theo_q, sample_q,
     pch = 19, cex = 0.5, col = "steelblue",
     main = "Normal Q-Q plot (sanity check)",
     xlab = sprintf("Theoretical Normal(%.3f, %.3f) quantiles", mu, sigma),
     ylab = "Sample quantiles of allelic ratio")
abline(0, 1, col = "red", lwd = 2)           # y = x reference line
legend("topleft",
       legend = c("quantile grid (500 pts)", "y = x (perfect Normal)"),
       col = c("steelblue", "red"), pch = c(19, NA), lwd = c(NA, 2), bty = "n")
}
qqplot_fun()
dev.off()

# Also write the Q-Q plot as a standalone PNG (viewable without a PDF renderer)
qq_png <- sub("\\.pdf$", "_normal_qq.png", outfile)
png(qq_png, width = 800, height = 800, res = 110)
qqplot_fun()
dev.off()
cat(sprintf("Saved Normal Q-Q PNG to: %s\n", normalizePath(qq_png)))

cat(sprintf("\nSaved plot to: %s\n", normalizePath(outfile)))
