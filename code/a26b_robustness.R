# A26b — robustness of the unexpected pair-level fixation->frontal-HbO
# coupling (falsify-first pass; descriptive, no new family): rank-based
# estimate, leave-one-out influence, first_side covariate, site RE,
# participant-cluster bootstrap CI, and the raw pairs for outlier inspection.

source(file.path("code", "helpers.R"))
A_DIR <- analysis_output_dir("a26_fixation_brain_coupling")

samples <- load_p27_samples()
samp <- samples %>%
  mutate(frontal_hbo = rowMeans(across(all_of(paste0("hbo_ch", 1:4))), na.rm = TRUE))
wide <- samp %>%
  select(sample_id, participant, site, first_side, side,
         fixation_s, laeq, frontal_hbo) %>%
  pivot_wider(names_from = side, values_from = c(fixation_s, laeq, frontal_hbo),
              names_glue = "{.value}__{side}") %>%
  mutate(d_fix = fixation_s__natural - fixation_s__composite,
         d_laeq = laeq__natural - laeq__composite,
         d_hbo = frontal_hbo__natural - frontal_hbo__composite) %>%
  filter(is.finite(d_fix), is.finite(d_hbo), is.finite(d_laeq))
stopifnot(nrow(wide) == 30)

d <- wide %>% mutate(y = as.numeric(scale(d_hbo)), x = as.numeric(scale(d_fix)),
                     a = as.numeric(scale(d_laeq)))

# (1) rank-based
sp <- suppressWarnings(cor.test(d$d_fix, d$d_hbo, method = "spearman"))
kd <- suppressWarnings(cor.test(d$d_fix, d$d_hbo, method = "kendall"))

# (2) leave-one-out influence on the primary model's x coefficient
loo <- map_dfr(seq_len(nrow(d)), function(i) {
  fit <- suppressMessages(suppressWarnings(
    lmer(y ~ x + a + (1 | participant), data = d[-i, ])))
  co <- summary(fit)$coefficients["x", ]
  tibble(dropped = as.character(d$sample_id[i]),
         beta = co[["Estimate"]], p = co[["Pr(>|t|)"]])
})

# (3) covariate / RE variants
fit_fs <- suppressMessages(suppressWarnings(
  lmer(y ~ x + a + first_side + (1 | participant), data = d)))
fit_site <- suppressMessages(suppressWarnings(
  lmer(y ~ x + a + (1 | participant) + (1 | site), data = d)))
co_fs <- summary(fit_fs)$coefficients["x", ]
co_site <- summary(fit_site)$coefficients["x", ]

# (4) participant-cluster bootstrap CI of the standardised slope
bstat <- function(dd) {
  fit <- suppressMessages(suppressWarnings(lm(y ~ x + a, data = dd)))
  coef(fit)[["x"]]
}
bs <- boot_by_participant(d, bstat)
ci <- quantile(bs, c(.025, .975), names = FALSE)

rob <- tibble(
  check = c("Spearman rho", "Kendall tau", "LOO beta range", "LOO max p",
            "first_side covariate", "site RE added",
            "cluster-bootstrap beta 95% CI"),
  value = c(sprintf("rho = %.3f, p = %.4f", sp$estimate, sp$p.value),
            sprintf("tau = %.3f, p = %.4f", kd$estimate, kd$p.value),
            sprintf("[%.3f, %.3f]", min(loo$beta), max(loo$beta)),
            sprintf("%.4f (dropping %s)", max(loo$p), loo$dropped[which.max(loo$p)]),
            sprintf("beta = %.3f, p = %.4f", co_fs[["Estimate"]], co_fs[["Pr(>|t|)"]]),
            sprintf("beta = %.3f, p = %.4f", co_site[["Estimate"]], co_site[["Pr(>|t|)"]]),
            sprintf("[%.3f, %.3f]", ci[1], ci[2])))
write_outcome(rob, A_DIR, "robustness_addendum.csv")
print(as.data.frame(rob))

# (5) raw pairs for outlier inspection
write_outcome(wide %>% select(sample_id, participant, site, first_side,
                              d_fix, d_laeq, d_hbo) %>% arrange(d_fix),
              A_DIR, "pairs_raw.csv")
cat("\nRaw pairs (sorted by d_fix):\n")
print(as.data.frame(wide %>% select(sample_id, d_fix, d_hbo) %>% arrange(d_fix)),
      digits = 3)
cat("\nA26b done\n")
