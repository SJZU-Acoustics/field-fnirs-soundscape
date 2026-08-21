# A05 — fNIRS signal anatomy.
# Inferential family (declared): 7 one-sample tests — does in-situ audiovisual
# exposure shift HbO from the masked resting baseline, per channel — BH; HbR
# convergent. Everything else descriptive: distributions/tails, channel
# intercorrelations, HbO-HbR coupling, regional coherence.

source(file.path("code", "helpers.R"))
A_DIR <- analysis_output_dir("a05_fnirs_signal_anatomy")

blocks  <- load_p27_blocks()
samples <- load_p27_samples()

# ---- 1. Activation vs baseline: one-sample tests on sample-side means ------
# unit = sample x side (80, minus invalid); participant-cluster LMM
act <- map_dfr(CHANNELS$ch, function(ch) {
  v <- paste0("hbo_", ch)
  d <- samples %>% filter(!is.na(.data[[v]]))
  fit <- suppressMessages(suppressWarnings(
    lmer(reformulate("1 + (1 | participant)", response = v), data = d)))
  co <- summary(fit)$coefficients[1, ]
  tibble(channel = ch, chromophore = "HbO", n_units = nrow(d),
         mean = mean(d[[v]]), est = co[["Estimate"]], se = co[["Std. Error"]],
         p = co[["Pr(>|t|)"]], singular = isSingular(fit))
}) %>% mutate(q = bh(p), sig_q05 = q < 0.05)
act_hbr <- map_dfr(CHANNELS$ch, function(ch) {
  v <- paste0("hbr_", ch)
  d <- samples %>% filter(!is.na(.data[[paste0("hbo_", ch)]]))
  fit <- suppressMessages(suppressWarnings(
    lmer(reformulate("1 + (1 | participant)", response = v), data = d)))
  co <- summary(fit)$coefficients[1, ]
  tibble(channel = ch, chromophore = "HbR", n_units = nrow(d),
         mean = mean(d[[v]], na.rm = TRUE), est = co[["Estimate"]],
         se = co[["Std. Error"]], p = co[["Pr(>|t|)"]], singular = isSingular(fit))
})
write_outcome(bind_rows(act, act_hbr), A_DIR, "activation_vs_baseline.csv")

# ---- 2. Distribution anatomy (block level, valid channels) -----------------
dist <- map_dfr(CHANNELS$ch, function(ch) {
  v <- paste0("hbo_", ch)
  x <- fnirs_channel_rows(blocks, ch)[[v]]
  tibble(channel = ch, n = length(x), mean = mean(x), sd = sd(x),
         median = median(x), q01 = quantile(x, .01), q99 = quantile(x, .99),
         skew = mean(((x - mean(x))/sd(x))^3),
         kurtosis = mean(((x - mean(x))/sd(x))^4),
         frac_abs_gt3sd = mean(abs(x - mean(x)) > 3 * sd(x)))
})
write_outcome(dist, A_DIR, "hbo_distribution_anatomy.csv")

# ---- 3. Channel intercorrelations at sample-side level ---------------------
sm <- samples %>% select(all_of(paste0("hbo_", CHANNELS$ch)))
cc <- cor(sm, use = "pairwise.complete.obs")
write_outcome(as_tibble(cc, rownames = "ch_i"), A_DIR, "hbo_channel_correlations.csv")
# region-pair summary
pairs_long <- as_tibble(cc, rownames = "ch_i") %>%
  pivot_longer(-ch_i, names_to = "ch_j", values_to = "r") %>%
  filter(ch_i < ch_j) %>%
  mutate(reg_i = CHANNELS$region[match(sub("hbo_", "", ch_i), CHANNELS$ch)],
         reg_j = CHANNELS$region[match(sub("hbo_", "", ch_j), CHANNELS$ch)],
         pair_type = ifelse(reg_i == reg_j, paste0("within_", reg_i), "between_region"))
write_outcome(pairs_long %>% group_by(pair_type) %>%
                summarise(mean_r = mean(r), n_pairs = n()), A_DIR, "regional_coherence.csv")

# ---- 4. HbO-HbR coupling per channel (sample-side level) -------------------
coup <- map_dfr(CHANNELS$ch, function(ch) {
  d <- samples %>% filter(!is.na(.data[[paste0("hbo_", ch)]]))
  tibble(channel = ch, n = nrow(d),
         r_hbo_hbr = cor(d[[paste0("hbo_", ch)]], d[[paste0("hbr_", ch)]],
                         use = "complete.obs"))
})
write_outcome(coup, A_DIR, "hbo_hbr_coupling.csv")

print(bind_rows(act, act_hbr) %>% select(channel, chromophore, n_units, est, se, p, q))
print(dist %>% select(channel, sd, skew, kurtosis, frac_abs_gt3sd))
print(coup)
cat("A05 done\n")
