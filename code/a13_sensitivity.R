# A13 — sensitivity and robustness of every headline estimate.
# No new inferential family (declared): everything here is stability evidence.

source(file.path("code", "helpers.R"))
A_DIR <- analysis_output_dir("a13_sensitivity_robustness")

blocks  <- load_p27_blocks()
samples <- load_p27_samples()

# ---- 1. Rank-based (Wilcoxon signed-rank) versions of the headlines --------
wilcox_of <- function(var) {
  p <- load_p27_paired(var); dcol <- paste0(var, "__diff")
  x <- p[[dcol]][is.finite(p[[dcol]])]
  wt <- wilcox.test(x, exact = FALSE)
  tibble(var = var, n = length(x), median_diff = median(x), p_wilcox = wt$p.value)
}
wil <- map_dfr(c("p_iso", "paq_monotonous", "hbo_ch1", "hbo_ch3", "hbo_ch4"), wilcox_of)
write_outcome(wil, A_DIR, "wilcoxon_headlines.csv")

# ---- 2. HbT convergence for the A06 survivors ------------------------------
for (ch in c("ch1", "ch3", "ch4")) {
  samples[[paste0("hbt_", ch)]] <- samples[[paste0("hbo_", ch)]] + samples[[paste0("hbr_", ch)]]
}
hbt <- map_dfr(paste0("hbt_", c("ch1", "ch3", "ch4")), function(v) {
  d <- samples %>% select(sample_id, participant, side, all_of(v)) %>%
    pivot_wider(names_from = side, values_from = all_of(v)) %>%
    mutate(diff = natural - composite) %>% filter(is.finite(diff))
  tt <- t.test(d$diff)
  tibble(var = v, n = nrow(d), mean_diff = mean(d$diff), d_z = mean(d$diff)/sd(d$diff),
         p = tt$p.value)
})
write_outcome(hbt, A_DIR, "hbt_convergence.csv")

# ---- 3. Acoustic-constancy claim under strict pairs ------------------------
# pairs where BOTH sides have >=1 genuinely recorded acoustic block: re-test
# the LAeq null and the P_ISO effect on that subset
ac_ok <- acoustic_rows(blocks) %>% distinct(sample_id, side) %>%
  count(sample_id) %>% filter(n == 2) %>% pull(sample_id)
strict <- map_dfr(c("laeq", "p_iso"), function(v) {
  p <- load_p27_paired(v) %>% filter(sample_id %in% ac_ok)
  dcol <- paste0(v, "__diff")
  x <- p[[dcol]][is.finite(p[[dcol]])]
  tt <- t.test(x)
  tibble(var = v, n_pairs = length(x), mean_diff = mean(x), d_z = mean(x)/sd(x), p = tt$p.value)
})
write_outcome(strict, A_DIR, "strict_acoustic_pairs.csv")

# ---- 4. Six-channel decoding (drop ch2, the most-invalid channel) ----------
pd <- load_p27_paired(c(HBO_COLS))
d6cols <- paste0(setdiff(HBO_COLS, "hbo_ch2"), "__diff")
d6 <- pd %>% filter(if_all(all_of(d6cols), is.finite))
X <- as.matrix(d6[, d6cols]); part <- as.character(d6$participant); n <- nrow(d6)
decode_acc <- function(Xm) {
  accs <- map_dbl(unique(part), function(p) {
    tr <- part != p
    Xtr <- rbind(Xm[tr, , drop = FALSE], -Xm[tr, , drop = FALSE])
    ytr <- c(rep(1, sum(tr)), rep(-1, sum(tr)))
    S <- cov(Xtr) + diag(1e-6, ncol(Xm))
    w <- solve(S, colMeans(Xtr[ytr == 1, , drop = FALSE]) -
                  colMeans(Xtr[ytr == -1, , drop = FALSE]))
    mean((Xm[!tr, , drop = FALSE] %*% w) > 0)
  })
  weights <- table(part)[unique(part)]
  sum(accs * as.numeric(weights)) / sum(weights)
}
acc_obs <- decode_acc(X)
set.seed(271)
acc_null <- replicate(BOOT_N, decode_acc(X * sample(c(-1, 1), n, replace = TRUE)))
dec <- tibble(variant = "6ch_no_ch2", n_pairs = n, accuracy = acc_obs,
              null_q95 = quantile(acc_null, .95), p_perm = mean(acc_null >= acc_obs))
write_outcome(dec, A_DIR, "decoding_6ch.csv")

# ---- 5. Leave-one-participant-out influence on the headlines ---------------
loo <- map_dfr(c("p_iso", "paq_monotonous", "hbo_ch1", "hbo_ch4"), function(v) {
  p <- load_p27_paired(v); dcol <- paste0(v, "__diff")
  d <- p %>% filter(is.finite(.data[[dcol]]))
  ests <- map_dfr(unique(as.character(d$participant)), function(pp) {
    x <- d[[dcol]][d$participant != pp]
    tibble(dropped = pp, est = mean(x), p = t.test(x)$p.value)
  })
  tibble(var = v, est_full = mean(d[[dcol]]), est_min = min(ests$est),
         est_max = max(ests$est), p_max = max(ests$p),
         sign_stable = all(sign(ests$est) == sign(mean(d[[dcol]]))))
})
write_outcome(loo, A_DIR, "loo_influence.csv")

# ---- 6. P_ISO position trend without the first exposure --------------------
d <- blocks %>% filter(exposure_position > 1) %>%
  mutate(yz = as.numeric(scale(p_iso)), pos = as.numeric(scale(exposure_position)))
fit <- p27_lmm(yz ~ side + pos, d, re = c("(1 | participant)", "(1 | site)", "(1 | sample_id)"))
co <- summary(fit)$coefficients["pos", ]
trend <- tibble(variant = "positions_2_to_6", beta = co[["Estimate"]],
                se = co[["Std. Error"]], p = co[["Pr(>|t|)"]])
write_outcome(trend, A_DIR, "position_trend_sensitivity.csv")

print(wil); print(hbt); print(strict); print(dec); print(loo); print(trend)
cat("A13 done\n")
