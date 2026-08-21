# A25 — post-review sensitivities (added after the AI-review triage of the
# first full draft, 2026-08-13; approved by Prof Zhang the same day).
# These are equivalence/stability checks on already-reported contrasts and
# carry no multiple-comparison correction (SI Methods states this). None was
# declared before its batch; all are labelled post-review in the manuscript.
#
# Components:
#  (1) TOST equivalence on the four acoustic indices (margin ±1.5 dB(A))
#  (2) site random-intercept sensitivity on the primary paired contrasts
#  (3) matched-sample decoding: PAQ classifier on the same 30 complete-channel
#      pairs as the HbO classifier + paired McNemar test of the accuracies
#  (4) matched-channel (complete-4) frontal-mean sensitivity
#  (5) first-session-only subset of the first-block between-group analysis
#  (6) 95% CI for the pair-level dose slope (A04 between_side estimator)
#  (7) corrected reliability-table statistics for Table S1 (units = side-units
#      with >= 2 recorded blocks; mean k WITHIN those units — fixes A01's
#      mean_k denominator and the caption's ">= 1" misdescription)

source(file.path("code", "helpers.R"))
A_DIR <- analysis_output_dir("a25_post_review_sensitivities")

blocks  <- load_p27_blocks()
samples <- load_p27_samples()

# ---- (1) TOST equivalence on acoustic indices ------------------------------
# Margin: +/- 1.5 dB(A), half the ~3 dB(A) level change conventionally taken
# as just noticeable for broadband environmental sound.
MARGIN <- 1.5
tost <- map_dfr(ACOUSTIC_COLS, function(v) {
  p <- load_p27_paired(v)
  d <- p[[paste0(v, "__diff")]]; d <- d[!is.na(d)]
  n <- length(d); m <- mean(d); se <- sd(d) / sqrt(n)
  t_lo <- (m - (-MARGIN)) / se   # H0: diff <= -margin
  t_hi <- (m - MARGIN) / se      # H0: diff >= +margin
  p_lo <- pt(t_lo, n - 1, lower.tail = FALSE)
  p_hi <- pt(t_hi, n - 1, lower.tail = TRUE)
  ci90 <- m + qt(c(.05, .95), n - 1) * se
  tibble(var = v, n = n, mean_diff = m, se = se,
         ci90_lo = ci90[1], ci90_hi = ci90[2],
         margin = MARGIN, p_tost = max(p_lo, p_hi))
})
write_outcome(tost, A_DIR, "tost_acoustic_equivalence.csv")
cat("TOST done\n"); print(as.data.frame(tost), digits = 3)

# ---- (2) site random-intercept sensitivity ---------------------------------
HEADLINES <- c("paq_monotonous", "p_iso", "hbo_ch1", "hbo_ch3", "hbo_ch4", "laeq")
site_re <- map_dfr(HEADLINES, function(v) {
  p <- load_p27_paired(v)
  dcol <- paste0(v, "__diff")
  d <- p %>% filter(!is.na(.data[[dcol]]))
  base <- side_contrast(v, boot = FALSE)
  fit <- suppressMessages(suppressWarnings(
    lmer(reformulate("1 + (1 | participant) + (1 | site)", response = dcol), data = d)))
  co <- summary(fit)$coefficients[1, ]
  tibble(var = v, n_pairs = nrow(d),
         est_participant_re = base$lmm_est, p_participant_re = base$lmm_p,
         est_plus_site_re = co[["Estimate"]], se_plus_site_re = co[["Std. Error"]],
         p_plus_site_re = co[["Pr(>|t|)"]], singular = isSingular(fit))
})
write_outcome(site_re, A_DIR, "site_re_sensitivity.csv")
cat("site RE done\n"); print(as.data.frame(site_re), digits = 3)

# ---- (3) matched-sample decoding + paired accuracy test --------------------
pd <- load_p27_paired(c(PAQ_COLS, HBO_COLS))
hbo_dcols <- paste0(HBO_COLS, "__diff")
paq_dcols <- paste0(PAQ_COLS, "__diff")
d30 <- pd %>% filter(if_all(all_of(hbo_dcols), is.finite))
stopifnot(nrow(d30) == 30)
part30 <- as.character(d30$participant)

# A09's Fisher-direction LOPO classifier, returning per-pair correctness.
decode_correct <- function(Xm, part) {
  n <- nrow(Xm); correct <- rep(NA, n)
  for (p in unique(part)) {
    tr <- part != p; te <- !tr
    Xtr <- rbind(Xm[tr, , drop = FALSE], -Xm[tr, , drop = FALSE])
    ytr <- c(rep(1, sum(tr)), rep(-1, sum(tr)))
    S <- cov(Xtr) + diag(1e-6, ncol(Xm))
    w <- solve(S, colMeans(Xtr[ytr == 1, , drop = FALSE]) -
                  colMeans(Xtr[ytr == -1, , drop = FALSE]))
    correct[te] <- (Xm[te, , drop = FALSE] %*% w) > 0
  }
  correct
}
# participant-weighted accuracy, as in A09
wacc <- function(correct, part) {
  accs <- tapply(correct, part, mean)
  w <- table(part)[names(accs)]
  sum(accs * as.numeric(w)) / sum(w)
}

X_hbo <- as.matrix(d30[, hbo_dcols])
X_paq <- as.matrix(d30[, paq_dcols])
c_hbo <- decode_correct(X_hbo, part30)
c_paq <- decode_correct(X_paq, part30)
acc_hbo <- wacc(c_hbo, part30); acc_paq <- wacc(c_paq, part30)

set.seed(271)
null_paq <- replicate(BOOT_N, {
  flip <- sample(c(-1, 1), nrow(X_paq), replace = TRUE)
  wacc(decode_correct(X_paq * flip, part30), part30)
})
p_paq30 <- mean(null_paq >= acc_paq)

# paired comparison on identical held-out pairs: exact McNemar on discordants
b <- sum(c_hbo & !c_paq); c_ <- sum(!c_hbo & c_paq)
p_mcnemar <- binom.test(b, b + c_, 0.5)$p.value

matched <- tibble(
  readout = c("side_from_dHbO_30", "side_from_dPAQ_30", "accuracy_difference"),
  n_pairs = 30,
  statistic = c(acc_hbo, acc_paq, acc_hbo - acc_paq),
  null_q95 = c(NA, quantile(null_paq, .95), NA),
  p = c(NA, p_paq30, p_mcnemar),
  note = c("A09 value, unchanged", "sign-flip permutation null",
           sprintf("McNemar exact on discordant pairs (%d vs %d)", b, c_)))
write_outcome(matched, A_DIR, "matched_decoding.csv")
cat("matched decoding done\n"); print(as.data.frame(matched), digits = 3)

# ---- (4) matched-channel (complete-4) frontal mean -------------------------
FRONTAL <- paste0("hbo_ch", c(1, 3, 2, 4))
pf <- load_p27_paired(FRONTAL)
c4 <- pf %>% filter(if_all(ends_with("__diff"), is.finite)) %>%
  mutate(frontal_diff = rowMeans(across(all_of(paste0(FRONTAL, "__diff")))),
         left_diff  = rowMeans(across(all_of(paste0(c("hbo_ch1", "hbo_ch3"), "__diff")))),
         right_diff = rowMeans(across(all_of(paste0(c("hbo_ch2", "hbo_ch4"), "__diff")))))
mk_row <- function(x, label) {
  tt <- t.test(x)
  tibble(measure = label, n = length(x), mean_diff = mean(x),
         d_z = mean(x) / sd(x), p = tt$p.value)
}
# varying-composition version for reference (valid channels per sample-side,
# as in the manuscript's regional definition)
samp_frontal <- samples %>%
  mutate(frontal = rowMeans(across(all_of(paste0("hbo_ch", 1:4))), na.rm = TRUE)) %>%
  select(sample_id, participant, side, frontal) %>%
  pivot_wider(names_from = side, values_from = frontal) %>%
  mutate(diff = natural - composite) %>% filter(is.finite(diff))
c4_res <- bind_rows(
  mk_row(samp_frontal$diff, "frontal mean, valid channels (manuscript)"),
  mk_row(c4$frontal_diff,  "frontal mean, matched complete-4 channels"),
  mk_row(c4$left_diff,     "left-frontal mean, complete-4"),
  mk_row(c4$right_diff,    "right-frontal mean, complete-4"),
  mk_row(c4$left_diff - c4$right_diff, "left minus right (diff-of-diffs), complete-4"))
write_outcome(c4_res, A_DIR, "matched_channel_frontal.csv")
cat("matched-channel frontal done\n"); print(as.data.frame(c4_res), digits = 3)

# ---- (5) first-session-only subset of the first-block analysis -------------
# Selection rule: first session by (session_date, session_time). The deposit
# omits session_time, so ordering falls back to session_date alone — identical
# here, because no participant's two sessions share a date with different
# times (each participant was run twice in the same time slot).
pairing <- load_p27_pairing()
by <- intersect(c("session_date", "session_time"), names(pairing))
first_sessions <- pairing %>%
  group_by(participant) %>%
  arrange(across(all_of(by)), .by_group = TRUE) %>%
  slice(1) %>% ungroup() %>% pull(sample_id)

pos1 <- blocks %>% filter(exposure_position == 1L)
frontal_pos1 <- pos1 %>% fnirs_rows() %>%
  rowwise() %>%
  mutate(frontal = mean(c_across(all_of(paste0("hbo_ch", 1:4)))[
    c_across(all_of(paste0("hbo_ch", 1:4, "_valid")))], na.rm = TRUE)) %>%
  ungroup() %>% select(sample_id, frontal)
pos1 <- pos1 %>% left_join(frontal_pos1, by = "sample_id")

first_block <- function(d, v) {
  dd <- d %>% filter(!is.na(.data[[v]]))
  g <- dd %>% group_by(first_side) %>%
    summarise(m = mean(.data[[v]]), n = n(), .groups = "drop")
  tt <- t.test(reformulate("first_side", response = v), data = dd)
  tibble(outcome = v,
         n_natural_first = g$n[g$first_side == "natural"],
         n_composite_first = g$n[g$first_side == "composite"],
         mean_natural_first = g$m[g$first_side == "natural"],
         mean_composite_first = g$m[g$first_side == "composite"],
         estimate = g$m[g$first_side == "natural"] - g$m[g$first_side == "composite"],
         p_welch = tt$p.value)
}
fs <- pos1 %>% filter(sample_id %in% first_sessions)
firstsess <- bind_rows(
  first_block(fs, "frontal") %>% mutate(outcome = "frontal HbO"),
  first_block(fs, "paq_monotonous"),
  first_block(fs, "p_iso")) %>%
  mutate(subset = "first sessions only (n = 27 participants)")
write_outcome(firstsess, A_DIR, "first_session_subset.csv")
cat("first-session subset done\n"); print(as.data.frame(firstsess), digits = 3)

# ---- (6) dose-slope 95% CI (A04 between_side estimator, replicated) --------
side_means <- blocks %>% group_by(sample_id, participant, side) %>%
  summarise(p_iso = mean(p_iso), .groups = "drop")
nat_means <- context_rows(blocks) %>% group_by(sample_id, side) %>%
  summarise(naturalness = mean(naturalness, na.rm = TRUE), .groups = "drop")
dw <- side_means %>% left_join(nat_means, by = c("sample_id", "side")) %>%
  pivot_wider(names_from = side, values_from = c(p_iso, naturalness)) %>%
  mutate(d_y = p_iso_natural - p_iso_composite,
         d_x = naturalness_natural - naturalness_composite) %>%
  filter(!is.na(d_x), !is.na(d_y)) %>%
  mutate(xz = as.numeric(scale(d_x)), yz = as.numeric(scale(d_y)))
fit <- suppressMessages(suppressWarnings(lmer(yz ~ xz + (1 | participant), data = dw)))
co <- summary(fit)$coefficients["xz", ]
ci <- suppressMessages(confint(fit, parm = "xz", method = "Wald"))
dose <- tibble(beta_std = co[["Estimate"]], se = co[["Std. Error"]],
               p = co[["Pr(>|t|)"]], ci95_lo = ci[1], ci95_hi = ci[2],
               n = nrow(dw), singular = isSingular(fit))
write_outcome(dose, A_DIR, "dose_slope_ci.csv")
cat("dose slope CI done\n"); print(as.data.frame(dose), digits = 3)

# ---- (7) corrected reliability-table statistics ----------------------------
unit_stats <- function(d, flagged, value_col, label) {
  g <- flagged %>% group_by(sample_id, side) %>%
    summarise(k = sum(!is.na(.data[[value_col]])), .groups = "drop")
  tibble(measure = label,
         units_ge1 = sum(g$k >= 1), units_ge2 = sum(g$k >= 2),
         mean_k_within_ge2 = mean(g$k[g$k >= 2]))
}
rel <- bind_rows(
  map_dfr(c("p_iso", "e_iso", PAQ_COLS, SRC_COLS),
          ~unit_stats(blocks, blocks, .x, .x)),
  map_dfr(ACOUSTIC_COLS, ~unit_stats(blocks, acoustic_rows(blocks), .x, .x)),
  map_dfr(CONTEXT_COLS[c(1, 2, 3)],
          ~unit_stats(blocks, context_rows(blocks), .x, .x)),
  unit_stats(blocks, fixation_rows(blocks), "fixation_s", "fixation_s"),
  map_dfr(CHANNELS$ch, function(ch)
    unit_stats(blocks, fnirs_channel_rows(blocks, ch), paste0("hbo_", ch),
               paste0("hbo_", ch))))
write_outcome(rel, A_DIR, "reliability_units_corrected.csv")
cat("reliability units done\n"); print(as.data.frame(rel), digits = 4)

cat("A25 done\n")
