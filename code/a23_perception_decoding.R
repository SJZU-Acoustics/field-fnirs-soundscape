# A23 — multivariate perception decoding: A09's estimator applied to appraisal.
# The HbO claim has a cross-validated multivariate readout (76.7%, p = 0.009);
# the perception claim has never had one. Also: does physiology add anything to
# perception for side discrimination (multimodal vector)?
# Declared readouts, each native permutation p (1,000 resamples, A09 estimator,
# symmetric-pair Fisher/LDA with LOPO):
#  (i) side from the 8-item dPAQ profile (all 40 pairs)
#  (ii) side from the multimodal dPAQ + 7-channel dHbO vector (complete-channel
#       pairs, n = 30 — same set as A09)

source(file.path("code", "helpers.R"))
A_DIR <- analysis_output_dir("a23_perception_decoding")

decode_acc <- function(Xm, part) {
  accs <- map_dbl(unique(part), function(p) {
    tr <- part != p; te <- !tr
    Xtr <- rbind(Xm[tr, , drop = FALSE], -Xm[tr, , drop = FALSE])
    ytr <- c(rep(1, sum(tr)), rep(-1, sum(tr)))
    S <- cov(Xtr) + diag(1e-6, ncol(Xm))
    w <- solve(S, colMeans(Xtr[ytr == 1, , drop = FALSE]) -
                  colMeans(Xtr[ytr == -1, , drop = FALSE]))
    mean((Xm[te, , drop = FALSE] %*% w) > 0)
  })
  weights <- table(part)[unique(part)]
  sum(accs * as.numeric(weights)) / sum(weights)
}

run_decode <- function(X, part, label) {
  n <- nrow(X)
  acc_obs <- decode_acc(X, part)
  set.seed(271)
  acc_null <- replicate(BOOT_N, decode_acc(X * sample(c(-1, 1), n, replace = TRUE), part))
  tibble(readout = label, n_pairs = n, accuracy = acc_obs,
         null_mean = mean(acc_null), null_q95 = quantile(acc_null, .95),
         p_perm = mean(acc_null >= acc_obs))
}

# (i) PAQ profile, all pairs
pd_paq <- load_p27_paired(PAQ_COLS)
paq_d <- paste0(PAQ_COLS, "__diff")
d1 <- pd_paq %>% filter(if_all(all_of(paq_d), is.finite))
X1 <- as.matrix(d1[, paq_d]); part1 <- as.character(d1$participant)

# (ii) multimodal, complete-channel pairs (same n as A09)
pd_mm <- load_p27_paired(c(PAQ_COLS, HBO_COLS))
mm_d <- c(paq_d, paste0(HBO_COLS, "__diff"))
d2 <- pd_mm %>% filter(if_all(all_of(mm_d), is.finite))
X2 <- as.matrix(d2[, mm_d]); part2 <- as.character(d2$participant)

res <- bind_rows(
  run_decode(X1, part1, "side_from_dPAQ_profile_LOPO"),
  run_decode(X2, part2, "side_from_dPAQ_plus_dHbO_LOPO"))
write_outcome(res, A_DIR, "perception_decoding_readouts.csv")

# descriptive: PAQ-only accuracy on the SAME complete-channel pairs, so the
# multimodal gain is read against a like-for-like baseline
d3 <- pd_mm %>% filter(if_all(all_of(mm_d), is.finite))
X3 <- as.matrix(d3[, paq_d])
des <- run_decode(X3, as.character(d3$participant), "paq_only_on_complete_pairs_descriptive")
write_outcome(des, A_DIR, "paq_on_complete_pairs_descriptive.csv")

print(res); print(des)
cat("A23 done\n")
