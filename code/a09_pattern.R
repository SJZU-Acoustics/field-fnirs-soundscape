# A09 — multivariate pattern readout (the estimator/feature-family change the
# protocol's exhaustion audit demands: cross-validated, all channels jointly).
# Declared readouts (../README.md), each native permutation p, 1,000 resamples:
#  (i) leave-participant-out decoding of side from the 7-channel dHbO vector
# (ii) leave-participant-out prediction of dP_ISO from the same vector
# Complete-channel pairs only.

source(file.path("code", "helpers.R"))
A_DIR <- analysis_output_dir("a09_multivariate_pattern")

pd <- load_p27_paired(c("p_iso", HBO_COLS))
dcols <- paste0(HBO_COLS, "__diff")
d <- pd %>% filter(if_all(all_of(dcols), is.finite))
n <- nrow(d)
cat("complete-channel pairs:", n, "\n")

X <- as.matrix(d[, dcols])          # natural - composite, per channel
y_piso <- d$p_iso__diff
part <- as.character(d$participant)

# ---- (i) side decoding -----------------------------------------------------
# Symmetric pairs: observation (X_i, +1) and (-X_i, -1). Classifier: Fisher
# direction w = solve(S) %*% mean-diff, equivalently LDA on the doubled set;
# LOPO: hold out all pairs of one participant.
decode_acc <- function(Xm) {
  accs <- map_dbl(unique(part), function(p) {
    tr <- part != p; te <- !tr
    Xtr <- rbind(Xm[tr, , drop = FALSE], -Xm[tr, , drop = FALSE])
    ytr <- c(rep(1, sum(tr)), rep(-1, sum(tr)))
    S <- cov(Xtr) + diag(1e-6, ncol(Xm))
    w <- solve(S, colMeans(Xtr[ytr == 1, , drop = FALSE]) -
                  colMeans(Xtr[ytr == -1, , drop = FALSE]))
    mean((Xm[te, , drop = FALSE] %*% w) > 0)
  })
  # participant-weighted mean accuracy over held-out pairs
  weights <- table(part)[unique(part)]
  sum(accs * as.numeric(weights)) / sum(weights)
}
acc_obs <- decode_acc(X)
set.seed(271)
acc_null <- replicate(BOOT_N, decode_acc(X * sample(c(-1, 1), n, replace = TRUE)))
p_decode <- mean(acc_null >= acc_obs)

# ---- (ii) dP_ISO prediction ------------------------------------------------
pred_r <- function(yv) {
  preds <- rep(NA_real_, n)
  for (p in unique(part)) {
    tr <- part != p; te <- !tr
    df_tr <- data.frame(y = yv[tr], X[tr, , drop = FALSE])
    fit <- lm(y ~ ., data = df_tr)
    preds[te] <- predict(fit, newdata = data.frame(X[te, , drop = FALSE]))
  }
  cor(preds, yv)
}
r_obs <- pred_r(y_piso)
set.seed(271)
r_null <- replicate(BOOT_N, pred_r(sample(y_piso)))
p_pred <- mean(r_null >= r_obs)

res <- tibble(
  readout = c("side_decoding_LOPO", "dpiso_prediction_LOPO"),
  n_pairs = n, statistic = c(acc_obs, r_obs),
  null_mean = c(mean(acc_null), mean(r_null)),
  null_q95 = c(quantile(acc_null, .95), quantile(r_null, .95)),
  p_perm = c(p_decode, p_pred))
write_outcome(res, A_DIR, "multivariate_readouts.csv")
print(res)
cat("A09 done\n")
