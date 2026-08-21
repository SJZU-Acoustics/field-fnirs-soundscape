# A26 — fixation–brain coupling (post-review round 2, 2026-08-14; approved by
# Prof Zhang the same day, prompted by RQ2-2 of the collector's pre-stage-1
# draft). Completes the multimodal coupling matrix: A07 brain–perception,
# A08 brain–context/acoustics, A10 fixation–appraisal/context — fixation→HbO
# had never been consumed.
#
# Declared family (exploration README, fixed before this ran): 2 tests, BH —
#  (1) pair level: dfrontal-HbO ~ dfixation + dLAeq (z-standardised,
#      participant RE; eye-complete pairs, valid-channel frontal mean)
#  (2) block level: frontal HbO ~ side + fixation_bs + fixation_ws +
#      fixation_wb + laeq_wb, design REs + (1 | side_unit); fixation_wb tested
# Convergent/descriptive (no correction): per-channel pair-level slopes (7),
# HbR frontal. Sensitivity: duplicate-vector exclusion, no-LAeq-control.

source(file.path("code", "helpers.R"))
A_DIR <- analysis_output_dir("a26_fixation_brain_coupling")

blocks  <- load_p27_blocks()
samples <- load_p27_samples()

FRONTAL_CH <- paste0("ch", 1:4)

# ---- pair-level table: dfixation, dLAeq, dfrontal (valid-channel mean), ----
# per-channel diffs and HbR frontal diff
samp <- samples %>%
  mutate(
    frontal_hbo = rowMeans(across(all_of(paste0("hbo_", FRONTAL_CH))), na.rm = TRUE),
    frontal_hbr = rowMeans(across(all_of(paste0("hbr_", FRONTAL_CH))), na.rm = TRUE)
  )
wide <- samp %>%
  select(sample_id, participant, site, side,
         fixation_s, laeq, frontal_hbo, frontal_hbr, all_of(HBO_COLS)) %>%
  pivot_wider(names_from = side,
              values_from = c(fixation_s, laeq, frontal_hbo, frontal_hbr,
                              all_of(HBO_COLS)),
              names_glue = "{.value}__{side}")
for (v in c("fixation_s", "laeq", "frontal_hbo", "frontal_hbr", HBO_COLS)) {
  wide[[paste0(v, "__diff")]] <-
    wide[[paste0(v, "__natural")]] - wide[[paste0(v, "__composite")]]
}

pairs_full <- wide %>%
  filter(is.finite(fixation_s__diff), is.finite(frontal_hbo__diff),
         is.finite(laeq__diff))
cat("Pair-level n (eye-complete, frontal, acoustic):", nrow(pairs_full), "\n")

pair_model <- function(d, control_laeq = TRUE, label) {
  d <- d %>%
    mutate(y = as.numeric(scale(frontal_hbo__diff)),
           x = as.numeric(scale(fixation_s__diff)),
           a = as.numeric(scale(laeq__diff)))
  f <- if (control_laeq) y ~ x + a + (1 | participant) else y ~ x + (1 | participant)
  fit <- suppressMessages(suppressWarnings(lmer(f, data = d)))
  co <- summary(fit)$coefficients["x", ]
  r_raw <- cor(d$frontal_hbo__diff, d$fixation_s__diff)
  tibble(model = label, n_pairs = nrow(d), r_bivariate = r_raw,
         beta_std = co[["Estimate"]], se = co[["Std. Error"]],
         p = co[["Pr(>|t|)"]], singular = isSingular(fit))
}

test1 <- pair_model(pairs_full, TRUE, "dfrontal ~ dfixation + dLAeq (primary)")

# sensitivity: exclude the two duplicate-vector samples entirely (n 30 -> 28)
pairs_nodup <- pairs_full %>% filter(!sample_id %in% EYE_DUP_LATER)
sens_pair <- bind_rows(
  pair_model(pairs_full,  FALSE, "no-LAeq-control variant"),
  pair_model(pairs_nodup, TRUE,  "duplicate vectors excluded"))

# ---- block-level model -----------------------------------------------------
bl <- blocks %>%
  filter(fixation_recorded %in% TRUE, fnirs_recorded %in% TRUE,
         acoustic_recorded %in% TRUE) %>%
  rowwise() %>%
  mutate(frontal_hbo = mean(c_across(all_of(paste0("hbo_", FRONTAL_CH)))[
           c_across(all_of(paste0("hbo_", FRONTAL_CH, "_valid")))], na.rm = TRUE),
         frontal_hbr = mean(c_across(all_of(paste0("hbr_", FRONTAL_CH)))[
           c_across(all_of(paste0("hbo_", FRONTAL_CH, "_valid")))], na.rm = TRUE)) %>%
  ungroup() %>%
  filter(is.finite(frontal_hbo)) %>%
  decompose_bw("fixation_s") %>%
  decompose_bw("laeq")
cat("Block-level rows (fixation & fNIRS & acoustic recorded, >=1 valid frontal):",
    nrow(bl), "over", n_distinct(bl$sample_id), "samples\n")

bl <- bl %>% mutate(across(c(fixation_s_bs, fixation_s_ws, fixation_s_wb,
                             laeq_wb), ~as.numeric(scale(.x))))
fit_bl <- p27_lmm(frontal_hbo ~ side + fixation_s_bs + fixation_s_ws +
                    fixation_s_wb + laeq_wb, bl,
                  re = c("(1 | participant)", "(1 | site)",
                         "(1 | sample_id)", "(1 | side_unit)"))
co_bl <- summary(fit_bl)$coefficients
test2 <- tibble(model = "block: frontal ~ side + fixation(bs/ws/wb) + laeq_wb",
                n_rows = nrow(bl), term = "fixation_s_wb",
                beta = co_bl["fixation_s_wb", "Estimate"],
                se = co_bl["fixation_s_wb", "Std. Error"],
                p = co_bl["fixation_s_wb", "Pr(>|t|)"],
                singular = isSingular(fit_bl))

# duplicate-vector sensitivity at block level (drop the duplicated side only)
bl_nodup <- bl %>%
  filter(!(sample_id %in% EYE_DUP_LATER &
             as.character(side) == EYE_DUP_SIDE[as.character(sample_id)]))
fit_bl2 <- p27_lmm(frontal_hbo ~ side + fixation_s_bs + fixation_s_ws +
                     fixation_s_wb + laeq_wb, bl_nodup,
                   re = c("(1 | participant)", "(1 | site)",
                          "(1 | sample_id)", "(1 | side_unit)"))
co_bl2 <- summary(fit_bl2)$coefficients
sens_block <- tibble(model = "block model, duplicate sides excluded",
                     n_rows = nrow(bl_nodup), term = "fixation_s_wb",
                     beta = co_bl2["fixation_s_wb", "Estimate"],
                     se = co_bl2["fixation_s_wb", "Std. Error"],
                     p = co_bl2["fixation_s_wb", "Pr(>|t|)"],
                     singular = isSingular(fit_bl2))

# ---- the declared family, BH over the 2 tested coefficients ----------------
family <- tibble(
  test = c("pair: dfrontal ~ dfixation + dLAeq",
           "block: frontal ~ ... + fixation_wb"),
  n = c(test1$n_pairs, test2$n_rows),
  estimate = c(test1$beta_std, test2$beta),
  se = c(test1$se, test2$se),
  p = c(test1$p, test2$p)) %>%
  mutate(q = bh(p))
write_outcome(family, A_DIR, "fixation_hbo_family.csv")
cat("\nDeclared family (BH):\n"); print(as.data.frame(family), digits = 3)

# ---- descriptive: per-channel pair-level slopes; HbR convergent ------------
per_ch <- map_dfr(HBO_COLS, function(v) {
  d <- wide %>%
    filter(is.finite(fixation_s__diff), is.finite(.data[[paste0(v, "__diff")]]),
           is.finite(laeq__diff)) %>%
    mutate(y = as.numeric(scale(.data[[paste0(v, "__diff")]])),
           x = as.numeric(scale(fixation_s__diff)),
           a = as.numeric(scale(laeq__diff)))
  fit <- suppressMessages(suppressWarnings(lmer(y ~ x + a + (1 | participant), data = d)))
  co <- summary(fit)$coefficients["x", ]
  tibble(channel = v, n_pairs = nrow(d), beta_std = co[["Estimate"]],
         se = co[["Std. Error"]], p_descriptive = co[["Pr(>|t|)"]])
})
hbr <- pairs_full %>%
  filter(is.finite(frontal_hbr__diff)) %>%
  mutate(y = as.numeric(scale(frontal_hbr__diff)),
         x = as.numeric(scale(fixation_s__diff)),
         a = as.numeric(scale(laeq__diff)))
fit_hbr <- suppressMessages(suppressWarnings(lmer(y ~ x + a + (1 | participant), data = hbr)))
co_hbr <- summary(fit_hbr)$coefficients["x", ]
descr <- bind_rows(
  per_ch,
  tibble(channel = "frontal HbR (convergent)", n_pairs = nrow(hbr),
         beta_std = co_hbr[["Estimate"]], se = co_hbr[["Std. Error"]],
         p_descriptive = co_hbr[["Pr(>|t|)"]]))
write_outcome(descr, A_DIR, "per_channel_descriptive.csv")
cat("\nPer-channel / HbR (descriptive):\n"); print(as.data.frame(descr), digits = 3)

sens <- bind_rows(
  sens_pair %>% mutate(level = "pair") %>%
    rename(estimate = beta_std) %>%
    select(level, model, n = n_pairs, estimate, se, p, singular),
  sens_block %>% mutate(level = "block") %>%
    rename(estimate = beta) %>%
    select(level, model, n = n_rows, estimate, se, p, singular))
write_outcome(sens, A_DIR, "sensitivities.csv")
cat("\nSensitivities:\n"); print(as.data.frame(sens), digits = 3)

cat("\nA26 done\n")
