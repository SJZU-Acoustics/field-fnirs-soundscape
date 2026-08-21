# A17 — formal frontal lateralisation of the natural-side HbO drop. A06 showed
# ch1 (L) and ch4 (R) both survive, but "bilateral" was never tested as a
# difference-of-differences. The montage leaves frontal the ONLY region where
# a hemisphere contrast exists (codebook §3), so run it formally.
# Family (declared): 3 tests — side contrast on left-frontal dHbO (ch1+ch3),
# on right-frontal dHbO (ch2+ch4), and the hemisphere diff-of-diffs
# (dLeft - dRight); pair level, participant RE; BH. Block-level side x
# hemisphere LMM as stability, not a second family.

source(file.path("code", "helpers.R"))
A_DIR <- analysis_output_dir("a17_frontal_lateralisation")

samples <- load_p27_samples() %>%
  rowwise() %>%
  mutate(hbo_L = mean(c(hbo_ch1, hbo_ch3), na.rm = TRUE),
         hbo_R = mean(c(hbo_ch2, hbo_ch4), na.rm = TRUE)) %>%
  ungroup()

pw <- samples %>%
  select(sample_id, participant, site, side, hbo_L, hbo_R) %>%
  pivot_wider(names_from = side, values_from = c(hbo_L, hbo_R)) %>%
  mutate(d_L = hbo_L_natural - hbo_L_composite,
         d_R = hbo_R_natural - hbo_R_composite,
         d_LR = d_L - d_R)

# pair-level test on a difference vector: t + participant-RE LMM + cluster boot
pair_test <- function(d, vcol, label) {
  dd <- d %>% filter(is.finite(.data[[vcol]]))
  v <- dd[[vcol]]
  tt <- t.test(v)
  fit <- suppressMessages(suppressWarnings(
    lmer(reformulate("1 + (1 | participant)", response = vcol), data = dd)))
  lmm <- summary(fit)$coefficients[1, ]
  set.seed(271)
  parts <- split(seq_len(nrow(dd)), dd$participant)
  bm <- replicate(BOOT_N, mean(v[unlist(parts[sample(names(parts), length(parts), TRUE)])]))
  tibble(test = label, n_pairs = nrow(dd), mean = mean(v), sd = sd(v),
         d_z = mean(v) / sd(v), t = unname(tt$statistic), p_t = tt$p.value,
         lmm_p = lmm[["Pr(>|t|)"]], lmm_singular = isSingular(fit),
         boot_lo = quantile(bm, .025, names = FALSE),
         boot_hi = quantile(bm, .975, names = FALSE))
}

fam <- bind_rows(pair_test(pw, "d_L", "left frontal (ch1+ch3)"),
                 pair_test(pw, "d_R", "right frontal (ch2+ch4)"),
                 pair_test(pw, "d_LR", "hemisphere diff-of-diffs (L-R)")) %>%
  mutate(q = bh(p_t))
write_outcome(fam, A_DIR, "frontal_lateralisation_tests.csv")

# stability: block-level LMM with side x hemisphere interaction
blocks <- load_p27_blocks()
blk <- blocks %>%
  filter(fnirs_recorded %in% TRUE) %>%
  rowwise() %>%
  mutate(hbo_L = ifelse(hbo_ch1_valid | hbo_ch3_valid,
                        mean(c(if (isTRUE(hbo_ch1_valid)) hbo_ch1,
                               if (isTRUE(hbo_ch3_valid)) hbo_ch3), na.rm = TRUE), NA_real_),
         hbo_R = ifelse(hbo_ch2_valid | hbo_ch4_valid,
                        mean(c(if (isTRUE(hbo_ch2_valid)) hbo_ch2,
                               if (isTRUE(hbo_ch4_valid)) hbo_ch4), na.rm = TRUE), NA_real_)) %>%
  ungroup() %>%
  select(sample_id, participant, site, side, block, side_unit, hbo_L, hbo_R) %>%
  pivot_longer(c(hbo_L, hbo_R), names_to = "hemisphere", values_to = "hbo",
               names_prefix = "hbo_") %>%
  filter(is.finite(hbo))
fit <- suppressMessages(suppressWarnings(lmer(
  hbo ~ side * hemisphere + (1 | participant) + (1 | site) +
    (1 | sample_id) + (1 | side_unit), data = blk, REML = TRUE)))
co <- summary(fit)$coefficients
stab <- tibble(term = rownames(co), est = co[, "Estimate"],
               se = co[, "Std. Error"], p = co[, "Pr(>|t|)"]) %>%
  filter(grepl(":", term))
write_outcome(stab, A_DIR, "block_level_side_x_hemisphere.csv")

print(fam); print(stab)
cat("A17 done\n")
