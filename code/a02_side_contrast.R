# A02 — the side contrast, formalised.
# Within-sample natural-minus-composite differences over the 40 pairs for every
# outcome, in the declared families:
#   F1 perception: 14 tests (8 PAQ + P_ISO + E_ISO + 4 sources), BH
#   F2 physical/context: 7 tests (4 acoustic + 3 context), BH
#   F3 fixation: 1 test, native p
# Estimator: side_contrast() = paired mean difference with participant random
# intercept, d_z, and 1,000-resample participant-cluster bootstrap CI.

source(file.path("code", "helpers.R"))
A_DIR <- analysis_output_dir("a02_side_contrast")

f1_vars <- c(PAQ_COLS, "p_iso", "e_iso", SRC_COLS)
f2_vars <- c(ACOUSTIC_COLS, CONTEXT_COLS)

f1 <- map_dfr(f1_vars, side_contrast) %>% mutate(family = "F1_perception", q = bh(p_t))
f2 <- map_dfr(f2_vars, side_contrast) %>% mutate(family = "F2_physical_context", q = bh(p_t))
f3 <- side_contrast("fixation_s") %>% mutate(family = "F3_fixation", q = p_t)

res <- bind_rows(f1, f2, f3) %>%
  mutate(sig_q05 = q < 0.05) %>%
  arrange(family, q)
write_outcome(res, A_DIR, "side_contrast_families.csv")

print(res %>% select(family, var, n_pairs, mean_diff, d_z, p_t, q, sig_q05, lmm_singular), n = 40)
cat("A02 done\n")
