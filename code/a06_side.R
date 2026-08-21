# A06 — the side contrast on cortical haemodynamics.
# Family (declared): 7 tests, HbO per channel, sample-level paired differences,
# BH. HbR convergent. Block-level LMM with exposure_position = stability.

source(file.path("code", "helpers.R"))
A_DIR <- analysis_output_dir("a06_fnirs_side_contrast")

blocks <- load_p27_blocks()

hbo <- map_dfr(paste0("hbo_", CHANNELS$ch), side_contrast) %>%
  mutate(chromophore = "HbO", q = bh(p_t), sig_q05 = q < 0.05)
hbr <- map_dfr(paste0("hbr_", CHANNELS$ch), side_contrast) %>%
  mutate(chromophore = "HbR", q = NA_real_, sig_q05 = NA)
res <- bind_rows(hbo, hbr)
write_outcome(res, A_DIR, "fnirs_side_contrast.csv")

# stability: block-level LMM side effect with exposure_position covariate
stab <- map_dfr(CHANNELS$ch, function(ch) {
  v <- paste0("hbo_", ch)
  d <- fnirs_channel_rows(blocks, ch) %>%
    mutate(yz = as.numeric(scale(.data[[v]])), pos = as.numeric(scale(exposure_position)))
  fit <- p27_lmm(yz ~ side + pos, d,
                 re = c("(1 | participant)", "(1 | site)", "(1 | sample_id)"))
  co <- as_tibble(summary(fit)$coefficients, rownames = "term") %>%
    filter(term != "(Intercept)") %>%
    mutate(channel = ch, singular = isSingular(fit))
  co
})
write_outcome(stab, A_DIR, "block_level_stability.csv")

print(res %>% select(chromophore, var, n_pairs, mean_diff, d_z, p_t, q))
print(stab %>% filter(term == "sidenatural") %>%
        select(channel, Estimate, `Std. Error`, `Pr(>|t|)`, singular))
cat("A06 done\n")
