# A07 — brain-perception coupling.
# Family (declared in ../README.md): per level, 14 tests = 2 axes x 7 HbO
# channels, BH within level. HbR convergent (not run as a family).
#   between_sample: sample means of HbO (valid channels) vs P/E means, n<=40
#   between_side:   paired diffs (the causal level), n<=40 pairs
#   within_side:    block deviations, fnirs_channel_rows, LMM

source(file.path("code", "helpers.R"))
A_DIR <- analysis_output_dir("a07_brain_perception_coupling")

blocks  <- load_p27_blocks()
samples <- load_p27_samples()
AXES <- c("p_iso", "e_iso")

# ---- between_sample --------------------------------------------------------
bs_tab <- samples %>% group_by(sample_id, participant) %>%
  summarise(across(all_of(c(AXES, HBO_COLS)), ~mean(.x, na.rm = TRUE)), .groups = "drop")
bs <- map_dfr(AXES, function(y) map_dfr(CHANNELS$ch, function(ch) {
  v <- paste0("hbo_", ch)
  d <- bs_tab %>% filter(is.finite(.data[[v]])) %>%
    mutate(xz = as.numeric(scale(.data[[y]])), yz = as.numeric(scale(.data[[v]])))
  fit <- suppressMessages(suppressWarnings(lmer(yz ~ xz + (1 | participant), data = d)))
  co <- summary(fit)$coefficients["xz", ]
  tibble(level = "between_sample", axis = y, channel = ch, n = nrow(d),
         beta_std = co[["Estimate"]], se = co[["Std. Error"]], p = co[["Pr(>|t|)"]],
         singular = isSingular(fit))
}))

# ---- between_side (paired diffs) -------------------------------------------
ws <- map_dfr(AXES, function(y) map_dfr(CHANNELS$ch, function(ch) {
  v <- paste0("hbo_", ch)
  pd <- load_p27_paired(c(y, v))
  d <- pd %>% filter(is.finite(.data[[paste0(v, "__diff")]]),
                     is.finite(.data[[paste0(y, "__diff")]])) %>%
    mutate(xz = as.numeric(scale(.data[[paste0(y, "__diff")]])),
           yz = as.numeric(scale(.data[[paste0(v, "__diff")]])))
  fit <- suppressMessages(suppressWarnings(lmer(yz ~ xz + (1 | participant), data = d)))
  co <- summary(fit)$coefficients["xz", ]
  tibble(level = "between_side", axis = y, channel = ch, n = nrow(d),
         beta_std = co[["Estimate"]], se = co[["Std. Error"]], p = co[["Pr(>|t|)"]],
         singular = isSingular(fit))
}))

# ---- within_side (block deviations) ----------------------------------------
wb <- map_dfr(AXES, function(y) map_dfr(CHANNELS$ch, function(ch) {
  v <- paste0("hbo_", ch)
  d <- fnirs_channel_rows(blocks, ch) %>% decompose_bw(y) %>%
    filter(is.finite(.data[[paste0(y, "_wb")]])) %>%
    mutate(xz = as.numeric(scale(.data[[paste0(y, "_wb")]])),
           yz = as.numeric(scale(.data[[v]])))
  fit <- p27_lmm(yz ~ xz, d,
                 re = c("(1 | participant)", "(1 | site)", "(1 | sample_id)", "(1 | side_unit)"))
  co <- summary(fit)$coefficients["xz", ]
  tibble(level = "within_side", axis = y, channel = ch, n = nrow(d),
         beta_std = co[["Estimate"]], se = co[["Std. Error"]], p = co[["Pr(>|t|)"]],
         singular = isSingular(fit))
}))

res <- bind_rows(bs, ws, wb) %>% group_by(level) %>% mutate(q = bh(p)) %>%
  ungroup() %>% mutate(sig_q05 = q < 0.05) %>% arrange(level, q)
write_outcome(res, A_DIR, "coupling_screens.csv")

# convergent HbR at the causal level only (descriptive)
hbr_ws <- map_dfr(AXES, function(y) map_dfr(CHANNELS$ch, function(ch) {
  v <- paste0("hbr_", ch)
  pd <- load_p27_paired(c(y, v))
  d <- pd %>% filter(is.finite(.data[[paste0(v, "__diff")]]),
                     is.finite(.data[[paste0(y, "__diff")]]))
  tibble(axis = y, channel = ch, n = nrow(d),
         r = cor(d[[paste0(y, "__diff")]], d[[paste0(v, "__diff")]]))
}))
write_outcome(hbr_ws, A_DIR, "hbr_convergent_between_side.csv")

print(res %>% filter(p < 0.10) %>% select(level, axis, channel, n, beta_std, p, q), n = 25)
cat("A07 done\n")
