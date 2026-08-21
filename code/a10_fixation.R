# A10 — the eye-tracking strand.
# Families (declared): 4 between-side coupling tests and 4 within-side coupling
# tests (dfixation vs P_ISO, E_ISO, naturalness, complexity), BH within family.
# Duplicate-vector and fill-down handling = sensitivity. Side contrast on
# fixation was A02-F3 (null); repeated here only in the sensitivity table.

source(file.path("code", "helpers.R"))
A_DIR <- analysis_output_dir("a10_fixation_strand")

blocks <- load_p27_blocks()

# ---- reliability & anatomy (descriptive) -----------------------------------
fx <- fixation_rows(blocks)
anat <- fx %>% summarise(
  n_blocks = n(), n_samples = n_distinct(sample_id),
  mean_s = mean(fixation_s), median_s = median(fixation_s), sd_s = sd(fixation_s),
  share_of_60s_median = median(fixation_s) / 60,
  frac_below_5s = mean(fixation_s < 5))
write_outcome(anat, A_DIR, "fixation_anatomy.csv")

# ---- between-side coupling family ------------------------------------------
targets <- c("p_iso", "e_iso", "naturalness", "complexity")
target_rows <- function(d, y) if (y %in% CONTEXT_COLS) context_rows(d) else d

bs_pairs <- function(drop_dups = FALSE) {
  d <- fx
  if (drop_dups) d <- d %>%
      filter(!(sample_id %in% EYE_DUP_LATER &
               as.character(side) == EYE_DUP_SIDE[as.character(sample_id)]))
  fx_sm <- d %>% group_by(sample_id, participant, side) %>%
    summarise(fixation_s = mean(fixation_s), .groups = "drop")
  out <- map_dfr(targets, function(y) {
    ym <- target_rows(blocks, y) %>% group_by(sample_id, side) %>%
      summarise(yv = mean(.data[[y]], na.rm = TRUE), .groups = "drop")
    dd <- fx_sm %>% left_join(ym, by = c("sample_id", "side")) %>%
      pivot_wider(names_from = side, values_from = c(fixation_s, yv)) %>%
      mutate(dfix = fixation_s_natural - fixation_s_composite,
             dy = yv_natural - yv_composite) %>%
      filter(is.finite(dfix), is.finite(dy))
    r <- cor(dd$dfix, dd$dy)
    tt <- cor.test(dd$dfix, dd$dy)
    tibble(target = y, n = nrow(dd), r = r, p = tt$p.value)
  })
  out
}
bs <- bs_pairs(FALSE) %>% mutate(q = bh(p), set = "all")
bs_sens <- bs_pairs(TRUE) %>% mutate(q = NA_real_, set = "dups_dropped")
write_outcome(bind_rows(bs, bs_sens), A_DIR, "between_side_coupling.csv")

# ---- within-side coupling family -------------------------------------------
wb <- map_dfr(targets, function(y) {
  d <- fx %>% { if (y %in% CONTEXT_COLS) filter(., context_recorded) else . } %>%
    decompose_bw("fixation_s") %>% filter(is.finite(fixation_s_wb)) %>%
    mutate(xz = as.numeric(scale(fixation_s_wb)), yz = as.numeric(scale(.data[[y]])))
  fit <- p27_lmm(yz ~ xz, d,
                 re = c("(1 | participant)", "(1 | site)", "(1 | sample_id)", "(1 | side_unit)"))
  co <- summary(fit)$coefficients["xz", ]
  tibble(target = y, n = nrow(d), beta_std = co[["Estimate"]], se = co[["Std. Error"]],
         p = co[["Pr(>|t|)"]], singular = isSingular(fit))
}) %>% mutate(q = bh(p))
write_outcome(wb, A_DIR, "within_side_coupling.csv")

# ---- sensitivity: A02-F3 side contrast without dups / strict recorded ------
paired_fix <- function(drop_dups) {
  d <- fx
  if (drop_dups) d <- d %>%
      filter(!(sample_id %in% EYE_DUP_LATER &
               as.character(side) == EYE_DUP_SIDE[as.character(sample_id)]))
  sm <- d %>% group_by(sample_id, side) %>%
    summarise(fixation_s = mean(fixation_s), .groups = "drop") %>%
    pivot_wider(names_from = side, values_from = fixation_s) %>%
    mutate(diff = natural - composite) %>% filter(is.finite(diff))
  tt <- t.test(sm$diff)
  tibble(drop_dups = drop_dups, n_pairs = nrow(sm), mean_diff = mean(sm$diff),
         d_z = mean(sm$diff)/sd(sm$diff), p = tt$p.value)
}
write_outcome(bind_rows(paired_fix(FALSE), paired_fix(TRUE)), A_DIR, "side_contrast_sensitivity.csv")

print(bind_rows(bs, bs_sens)); print(wb)
cat("A10 done\n")
