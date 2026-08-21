# A21 — naive first-exposure replication. A16 showed the effect is present from
# block 1; the strongest version of that control is exposure position 1 ONLY:
# participants who have never been exposed to either side, a fully between-
# subject contrast (natural-first vs composite-first), immune to any within-
# session carryover by construction.
# Family (declared): 3 tests — P_ISO, monotony, frontal HbO at position 1,
# y ~ first_side + (1 | participant); BH. Welch descriptive alongside.

source(file.path("code", "helpers.R"))
A_DIR <- analysis_output_dir("a21_naive_first_exposure")

blocks <- load_p27_blocks()
firsts <- blocks %>% filter(exposure_position == 1L)

# frontal HbO at position 1, valid channels only
hbof1 <- firsts %>%
  filter(fnirs_recorded %in% TRUE) %>%
  rowwise() %>%
  mutate(hbo_frontal = mean(c_across(all_of(paste0("hbo_ch", 1:4)))[
    c_across(all_of(paste0("hbo_ch", 1:4, "_valid")))], na.rm = TRUE)) %>%
  ungroup() %>% filter(is.finite(hbo_frontal)) %>%
  select(sample_id, participant, site, first_side, hbo_frontal)

run_test <- function(d, y, label) {
  f <- as.formula(paste0(y, " ~ first_side + (1 | participant)"))
  fit <- suppressMessages(suppressWarnings(lmer(f, data = d, REML = TRUE)))
  co <- summary(fit)$coefficients
  term <- grep("first_side", rownames(co), value = TRUE)
  wt <- t.test(d[[y]] ~ d$first_side)   # Welch descriptive
  gm <- d %>% group_by(first_side) %>% summarise(m = mean(.data[[y]]), .groups = "drop")
  tibble(outcome = label, n = nrow(d),
         mean_natural_first = gm$m[gm$first_side == "natural"],
         mean_composite_first = gm$m[gm$first_side == "composite"],
         welch_p = wt$p.value,
         lmm_est = co[term, "Estimate"], lmm_se = co[term, "Std. Error"],
         lmm_p = co[term, "Pr(>|t|)"], singular = isSingular(fit))
}

res <- bind_rows(
  run_test(firsts, "p_iso", "P_ISO"),
  run_test(firsts, "paq_monotonous", "monotony"),
  run_test(hbof1, "hbo_frontal", "frontal HbO")
) %>% mutate(q = bh(lmm_p))
write_outcome(res, A_DIR, "first_exposure_between_groups.csv")

print(res)
cat("A21 done\n")
