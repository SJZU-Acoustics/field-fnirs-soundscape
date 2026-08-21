# A15 — the person level (exhaustion-audit gap-fill).
# Gate first (house rule): reliability of the per-person side effect, from the
# 13 participants measured at two sites — correlation between their two
# per-sample deltas. Reported whatever it says. Then one family of 4 moderation
# tests = {sex, age} x {dP_ISO, dfrontal-HbO}, BH. No subgroup selection.

source(file.path("code", "helpers.R"))
A_DIR <- analysis_output_dir("a15_person_level")

samples <- load_p27_samples()
parts   <- load_p27_participants()

frontal_cols <- paste0("hbo_ch", 1:4)
sm <- samples %>%
  mutate(hbo_frontal = rowMeans(across(all_of(frontal_cols)), na.rm = TRUE)) %>%
  select(sample_id, participant, side, p_iso, hbo_frontal) %>%
  pivot_wider(names_from = side, values_from = c(p_iso, hbo_frontal)) %>%
  mutate(d_piso = p_iso_natural - p_iso_composite,
         d_hbof = hbo_frontal_natural - hbo_frontal_composite)

# ---- reliability gate ------------------------------------------------------
two <- sm %>% group_by(participant) %>% filter(n() == 2) %>%
  mutate(occ = row_number()) %>% ungroup()
rel <- map_dfr(c("d_piso", "d_hbof"), function(v) {
  w <- two %>% select(participant, occ, all_of(v)) %>%
    pivot_wider(names_from = occ, values_from = all_of(v), names_prefix = "occ") %>%
    filter(is.finite(occ1), is.finite(occ2))
  ct <- suppressWarnings(cor.test(w$occ1, w$occ2))
  tibble(measure = v, n_participants = nrow(w), r = unname(ct$estimate), p = ct$p.value)
})
write_outcome(rel, A_DIR, "per_person_effect_reliability.csv")

# ---- moderation family (4 tests, BH) ---------------------------------------
# Disclosure control: the deposit bands age (5-year bands), so the two
# age-moderation tests cannot be recomputed from it. They are skipped here
# with a message; the sex tests run, and SI Table S17's age rows stand as
# printed in the article's Supplementary Information.
HAS_AGE <- "age" %in% names(parts)
if (!HAS_AGE) {
  message("A15: exact age is withheld in the deposit (5-year age bands); ",
          "skipping the two age-moderation tests (sex tests run). ",
          "SI Table S17's age rows are not regenerated.")
}
d <- sm %>% left_join(parts %>% select(participant, sex, any_of("age")), by = "participant")
traits <- if (HAS_AGE) c("sex", "age") else "sex"
mod <- map_dfr(c("d_piso", "d_hbof"), function(y) map_dfr(traits, function(x) {
  dd <- d %>% filter(is.finite(.data[[y]]))
  dd$xv <- if (x == "sex") as.numeric(dd$sex == "female") else as.numeric(scale(dd$age))
  fit <- suppressMessages(suppressWarnings(
    lmer(reformulate(c("xv", "(1 | participant)"), response = y), data = dd)))
  co <- summary(fit)$coefficients["xv", ]
  tibble(outcome = y, trait = x, n = nrow(dd), est = co[["Estimate"]],
         se = co[["Std. Error"]], p = co[["Pr(>|t|)"]], singular = isSingular(fit))
})) %>% mutate(q = bh(p))
write_outcome(mod, A_DIR, "trait_moderation.csv")

print(rel); print(mod)
cat("A15 done\n")
