# A12 — generality of the side effect across the 25 sites.
# Families (declared): (i) 4 site-level dose-response tests (site-mean dP_ISO
# and dHbO_frontal vs site dnaturalness; and vs distance_to_main_road), BH;
# (ii) 2 site-heterogeneity permutation tests (P_ISO, frontal HbO): statistic =
# between-site dispersion of site-mean pair differences, null = site labels
# shuffled across samples (1,000 resamples), assumption-light per house rule 7.

source(file.path("code", "helpers.R"))
A_DIR <- analysis_output_dir("a12_site_generality")

samples <- load_p27_samples()
sites   <- load_p27_sites()

# frontal HbO mean at sample-side level (A05: within-frontal coherence 0.49)
frontal_cols <- paste0("hbo_ch", 1:4)
samples <- samples %>%
  mutate(hbo_frontal = rowMeans(across(all_of(frontal_cols)), na.rm = TRUE))

# pair diffs
pw <- samples %>% select(sample_id, participant, site, side, p_iso, naturalness, hbo_frontal) %>%
  pivot_wider(names_from = side, values_from = c(p_iso, naturalness, hbo_frontal)) %>%
  mutate(d_piso = p_iso_natural - p_iso_composite,
         d_nat  = naturalness_natural - naturalness_composite,
         d_hbof = hbo_frontal_natural - hbo_frontal_composite)

# per-site summary (sites carry 1-2 samples)
site_sum <- pw %>% group_by(site) %>%
  summarise(n_samples = n(), d_piso = mean(d_piso, na.rm = TRUE),
            d_nat = mean(d_nat, na.rm = TRUE), d_hbof = mean(d_hbof, na.rm = TRUE),
            .groups = "drop") %>%
  left_join(sites %>% mutate(site = factor(site)) %>%
              select(site, distance_to_main_road_m), by = "site")
write_outcome(site_sum, A_DIR, "per_site_effects.csv")

# ---- family (i): dose-response over sites ----------------------------------
tests <- list(
  c("d_piso", "d_nat"), c("d_hbof", "d_nat"),
  c("d_piso", "distance_to_main_road_m"), c("d_hbof", "distance_to_main_road_m"))
fam1 <- map_dfr(tests, function(tp) {
  d <- site_sum %>% filter(is.finite(.data[[tp[1]]]), is.finite(.data[[tp[2]]]))
  ct <- cor.test(d[[tp[1]]], d[[tp[2]]])
  tibble(outcome = tp[1], predictor = tp[2], n_sites = nrow(d),
         r = unname(ct$estimate), p = ct$p.value)
}) %>% mutate(q = bh(p))
write_outcome(fam1, A_DIR, "site_dose_response.csv")

# ---- family (ii): heterogeneity permutation --------------------------------
het_perm <- function(dv) {
  d <- pw %>% filter(is.finite(.data[[dv]]))
  stat <- function(site_labels) {
    tibble(s = site_labels, x = d[[dv]]) %>% group_by(s) %>%
      summarise(m = mean(x), n = n(), .groups = "drop") %>%
      summarise(sum(n * (m - mean(d[[dv]]))^2)) %>% pull()
  }
  obs <- stat(d$site)
  set.seed(271)
  nulls <- replicate(BOOT_N, stat(sample(d$site)))
  tibble(outcome = dv, n_pairs = nrow(d), n_sites = n_distinct(d$site),
         stat_obs = obs, null_q95 = quantile(nulls, .95), p_perm = mean(nulls >= obs))
}
fam2 <- bind_rows(het_perm("d_piso"), het_perm("d_hbof"))
write_outcome(fam2, A_DIR, "site_heterogeneity_permutation.csv")

# descriptive: sign consistency across sites
signs <- site_sum %>% summarise(
  piso_pos = sum(d_piso > 0, na.rm = TRUE), piso_n = sum(is.finite(d_piso)),
  hbof_neg = sum(d_hbof < 0, na.rm = TRUE), hbof_n = sum(is.finite(d_hbof)))
write_outcome(signs, A_DIR, "sign_consistency.csv")

print(fam1); print(fam2); print(signs)
cat("A12 done\n")
