# A19 — faced-scene content as moderator of the pair-level effect size.
# A14 screened DELTAS of 8 aggregated class groups (null); the synthesis named
# "specific water visibility" an unmeasured candidate driver of the site
# texture — but water IS measured in the element layer. What has never been
# tested: does the content of the scene actually FACED (natural-side water /
# grass / sky proportion) order the appraisal and physiology shifts?
# Family (declared): 6 tests = 3 natural-side content measures (water = 水体+
# 河流+海洋; grass = 草地; sky = 天空, sample means over context-recorded
# blocks) x 2 outcomes (dP_ISO, dFrontal-HbO), pair level, BH.
# Water-visible split descriptive only.

source(file.path("code", "helpers.R"))
A_DIR <- analysis_output_dir("a19_scene_content_moderators")

elements <- load_p27_elements()
samples  <- load_p27_samples()

# natural-side content per sample, context-recorded blocks only
content <- elements %>%
  filter(context_recorded %in% TRUE, side == "natural") %>%
  mutate(water = 水体 + 河流 + 海洋, grass = 草地, sky = 天空) %>%
  group_by(sample_id) %>%
  summarise(water_nat = mean(water), grass_nat = mean(grass),
            sky_nat = mean(sky), n_ctx_blocks = n(), .groups = "drop")

# outcomes: pair differences
pw <- samples %>%
  rowwise() %>%
  mutate(hbo_frontal = mean(c_across(all_of(paste0("hbo_ch", 1:4))), na.rm = TRUE)) %>%
  ungroup() %>%
  select(sample_id, participant, site, side, p_iso, hbo_frontal) %>%
  pivot_wider(names_from = side, values_from = c(p_iso, hbo_frontal)) %>%
  mutate(d_piso = p_iso_natural - p_iso_composite,
         d_hbof = hbo_frontal_natural - hbo_frontal_composite) %>%
  left_join(content, by = "sample_id")

mods <- c("water_nat", "grass_nat", "sky_nat")
outs <- c("d_piso", "d_hbof")
res <- bind_rows(lapply(mods, function(m) bind_rows(lapply(outs, function(o) {
  d <- pw %>% filter(is.finite(.data[[m]]), is.finite(.data[[o]]))
  ct <- cor.test(d[[m]], d[[o]])
  tibble(moderator = m, outcome = o, n_pairs = nrow(d),
         r = unname(ct$estimate), p = ct$p.value)
})))) %>% mutate(q = bh(p))
write_outcome(res, A_DIR, "content_moderator_screen.csv")

# descriptive: water-visible split (water > 1% of pixels on the natural side)
split <- pw %>% filter(is.finite(water_nat)) %>%
  mutate(water_visible = water_nat > 0.01) %>%
  group_by(water_visible) %>%
  summarise(n = n(), water_med = median(water_nat),
            d_piso_mean = mean(d_piso, na.rm = TRUE),
            d_piso_sd = sd(d_piso, na.rm = TRUE),
            d_hbof_mean = mean(d_hbof, na.rm = TRUE),
            d_hbof_sd = sd(d_hbof, na.rm = TRUE), .groups = "drop")
write_outcome(split, A_DIR, "water_visible_split_descriptive.csv")

# descriptive stability for the grass survivor: is natural-side grass a proxy
# for the pair's naturalness gap (A12's null moderator) or for composite-side
# grass, and does it survive adjusting for dnaturalness?
d_nat_df <- load_p27_samples() %>%
  select(sample_id, side, naturalness) %>%
  pivot_wider(names_from = side, values_from = naturalness) %>%
  mutate(d_nat = natural - composite) %>%
  select(sample_id, d_nat)
grass_comp <- elements %>%
  filter(context_recorded %in% TRUE, side == "composite") %>%
  group_by(sample_id) %>% summarise(grass_comp = mean(草地), .groups = "drop")
stab <- pw %>%
  left_join(d_nat_df, by = "sample_id") %>%
  left_join(grass_comp, by = "sample_id") %>%
  filter(is.finite(grass_nat), is.finite(d_piso))
grass_checks <- tibble(
  check = c("cor(grass_nat, d_nat)", "cor(grass_nat, grass_composite)",
            "partial cor(grass_nat, d_piso | d_nat)", "n_pairs"),
  value = c(cor(stab$grass_nat, stab$d_nat),
            cor(stab$grass_nat, stab$grass_comp, use = "complete.obs"),
            cor(residuals(lm(grass_nat ~ d_nat, data = stab)),
                residuals(lm(d_piso ~ d_nat, data = stab))),
            nrow(stab)))
write_outcome(grass_checks, A_DIR, "grass_stability_descriptive.csv")

print(res); print(split)
cat("A19 done\n")
