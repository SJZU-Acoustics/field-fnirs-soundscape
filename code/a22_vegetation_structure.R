# A22 — refinement of A19's grass moderator: is it grass specifically or the
# faced scene's vegetation structure? Does it hold at the site grain (where
# A12's heterogeneity lives)? Does it reach the lead item (monotony)?
# Family (declared): 4 tests, BH — (i) natural-side trees (树木) x dP_ISO;
# (ii) natural-side other vegetation (植物+花卉+棕榈树) x dP_ISO;
# (iii) site-level site-mean grass_nat x site-mean dP_ISO (n = 25);
# (iv) grass_nat x dMonotony. Joint grass+trees model descriptive.

source(file.path("code", "helpers.R"))
A_DIR <- analysis_output_dir("a22_vegetation_structure")

elements <- load_p27_elements()
samples  <- load_p27_samples()

content <- elements %>%
  filter(context_recorded %in% TRUE, side == "natural") %>%
  mutate(grass = 草地, trees = 树木, other_veg = 植物 + 花卉 + 棕榈树) %>%
  group_by(sample_id) %>%
  summarise(grass_nat = mean(grass), trees_nat = mean(trees),
            otherveg_nat = mean(other_veg), .groups = "drop")

pw <- samples %>%
  select(sample_id, site, side, p_iso, paq_monotonous) %>%
  pivot_wider(names_from = side, values_from = c(p_iso, paq_monotonous)) %>%
  mutate(d_piso = p_iso_natural - p_iso_composite,
         d_mono = paq_monotonous_natural - paq_monotonous_composite) %>%
  left_join(content, by = "sample_id")

ct_row <- function(x, y, label, d) {
  ok <- is.finite(d[[x]]) & is.finite(d[[y]])
  ct <- cor.test(d[[x]][ok], d[[y]][ok])
  tibble(test = label, n = sum(ok), r = unname(ct$estimate), p = ct$p.value)
}

site_lvl <- pw %>% group_by(site) %>%
  summarise(grass_nat = mean(grass_nat, na.rm = TRUE),
            d_piso = mean(d_piso, na.rm = TRUE), .groups = "drop") %>%
  filter(is.finite(grass_nat))

res <- bind_rows(
  ct_row("trees_nat", "d_piso", "trees_nat x dP_ISO (pair)", pw),
  ct_row("otherveg_nat", "d_piso", "otherveg_nat x dP_ISO (pair)", pw),
  ct_row("grass_nat", "d_piso", "grass_nat x dP_ISO (SITE level)", site_lvl),
  ct_row("grass_nat", "d_mono", "grass_nat x dMonotony (pair)", pw)
) %>% mutate(q = bh(p))
write_outcome(res, A_DIR, "vegetation_structure_tests.csv")

# descriptive: joint grass + trees model on dP_ISO
jfit <- lm(d_piso ~ grass_nat + trees_nat, data = pw)
joint <- broom::tidy(jfit) %>% filter(term != "(Intercept)")
write_outcome(joint, A_DIR, "joint_grass_trees_descriptive.csv")

print(res); print(joint)
cat("A22 done\n")
