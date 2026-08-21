# Figure 4 — Generality and boundary: per-site effects on the park map and as
# aligned forests, the dose-response null, and the lawn moderator.
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(readr); library(patchwork)})
source(file.path("code", "helpers.R"))
source(file.path("code", "style.R"))

per_site <- read_csv(OUTC("analysis_12_site_generality", "per_site_effects.csv"), show_col_types = FALSE)
sites    <- load_p27_sites()  # read for parity with the working pipeline; not used downstream
samples  <- load_p27_samples()
elems    <- load_p27_elements()

## ---- (a) aligned per-site forests ------------------------------------------
ord <- per_site %>% arrange(d_piso) %>% pull(site)
forest <- per_site %>%
  mutate(site_f = factor(site, levels = ord)) %>%
  select(site_f, n_samples, d_piso, d_hbof)
write_csv(forest, file.path(LOCKDIR, "fig4a_forest.csv"))

forest_panel <- function(var, xlab, tag = NULL, show_y = TRUE) {
  p <- ggplot(forest, aes(.data[[var]], site_f)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "grey40") +
    geom_segment(aes(x = 0, xend = .data[[var]], yend = site_f), linewidth = 0.35,
                 colour = "grey60") +
    geom_point(aes(fill = .data[[var]] > 0), shape = 21, size = 1.7, stroke = 0.35,
               show.legend = FALSE) +
    scale_fill_manual(values = c(`TRUE` = COL_NATURAL, `FALSE` = COL_COMPOSITE)) +
    labs(x = xlab, y = if (show_y) "Site (ordered by site-mean Δ ISO pleasantness)" else NULL,
         tag = tag) +
    theme_pub(base_size = 8, axis_title_size = 9) +
    theme(plot.tag = element_text(size = 10, face = "bold"), plot.tag.position = c(0, 1))
  if (!show_y) p <- p + theme(axis.text.y = element_blank())
  p
}
p_forest <- forest_panel("d_piso", "Site-mean Δ ISO pleasantness", tag = "a") |
  forest_panel("d_hbof", expression("Site-mean "*Delta*" frontal HbO ("*mu*"mol/L)"),
               show_y = FALSE)

## ---- (b) dose-response null -------------------------------------------------
pairs <- samples %>%
  select(sample_id, side, p_iso, naturalness) %>%
  pivot_wider(names_from = side, values_from = c(p_iso, naturalness)) %>%
  mutate(d_piso = p_iso_natural - p_iso_composite,
         d_nat  = naturalness_natural - naturalness_composite)
write_csv(pairs, file.path(LOCKDIR, "fig4b_dose_pairs.csv"))

p_dose <- ggplot(pairs, aes(d_nat, d_piso)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "grey40") +
  geom_smooth(method = "lm", se = TRUE, colour = okabe_ito[1], fill = "grey85", linewidth = 0.6) +
  geom_point(shape = 21, fill = "grey85", size = 1.8, stroke = 0.4) +
  annotate("text", x = -0.12, y = 0.52, hjust = 0, size = 2.8, family = "Helvetica",
           label = "β = +0.11, p = 0.50") +
  labs(x = "Δ Viewed naturalness (proportion of scene)", y = "Δ ISO pleasantness", tag = "b") +
  theme_pub() +
  theme(plot.tag = element_text(size = 10, face = "bold"), plot.tag.position = c(0, 1))

## ---- (c) lawn moderator at the site grain -----------------------------------
grass <- elems %>%
  filter(side == "natural", context_recorded %in% c(TRUE, "True", "TRUE")) %>%
  group_by(sample_id) %>%
  summarise(grass_nat = mean(`草地`), .groups = "drop") %>%
  left_join(samples %>% filter(side == "natural") %>% select(sample_id, site), by = "sample_id") %>%
  left_join(pairs %>% select(sample_id, d_piso), by = "sample_id") %>%
  group_by(site) %>%
  summarise(grass_nat = mean(grass_nat), d_piso = mean(d_piso), .groups = "drop")
write_csv(grass, file.path(LOCKDIR, "fig4c_grass_site.csv"))
ct <- cor.test(grass$grass_nat, grass$d_piso)
cat(sprintf("site grass r = %.3f p = %.4f (A22: r = 0.525, p = 0.0070)\n", ct$estimate, ct$p.value))

p_grass <- ggplot(grass, aes(grass_nat, d_piso)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "grey40") +
  geom_smooth(method = "lm", se = TRUE, colour = okabe_ito[1], fill = "grey85", linewidth = 0.6) +
  geom_point(shape = 21, fill = COL_NATURAL, size = 2.0, stroke = 0.4) +
  annotate("text", x = 0.02, y = 0.45, hjust = 0, size = 2.8, family = "Helvetica",
           label = "r = +0.53, q = 0.028") +
  labs(x = "Faced lawn proportion, natural side", y = "Site-mean Δ ISO pleasantness", tag = "c") +
  theme_pub() +
  theme(plot.tag = element_text(size = 10, face = "bold"), plot.tag.position = c(0, 1))

fig <- p_forest / (p_dose | p_grass) + plot_layout(heights = c(1.15, 1))
save_fig(fig, file.path(FIGDIR, "fig4_generality.png"), 178, 150)
cat("fig4 done\n")
