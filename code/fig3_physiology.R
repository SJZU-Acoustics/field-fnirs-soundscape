# Figure 3 — Physiology: per-channel HbO side contrasts with HbT convergence,
# multivariate decoding vs permutation null, naive first-exposure replication.
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(readr); library(patchwork)})
source(file.path("code", "helpers.R"))
source(file.path("code", "style.R"))

## ---- (a) channel forest with region header rows -----------------------------
fn  <- read_csv(OUTC("analysis_06_fnirs_side_contrast", "fnirs_side_contrast.csv"), show_col_types = FALSE)
hbt <- read_csv(OUTC("analysis_13_sensitivity_robustness", "hbt_convergence.csv"), show_col_types = FALSE)

region_map <- c(ch1 = "Frontal left", ch3 = "Frontal left", ch2 = "Frontal right",
                ch4 = "Frontal right", ch5 = "Temporal left", ch7 = "Temporal left",
                ch10 = "Occipital right")
hbo <- fn %>%
  filter(chromophore == "HbO") %>%
  mutate(ch = sub("hbo_", "", var), region = region_map[ch],
         sig = ifelse(sig_q05, "q < 0.05", "not significant"),
         label = sub("ch", "channel ", ch))
hbt2 <- hbt %>% mutate(ch = sub("hbt_", "", var), label = sub("ch", "channel ", ch))
write_csv(hbo,  file.path(LOCKDIR, "fig3a_hbo.csv"))
write_csv(hbt2, file.path(LOCKDIR, "fig3a_hbt.csv"))

lev_top_to_bottom <- c("**Frontal left**",   "channel 1", "channel 3",
                       "**Frontal right**",  "channel 2", "channel 4",
                       "**Temporal left**",  "channel 5", "channel 7",
                       "**Occipital right**", "channel 10")
hbo  <- hbo  %>% mutate(label = factor(label, levels = rev(lev_top_to_bottom)))
hbt2 <- hbt2 %>% mutate(label = factor(label, levels = rev(lev_top_to_bottom)))

p_ch <- ggplot(hbo, aes(mean_diff, label)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "grey40") +
  geom_errorbar(aes(xmin = boot_lo, xmax = boot_hi), width = 0, linewidth = 0.45,
                orientation = "y") +
  geom_point(aes(fill = sig), shape = 21, size = 2.2, stroke = 0.4) +
  geom_point(data = hbt2, aes(mean_diff, label), shape = 24, size = 1.9, stroke = 0.4,
             fill = okabe_ito[1], inherit.aes = FALSE) +
  scale_fill_manual(values = c("q < 0.05" = "black", "not significant" = "white"),
                    breaks = c("q < 0.05", "not significant")) +
  scale_y_discrete(drop = FALSE) +
  scale_x_continuous(limits = c(-0.78, 0.40), breaks = c(-0.6, -0.3, 0, 0.3)) +
  labs(x = expression(Delta*"HbO, natural − composite ("*mu*"mol/L)"), y = NULL, tag = "a") +
  theme_pub() +
  theme(axis.text.y = ggtext::element_markdown(size = 8.5, colour = "black"),
        axis.ticks.y = element_blank(),
        legend.position = "inside", legend.position.inside = c(0.98, 0.98),
        legend.justification.inside = c(1, 1),
        plot.tag = element_text(size = 10, face = "bold"), plot.tag.position = c(0, 1))

## ---- (b) decoding -----------------------------------------------------------
d9  <- read_csv(OUTC("analysis_09_multivariate_pattern", "multivariate_readouts.csv"), show_col_types = FALSE)
d13 <- read_csv(OUTC("analysis_13_sensitivity_robustness", "decoding_6ch.csv"), show_col_types = FALSE)
d23 <- read_csv(OUTC("analysis_23_perception_decoding", "perception_decoding_readouts.csv"), show_col_types = FALSE)
dec <- bind_rows(
  d9  %>% filter(readout == "side_decoding_LOPO") %>%
    transmute(set = "ΔHbO\n(7 channels)", acc = statistic, null_q95 = null_q95, p = p_perm, n = n_pairs),
  d13 %>% slice(1) %>%
    transmute(set = "ΔHbO\n(6 channels)", acc = accuracy, null_q95 = null_q95, p = p_perm, n = n_pairs),
  d23 %>% filter(readout == "side_from_dPAQ_profile_LOPO") %>%
    transmute(set = "ΔPAQ\n(8 items)", acc = accuracy, null_q95 = null_q95, p = p_perm, n = n_pairs)
) %>% mutate(set = factor(set, levels = c("ΔHbO\n(7 channels)", "ΔHbO\n(6 channels)", "ΔPAQ\n(8 items)")),
             sig = ifelse(p < 0.05, "significant", "not significant"),
             lab_y = pmax(acc, null_q95) + 0.035)
write_csv(dec, file.path(LOCKDIR, "fig3b_decoding.csv"))

p_dec <- ggplot(dec, aes(set, acc)) +
  geom_col(aes(fill = sig), colour = "black", linewidth = 0.35, width = 0.62, show.legend = FALSE) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.35, colour = "grey40") +
  geom_errorbar(aes(ymin = null_q95, ymax = null_q95), width = 0.62, linewidth = 0.5,
                colour = okabe_ito[4]) +
  geom_text(aes(y = lab_y, label = sprintf("p = %.3f", p)), vjust = 0, size = 2.7,
            family = "Helvetica") +
  scale_fill_manual(values = c("significant" = "grey55", "not significant" = "white")) +
  scale_y_continuous(limits = c(0, 0.95), expand = expansion(mult = c(0, 0)),
                     breaks = c(0, 0.25, 0.5, 0.75)) +
  labs(x = NULL, y = "Leave-one-participant-out accuracy", tag = "b") +
  theme_pub() +
  theme(axis.text.x = element_text(size = 8, lineheight = 0.9),
        plot.tag = element_text(size = 10, face = "bold"), plot.tag.position = c(0, 1))

## ---- (c) naive first exposure -----------------------------------------------
blocks  <- load_p27_blocks() %>%
  mutate(side = as.character(side), sample_id = as.character(sample_id))
pairing <- load_p27_pairing()
valid   <- load_p27_validity()

vw <- valid %>%
  mutate(valid = valid %in% c(TRUE, "True", "TRUE"),
         channel = paste0("hbo_ch", channel)) %>%
  filter(channel %in% c("hbo_ch1", "hbo_ch2", "hbo_ch3", "hbo_ch4")) %>%
  select(sample_id, side, channel, valid) %>%
  pivot_wider(names_from = channel, values_from = valid, names_prefix = "v_")

first <- blocks %>%
  filter(block == 1, as.character(side) == first_side) %>%
  left_join(vw, by = c("sample_id", "side")) %>%
  mutate(fnirs_recorded = fnirs_recorded %in% c(TRUE, "True", "TRUE"))
frontal <- first %>%
  filter(fnirs_recorded) %>%
  rowwise() %>%
  mutate(hbo_frontal = mean(c(ifelse(v_hbo_ch1, hbo_ch1, NA), ifelse(v_hbo_ch2, hbo_ch2, NA),
                              ifelse(v_hbo_ch3, hbo_ch3, NA), ifelse(v_hbo_ch4, hbo_ch4, NA)),
                            na.rm = TRUE)) %>%
  ungroup()

naive <- bind_rows(
  frontal %>% transmute(sample_id, group = first_side, outcome = "hbo", value = hbo_frontal),
  first %>% transmute(sample_id, group = first_side, outcome = "mono", value = paq_monotonous),
  first %>% transmute(sample_id, group = first_side, outcome = "piso", value = p_iso)
) %>% mutate(group = ifelse(group == "natural", "Natural\nfirst", "Composite\nfirst"))
write_csv(naive, file.path(LOCKDIR, "fig3c_naive.csv"))
summ <- naive %>% group_by(outcome, group) %>%
  summarise(m = mean(value, na.rm = TRUE), se = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))),
            .groups = "drop")
print(summ)

# The two opening sides carry the figure-wide side colours (natural green,
# composite orange), on the session points and the group-mean squares alike.
COL_FIRST <- c(`Natural\nfirst` = COL_NATURAL, `Composite\nfirst` = COL_COMPOSITE)
naive_panel <- function(key, ylab, tag = NULL) {
  d <- naive %>% filter(outcome == key); s <- summ %>% filter(outcome == key)
  ggplot(d, aes(group, value)) +
    geom_jitter(aes(fill = group), width = 0.09, height = 0, shape = 21, size = 1.5,
                stroke = 0.35, alpha = 0.55, show.legend = FALSE) +
    geom_errorbar(data = s, aes(group, m, ymin = m - se, ymax = m + se),
                  width = 0.18, linewidth = 0.5, inherit.aes = FALSE) +
    geom_point(data = s, aes(group, m, fill = group), shape = 22, size = 2.4,
               inherit.aes = FALSE, show.legend = FALSE) +
    scale_fill_manual(values = COL_FIRST) +
    labs(x = NULL, y = ylab, tag = tag) +
    theme_pub(base_size = 8, axis_title_size = 9) +
    theme(plot.tag = element_text(size = 10, face = "bold"), plot.tag.position = c(0, 1))
}
# Tag sits at the left of p_n1's cell by default; after free(type = "label")
# that cell still includes the empty gutter under panel a's labels, so the
# tag would be stranded at the figure's left edge. Place it in plot
# coordinates just to the left of this panel's own y-axis title.
p_n1 <- naive_panel("hbo",  expression("Frontal HbO ("*mu*"mol/L)"), tag = "c") +
  theme(plot.tag.position = c(0.20, 1), plot.tag.location = "plot")
p_n2 <- naive_panel("mono", "Perceived monotony (1–5)")
p_n3 <- naive_panel("piso", "ISO pleasantness (−1 to 1)")

# Panel a's region/channel labels would otherwise strand panel c's y-axis titles
# at the figure's left edge. Free the whole bottom row as one patchwork
# (type = "label") so the three naive panels stay aligned with each other
# and with panel a's plot area, while their titles hug their own axes.
# type = "panel" unsticks the title by shifting the bottom-left plot area
# left. Nested free(p_n1) | p_n2 | p_n3 crashes patchwork 1.3.2 on save.
# Omit side=.
# Seeded so the jittered panel-c point clouds reproduce exactly (the working
# figure's jitter was an unseeded draw; the point positions within their
# group strips differ, nothing else).
set.seed(271)
fig <- wrap_plots(p_ch | p_dec,
                  free(p_n1 | p_n2 | p_n3, type = "label"),
                  ncol = 1, heights = c(1.25, 1))
save_fig(fig, file.path(FIGDIR, "fig3_physiology.png"), 178, 140)
cat("fig3 done\n")
