# Figure 2 — Perception: side contrasts (items + sources; ISO axes) and the
# circumplex shift profile with its fitted pleasantness-axis displacement.
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(readr); library(patchwork)})
source(file.path("code", "style.R"))

fam <- read_csv(OUTC("analysis_02_side_contrast", "side_contrast_families.csv"), show_col_types = FALSE)

# paq_engaging is the data-column name; the manuscript item label is the ISO
# attribute "vibrant" (collector's correction, 2026-08-19).
lab_map <- c(paq_pleasant = "pleasant", paq_eventful = "eventful", paq_engaging = "vibrant",
             paq_chaotic = "chaotic", paq_annoying = "annoying", paq_monotonous = "monotonous",
             paq_uneventful = "uneventful", paq_calm = "calm",
             src_traffic = "traffic sounds", src_human = "human sounds",
             src_natural = "natural sounds", src_other = "other sounds")

items <- fam %>%
  filter(var %in% names(lab_map)) %>%
  mutate(label = unname(lab_map[var]),
         group = ifelse(grepl("^src_", var), "Source salience", "Affective quality"),
         sig = ifelse(sig_q05, "q < 0.05", "not significant"))
write_csv(items, file.path(LOCKDIR, "fig2a_items.csv"))

# House forest recipe: group labels as bold header pseudo-rows inside the y axis.
aq  <- items %>% filter(group == "Affective quality") %>% arrange(desc(mean_diff))
src <- items %>% filter(group == "Source salience") %>% arrange(desc(mean_diff))
lev_top_to_bottom <- c("**Affective quality**", aq$label, " ", "**Source salience**", src$label)
items <- items %>% mutate(label = factor(label, levels = rev(lev_top_to_bottom)))

p_items <- ggplot(items, aes(mean_diff, label)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "grey40") +
  geom_errorbar(aes(xmin = boot_lo, xmax = boot_hi), width = 0, linewidth = 0.45,
                orientation = "y") +
  geom_point(aes(fill = sig), shape = 21, size = 2.1, stroke = 0.4, show.legend = FALSE) +
  scale_fill_manual(values = c("q < 0.05" = "black", "not significant" = "white")) +
  scale_y_discrete(drop = FALSE) +
  labs(x = "Difference (scale steps)", y = NULL, tag = "a") +
  theme_pub() +
  theme(axis.text.y = ggtext::element_markdown(size = 9, colour = "black"),
        axis.ticks.y = element_blank(),
        plot.tag = element_text(size = 10, face = "bold"), plot.tag.position = c(0, 1))

axes <- fam %>%
  filter(var %in% c("p_iso", "e_iso")) %>%
  mutate(label = ifelse(var == "p_iso", "ISO pleasantness", "ISO eventfulness"),
         sig = ifelse(sig_q05, "q < 0.05", "not significant")) %>%
  mutate(label = factor(label, levels = c("ISO eventfulness", "ISO pleasantness")))
write_csv(axes, file.path(LOCKDIR, "fig2b_axes.csv"))

p_axes <- ggplot(axes, aes(mean_diff, label)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "grey40") +
  geom_errorbar(aes(xmin = boot_lo, xmax = boot_hi), width = 0, linewidth = 0.45,
                orientation = "y") +
  geom_point(aes(fill = sig), shape = 21, size = 2.3, stroke = 0.4, show.legend = FALSE) +
  scale_fill_manual(values = c("q < 0.05" = "black", "not significant" = "white")) +
  scale_x_continuous(limits = c(-0.06, 0.14), breaks = c(-0.05, 0, 0.05, 0.10)) +
  labs(x = "Difference (ISO units)", y = NULL, tag = "b") +
  theme_pub() +
  theme(plot.tag = element_text(size = 10, face = "bold"), plot.tag.position = c(0, 1))

## ---- (c) circumplex shift profile ------------------------------------------
prof <- read_csv(OUTC("analysis_03_perception_structure", "item_shift_profile.csv"), show_col_types = FALSE)
fit  <- read_csv(OUTC("analysis_03_perception_structure", "shift_profile_fit.csv"), show_col_types = FALSE)
a_cos <- fit$est[fit$term == "cos_t"]; a_sin <- fit$est[fit$term == "sin_t"]
curve_df <- tibble(angle = seq(0, 360, 2),
                   fitted = a_cos * cos(angle * pi / 180) + a_sin * sin(angle * pi / 180))
prof <- prof %>% mutate(label = sub("^paq_", "", item),
                        label = ifelse(label == "engaging", "vibrant", label))
write_csv(prof, file.path(LOCKDIR, "fig2c_profile.csv"))

# The eight item labels are placed by explicit per-item nudges into pockets
# clear of the fitted curve and every error bar (Prof Zhang's round,
# 2026-08-20); obstacle rows (empty labels repel but do not draw) trace the
# curve and the bar tips so ggrepel cannot drift a label back onto either.
nudges <- tibble::tribble(
  ~item,             ~nx, ~ny,
  "paq_pleasant",     20,  0.155,
  "paq_eventful",     40,  0.070,
  "paq_engaging",     57, -0.160,
  "paq_chaotic",     -36, -0.075,
  "paq_annoying",     12,  0.115,
  "paq_monotonous",   36, -0.048,
  "paq_uneventful",   33, -0.015,
  "paq_calm",         31, -0.165)
repel_df <- bind_rows(
  prof %>% select(item, angle_deg, mean_shift, label) %>%
    left_join(nudges, by = "item") %>% select(-item),
  curve_df %>% filter(angle %% 10 == 0) %>%
    transmute(angle_deg = angle, mean_shift = fitted, label = "", nx = 0, ny = 0),
  prof %>% transmute(angle_deg, mean_shift = mean_shift - se, label = "", nx = 0, ny = 0),
  prof %>% transmute(angle_deg, mean_shift = mean_shift + se, label = "", nx = 0, ny = 0))
p_prof <- ggplot(prof, aes(angle_deg, mean_shift)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "grey40") +
  geom_line(data = curve_df, aes(angle, fitted), colour = okabe_ito[1], linewidth = 0.6) +
  geom_errorbar(aes(ymin = mean_shift - se, ymax = mean_shift + se), width = 0, linewidth = 0.45) +
  geom_point(shape = 21, fill = "grey85", size = 2.1, stroke = 0.4) +
  ggrepel::geom_text_repel(data = repel_df, aes(label = label),
                           nudge_x = repel_df$nx, nudge_y = repel_df$ny,
                           size = 2.6, family = "Helvetica", seed = 3,
                           min.segment.length = 0.2, segment.size = 0.2,
                           segment.colour = "grey55",
                           box.padding = 0.25, point.padding = 0.20, max.overlaps = Inf) +
  scale_x_continuous(breaks = c(0, 90, 180, 270, 360)) +
  scale_y_continuous(expand = expansion(mult = 0.12)) +
  labs(x = "Item angle on the circumplex (°)",
       y = "Item shift (scale steps)", tag = "c") +
  theme_pub() +
  theme(plot.margin = margin(5.5, 12, 5.5, 5.5),
        plot.tag = element_text(size = 10, face = "bold"), plot.tag.position = c(0, 1))

# Panel b's categorical y labels ("ISO pleasantness") would otherwise strand
# panel c's y-axis title far from its numeric ticks. type = "label" unsticks
# the title while keeping b and c plot-area alignment. Nested
# p_items | (p_axes / free(p_prof)) and wrap_plots(..., design = "AB\nAC")
# either crash patchwork 1.3.2 on save or ignore free(); isolate the right
# stack and wrap_elements() it so the left forest is not nested with free().
# Omit side=; do not nest two free() plots.
right <- p_axes / free(p_prof, type = "label") + plot_layout(heights = c(1, 2.6))
# Panel c widened again (right column 1 : 1.42, Prof Zhang's round 2026-08-20)
# so the circumplex item labels have room clear of the fitted line; the left
# forest keeps enough width for its bold header rows.
fig <- (p_items | wrap_elements(full = right, clip = FALSE)) +
  plot_layout(widths = c(1, 1.42))
save_fig(fig, file.path(FIGDIR, "fig2_perception.png"), 178, 110)
cat("fig2 done\n")
