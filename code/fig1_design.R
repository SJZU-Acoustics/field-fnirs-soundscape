# Figure 1 — Design: site map, the two faced scenes, session protocol,
# fNIRS montage, acoustic constancy.
#
# Panels b and d absorb the design records drawn by the data collector for her
# Chinese draft (site plan, procedure schematic, optode montage): the montage
# geometry is digitised from the delivered figure by
# data/code/digitise_fnirs_montage.py, and the map context is OpenStreetMap
# geometry frozen by data/code/build_site_basemap.py (map data (c) OpenStreetMap
# contributors, ODbL).  Site coordinates are the WGS-84 conversion of the
# delivered GCJ-02 record (DATA_AUDIT_08.md).
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(readr); library(patchwork)})
source(file.path("code", "helpers.R"))
source(file.path("code", "style.R"))

sites   <- load_p27_sites_wgs84()
samples <- load_p27_samples()
basemap <- read_csv(file.path("data", "figure_inputs", "p27_basemap_osm.csv"), show_col_types = FALSE)
montage <- read_csv(file.path("data", "figure_inputs", "p27_montage_layout.csv"), show_col_types = FALSE)

COL_WATER <- "#C9DCE8"
COL_PARK  <- "#DCE9D5"
COL_BUILT <- "#DBDBDB"
COL_ROAD  <- "#9A9A9A"

## ---- (a) site map with its geographic context --------------------------------
n_samples_site <- samples %>% filter(side == "natural") %>% count(site) %>%
  mutate(site = as.integer(as.character(site)))
sites <- sites %>%
  left_join(n_samples_site, by = "site") %>%
  rename(longitude = longitude_wgs84, latitude = latitude_wgs84)

lat0 <- mean(sites$latitude)
asp  <- 1 / cos(lat0 * pi / 180)                      # degrees lat per degree lon
m_per_deg_lon <- 111320 * cos(lat0 * pi / 180)

# Window: the full point transect, a strip of river to the north and the road
# corridor to the south.
xlim <- c(min(sites$longitude) - 0.0021, max(sites$longitude) + 0.0021)
ylim <- c(min(sites$latitude)  - 0.0011, max(sites$latitude)  + 0.0011)

# Widen the window to the row's drawn aspect so the fixed-ratio map fills the
# figure's full width, aligning panel a with the rows beneath it (rather than
# being centred with white side margins). 3.0 ~ row-1 cell width : height at
# heights c(1.00, 0.95, 1.30) over the 178 x 181 mm canvas, net of plot margins
# (row-1 absolute height is unchanged from the earlier 178 mm square build).
target_asp <- 3.1
cur_asp <- diff(xlim) / (asp * diff(ylim))
if (cur_asp < target_asp) {
  extra <- (target_asp * asp * diff(ylim) - diff(xlim)) / 2
  xlim <- xlim + c(-extra, extra)
}

lay <- function(l) basemap %>% filter(layer == l) %>% arrange(feature_id, seq)
poly_layer <- function(l, fill, colour = NA, linewidth = 0) {
  d <- lay(l)
  if (nrow(d) == 0) return(NULL)
  geom_polygon(data = d, aes(lon, lat, group = feature_id),
               fill = fill, colour = colour, linewidth = linewidth)
}
line_layer <- function(l, colour, linewidth) {
  d <- lay(l)
  if (nrow(d) == 0) return(NULL)
  geom_path(data = d, aes(lon, lat, group = feature_id), colour = colour, linewidth = linewidth)
}

sb_len <- 200 / m_per_deg_lon                          # 200 m scale bar
sb_x0  <- xlim[1] + diff(xlim) * 0.02
sb_y   <- ylim[1] + diff(ylim) * 0.16
na_x   <- xlim[2] - diff(xlim) * 0.025
na_y0  <- ylim[1] + diff(ylim) * 0.16

p_map <- ggplot() +
  annotate("rect", xmin = xlim[1], xmax = xlim[2], ymin = ylim[1], ymax = ylim[2],
           fill = "#F4F1EC", colour = NA) +
  poly_layer("park", COL_PARK) +
  poly_layer("water", COL_WATER) +
  poly_layer("water_hole", COL_PARK) +
  poly_layer("building", COL_BUILT) +
  line_layer("footway", "white", 0.25) +
  line_layer("road", COL_ROAD, 0.5) +
  geom_point(data = sites, aes(longitude, latitude, size = factor(n)),
             shape = 21, fill = "grey20", colour = "white", stroke = 0.35) +
  ggrepel::geom_text_repel(data = sites, aes(longitude, latitude, label = site),
                           size = 2.5, family = "Helvetica", seed = 27,
                           bg.color = "white", bg.r = 0.10,
                           min.segment.length = 0.15, segment.size = 0.15,
                           max.overlaps = Inf, box.padding = 0.16, point.padding = 0.12) +
  annotate("text", x = xlim[1] + diff(xlim) * 0.42, y = ylim[2] - diff(ylim) * 0.07,
           label = "Hun River", size = 2.9, family = "Helvetica",
           fontface = "italic", colour = "#5B7C93") +
  annotate("text", x = xlim[2] - diff(xlim) * 0.02, y = ylim[1] + diff(ylim) * 0.09,
           hjust = 1, label = "city fabric", size = 2.7, family = "Helvetica", colour = "grey40") +
  annotate("segment", x = sb_x0, xend = sb_x0 + sb_len, y = sb_y, yend = sb_y, linewidth = 0.7) +
  annotate("text", x = sb_x0 + sb_len / 2, y = sb_y - diff(ylim) * 0.075, label = "200 m",
           size = 2.7, family = "Helvetica") +
  annotate("segment", x = na_x, xend = na_x, y = na_y0, yend = na_y0 + diff(ylim) * 0.14,
           arrow = arrow(length = unit(0.12, "cm"), type = "closed"), linewidth = 0.45) +
  annotate("text", x = na_x, y = na_y0 + diff(ylim) * 0.23, label = "N",
           size = 3.1, family = "Helvetica") +
  scale_size_manual(name = NULL, values = c(`1` = 1.5, `2` = 2.6),
                    labels = c("1 paired sample", "2 paired samples")) +
  coord_fixed(ratio = asp, xlim = xlim, ylim = ylim, expand = FALSE) +
  theme_void(base_family = "Helvetica") +
  theme(legend.position = "inside", legend.position.inside = c(0.003, 0.99),
        legend.justification.inside = c(0, 1),
        legend.text = element_text(size = 8),
        legend.key.spacing.y = unit(-0.03, "cm"),
        legend.margin = margin(2, 5, 2, 2),
        # opaque chip so the bridge lines crossing the river do not run
        # through the legend keys
        legend.background = element_rect(fill = alpha("white", 0.80), colour = NA),
        legend.key = element_blank(),
        plot.background = element_rect(fill = "white", colour = NA))

## ---- (b) the two faced scenes at one seat (section through the walkway) ------
GY <- 0.42                                            # ground line
circle <- function(x0, y0, r, id, n = 64) {
  tibble(th = seq(0, 2 * pi, length.out = n)) %>%
    mutate(x = x0 + r * cos(th), y = y0 + r * sin(th), id = id)
}
arc <- function(x0, y0, r, a0, a1, id, n = 40) {
  tibble(th = seq(a0, a1, length.out = n) * pi / 180) %>%
    mutate(x = x0 + r * cos(th), y = y0 + r * sin(th), id = id)
}
# A seated person in section: d = -1 faces left (natural), +1 faces right.
# Shoulders lean back onto the chair, so the two bodies read as back to back.
S <- 1.22                                             # body scale
person_limbs <- function(hx, d, id) {
  hip <- c(hx, GY + 0.52 * S); knee <- c(hx + 0.40 * S * d, GY + 0.52 * S)
  foot <- c(hx + 0.44 * S * d, GY + 0.03); sho <- c(hx - 0.06 * S * d, GY + 1.10 * S)
  hand <- c(hx + 0.32 * S * d, GY + 0.60 * S)
  tibble(id = id,
         x    = c(hip[1], hip[1],  knee[1], sho[1]),
         y    = c(hip[2], hip[2],  knee[2], sho[2]),
         xend = c(sho[1], knee[1], foot[1], hand[1]),
         yend = c(sho[2], knee[2], foot[2], hand[2]))
}
person_head <- function(hx, d) tibble(x = hx + 0.02 * S * d, y = GY + 1.32 * S)
chair_path <- function(hx, d, id) tibble(
  id = id,
  x = c(hx - 0.20 * S * d, hx - 0.20 * S * d, hx + 0.42 * S * d, hx + 0.42 * S * d),
  y = c(GY + 1.02 * S,     GY + 0.46 * S,     GY + 0.46 * S,     GY))

A_X <- -0.30; B_X <- 0.30
limbs <- bind_rows(person_limbs(A_X, -1, "A"), person_limbs(B_X, 1, "B"))
chairs <- bind_rows(chair_path(A_X, -1, "A"), chair_path(B_X, 1, "B"),
                    tibble(id = "Aleg", x = c(A_X + 0.20 * S, A_X + 0.20 * S),
                           y = c(GY + 0.46 * S, GY)),
                    tibble(id = "Bleg", x = c(B_X - 0.20 * S, B_X - 0.20 * S),
                           y = c(GY + 0.46 * S, GY)))
canopies <- bind_rows(circle(-2.62, GY + 1.62, 0.50, "t1"), circle(-2.18, GY + 1.24, 0.34, "t2"),
                      circle(-1.62, GY + 1.42, 0.42, "t3"), circle(-1.18, GY + 1.06, 0.28, "t4"))
trunks <- tibble(x = c(-2.62, -2.18, -1.62, -1.18), xend = c(-2.62, -2.18, -1.62, -1.18),
                 y = GY, yend = c(GY + 1.28, GY + 1.04, GY + 1.12, GY + 0.88))
# Buildings start well right of the seat: in the real scene people sit at a
# distance from the buildings (Prof Zhang's round, 2026-08-20).
blds <- tibble(id = c("b2", "b3"), x0 = c(1.95, 2.70), x1 = c(2.60, 3.28),
               y1 = c(GY + 2.15, GY + 1.35))
windows <- blds %>% rowwise() %>%
  do(tibble(x = seq(.$x0 + 0.15, .$x1 - 0.10, by = 0.21),
            y0 = GY + 0.25, y1 = .$y1 - 0.18)) %>% ungroup()
# The composite side carries natural elements at close range too (sparse
# vegetation near the seat, buildings behind): two small foreground trees,
# in the open ground between the seat and the set-back building row.
canopies_c <- bind_rows(circle(1.06, GY + 0.96, 0.24, "c1"),
                        circle(1.52, GY + 0.80, 0.19, "c2"))
trunks_c <- tibble(x = c(1.06, 1.52), xend = c(1.06, 1.52),
                   y = GY, yend = c(GY + 0.78, GY + 0.66))
p_scene <- ggplot() +
  annotate("rect", xmin = -3.35, xmax = 0, ymin = GY, ymax = 3.25, fill = COL_PARK, alpha = 0.55) +
  annotate("rect", xmin = 0, xmax = 3.35, ymin = GY, ymax = 3.25, fill = COL_BUILT, alpha = 0.45) +
  annotate("rect", xmin = -3.35, xmax = -2.10, ymin = GY, ymax = GY + 0.26, fill = COL_WATER) +
  geom_polygon(data = canopies, aes(x, y, group = id), fill = "#8FBF8A", colour = NA) +
  geom_segment(data = trunks, aes(x, y, xend = xend, yend = yend),
               colour = "#7A6A55", linewidth = 0.6) +
  geom_rect(data = blds, aes(xmin = x0, xmax = x1, ymin = GY, ymax = y1),
            fill = "#BFBFBF", colour = "grey35", linewidth = 0.3) +
  geom_segment(data = windows, aes(x = x, xend = x, y = y0, yend = y1),
               colour = "white", linewidth = 0.35) +
  geom_segment(data = trunks_c, aes(x, y, xend = xend, yend = yend),
               colour = "#7A6A55", linewidth = 0.5) +
  geom_polygon(data = canopies_c, aes(x, y, group = id), fill = "#8FBF8A", colour = NA) +
  annotate("segment", x = -3.35, xend = 3.35, y = GY, yend = GY, linewidth = 0.5, colour = "grey25") +
  geom_path(data = chairs, aes(x, y, group = id), colour = "grey35", linewidth = 0.5) +
  geom_segment(data = limbs, aes(x, y, xend = xend, yend = yend, colour = id),
               linewidth = 1.5, lineend = "round") +
  geom_point(data = mutate(person_head(A_X, -1), id = "A"), aes(x, y, colour = id), size = 3.1) +
  geom_point(data = mutate(person_head(B_X, 1), id = "B"), aes(x, y, colour = id), size = 3.1) +
  # A wears the cap and the eye-tracker; B wears the recorder and keeps eyes shut.
  # The cap is the montage's source colour (panel d), not A's body green — on the
  # green head a green cap disappears.
  geom_path(data = arc(A_X - 0.02 * S, GY + 1.32 * S, 0.245, 20, 210, "cap"), aes(x, y),
            colour = "#7A1B45", linewidth = 1.2) +
  annotate("segment", x = A_X - 0.02 * S - 0.24, xend = A_X - 0.02 * S - 0.24,
           y = GY + 1.32 * S - 0.09, yend = GY + 1.32 * S + 0.09,
           colour = "grey20", linewidth = 0.8) +
  geom_path(data = arc(B_X + 0.02 * S, GY + 1.32 * S, 0.235, -15, 195, "band"), aes(x, y),
            colour = "grey30", linewidth = 0.9) +
  annotate("point", x = B_X + 0.02 * S, y = GY + 1.32 * S + 0.235, size = 1.2, colour = "grey30") +
  # what each faces
  annotate("segment", x = A_X - 0.34, xend = -0.95, y = GY + 1.32 * S, yend = GY + 1.32 * S,
           arrow = arrow(length = unit(0.11, "cm"), type = "closed"),
           linewidth = 0.5, colour = COL_NATURAL) +
  annotate("segment", x = B_X + 0.34, xend = 0.95, y = GY + 1.32 * S, yend = GY + 1.32 * S,
           arrow = arrow(length = unit(0.11, "cm"), type = "closed"),
           linewidth = 0.5, colour = COL_COMPOSITE, linetype = "22") +
  # header and footer text kept clear of the drawn objects (canopy tops 2.54,
  # building top 2.57; feet at GY + 0.03) — Prof Zhang's round, 2026-08-20
  annotate("text", x = -3.28, y = 3.02, hjust = 0, size = 2.8, family = "Helvetica",
           fontface = "bold", colour = "grey20", label = "natural side") +
  annotate("text", x = 3.28, y = 3.02, hjust = 1, size = 2.8, family = "Helvetica",
           fontface = "bold", colour = "grey20", label = "composite side") +
  annotate("text", x = 0, y = 2.78, size = 2.5, family = "Helvetica", fontface = "italic",
           colour = "grey30", label = "one seat, one sound field") +
  annotate("text", x = -0.72, y = GY - 0.46, hjust = 1, size = 2.7, family = "Helvetica",
           lineheight = 0.9, label = "A  faces the scene\nfNIRS + eye-tracker") +
  annotate("text", x = 0.72, y = GY - 0.46, hjust = 0, size = 2.7, family = "Helvetica",
           lineheight = 0.9, label = "B  eyes closed\nbinaural recorder") +
  scale_colour_manual(values = c(A = COL_NATURAL, B = COL_COMPOSITE), guide = "none") +
  # clip off: A's footer line starts ~0.1 unit left of the x window at this
  # panel scale (fixed-pt text, smaller mm-per-unit) and must not lose its "f"
  coord_fixed(ratio = 1, xlim = c(-3.35, 3.35), ylim = c(-0.90, 3.25),
              expand = FALSE, clip = "off") +
  theme_void(base_family = "Helvetica") +
  theme(plot.background = element_rect(fill = "white", colour = NA))

## ---- (c) session protocol, both partners ------------------------------------
# Each overview block carries the same four segments as the magnified block
# below it — baseline, exposure, ratings, and the 120-s transition to the next
# block (drawn after blocks 1–5; no transition follows the sixth).
unit_w <- 2.85; base_w <- 0.9; expo_w <- 0.9; quest_w <- 0.6; gap_w <- 0.45
units <- tibble(idx = 1:6, side = rep(c("Natural side", "Composite side"), 3)) %>%
  mutate(x0 = (idx - 1) * unit_w)
rects <- units %>%
  rowwise() %>%
  do(tibble(idx = .$idx, side = .$side,
            part = c("baseline", "exposure", "questionnaire", "gap"),
            x0 = .$x0 + c(0, base_w, base_w + expo_w, base_w + expo_w + quest_w),
            w  = c(base_w, expo_w, quest_w, gap_w))) %>%
  ungroup() %>%
  filter(!(part == "gap" & idx == 6)) %>%
  # Transitions are genuine gaps, so they are white; the ratings segment is the
  # activity and carries the light grey (Prof Zhang's round, 2026-08-20).
  mutate(fill = case_when(part == "baseline" ~ "grey80",
                          part == "questionnaire" ~ "grey90",
                          part == "gap" ~ "white",
                          side == "Natural side" ~ COL_NATURAL,
                          TRUE ~ COL_COMPOSITE))
STRIP_Y <- c(3.15, 3.85)                                   # six-block overview
UNIT_Y  <- c(1.15, 1.95)                                   # one block, expanded
mag <- tibble(part = c("baseline", "exposure", "ratings", "gap"),
              x0 = c(0, 5.0, 10.0, 12.9), x1 = c(5.0, 10.0, 12.9, 15.6),
              y0 = UNIT_Y[1], y1 = UNIT_Y[2],
              fill = c("grey80", "white", "grey90", "white"),
              label = c("masked\nbaseline 60 s", "exposure\n60 s", "block\nratings",
                        "120 s to\nnext block"))
lanes <- tibble(y0 = c(0.60, 0.10), y1 = c(1.02, 0.52),
                x0 = c(0, 5.0), x1 = c(10.0, 10.0),
                fill = c("#DCE6F0", "#D9E7E2"),
                label = c("A   fNIRS + eye-tracking", "B   binaural recording"))
p_proto <- ggplot() +
  geom_rect(data = rects, aes(xmin = x0, xmax = x0 + w, ymin = STRIP_Y[1], ymax = STRIP_Y[2],
                              fill = I(fill)), colour = "black", linewidth = 0.3) +
  geom_text(data = units, aes(x = x0 + (base_w + expo_w + quest_w) / 2,
                              y = STRIP_Y[2] + 0.27,
                              label = ifelse(side == "Natural side", "N", "C")),
            size = 3.0, family = "Helvetica", fontface = "bold") +
  # (the six-block / role-swap sentence lives in the caption, not the panel)
  annotate("segment", x = 0, xend = 0, y = STRIP_Y[1], yend = UNIT_Y[2],
           linetype = "dashed", linewidth = 0.25, colour = "grey45") +
  annotate("segment", x = base_w + expo_w + quest_w + gap_w, xend = 15.6,
           y = STRIP_Y[1], yend = UNIT_Y[2],
           linetype = "dashed", linewidth = 0.25, colour = "grey45") +
  geom_rect(data = mag, aes(xmin = x0, xmax = x1, ymin = y0, ymax = y1, fill = I(fill)),
            colour = "black", linewidth = 0.35) +
  geom_text(data = mag, aes(x = (x0 + x1) / 2, y = (y0 + y1) / 2, label = label),
            size = 2.7, family = "Helvetica", lineheight = 0.85) +
  geom_rect(data = lanes, aes(xmin = x0, xmax = x1, ymin = y0, ymax = y1, fill = I(fill)),
            colour = NA) +
  # A's lane continues, hollow, through the block-ratings segment: no recording
  # runs there, but A is the one giving the ratings (Prof Zhang, 2026-08-20).
  annotate("rect", xmin = 10.0, xmax = 12.9, ymin = 0.60, ymax = 1.02,
           fill = NA, colour = "#8FA9C6", linewidth = 0.45) +
  geom_text(data = lanes, aes(x = x0 + 0.18, y = (y0 + y1) / 2, label = label),
            size = 2.7, family = "Helvetica", hjust = 0) +
  scale_x_continuous(limits = c(-0.25, 16.90)) +
  scale_y_continuous(limits = c(-0.05, 4.30)) +
  theme_void(base_family = "Helvetica") +
  theme(plot.background = element_rect(fill = "white", colour = NA))

## ---- (d) fNIRS montage ------------------------------------------------------
src <- montage %>% filter(kind == "source")
det <- montage %>% filter(kind == "detector") %>%
  mutate(region = factor(region, c("frontal", "temporal", "occipital")),
         status = ifelse(retained == 1, "retained", "excluded"))
lmk <- montage %>% filter(kind == "landmark")
# One channel per source-detector pair: the region's own source (temporal by side).
links <- det %>%
  mutate(src_label = case_when(region == "frontal" ~ "S1",
                               region == "temporal" & hemisphere == "left" ~ "S2",
                               region == "temporal" & hemisphere == "right" ~ "S3",
                               TRUE ~ "S4")) %>%
  left_join(src %>% select(src_label = label, sx = x, sy = y), by = "src_label")
head_circle <- tibble(th = seq(0, 2 * pi, length.out = 200)) %>%
  mutate(x = cos(th), y = sin(th))
nose <- tibble(x = c(-0.085, 0, 0.085), y = c(0.996, 1.10, 0.996))
ear <- function(side) {
  th <- seq(-70, 70, length.out = 40) * pi / 180
  tibble(x = side * (1 + 0.055 * cos(th) - 0.020), y = 0.16 * sin(th) + 0.02)
}
REGION_FILL <- c(frontal = "#7FB3AE", temporal = "#A9BDD9", occipital = "#C7C3A8")
p_mont <- ggplot() +
  geom_path(data = head_circle, aes(x, y), linewidth = 0.4, colour = "grey35") +
  geom_path(data = nose, aes(x, y), linewidth = 0.4, colour = "grey35") +
  geom_path(data = ear(-1), aes(x, y), linewidth = 0.4, colour = "grey35") +
  geom_path(data = ear(1), aes(x, y), linewidth = 0.4, colour = "grey35") +
  geom_point(data = lmk, aes(x, y), size = 0.7, colour = "grey78") +
  geom_text(data = lmk %>% filter(label %in% c("Fz", "Cz", "Pz", "T7", "T8")),
            aes(x, y, label = label), size = 2.1, family = "Helvetica",
            colour = "grey55", vjust = -0.9) +
  geom_segment(data = links, aes(x = sx, y = sy, xend = x, yend = y),
               colour = "grey45", linewidth = 0.3) +
  geom_point(data = det, aes(x, y, fill = region, alpha = status),
             shape = 21, size = 4.4, colour = "black", stroke = 0.35) +
  geom_point(data = det %>% filter(status == "excluded"), aes(x, y),
             shape = 4, size = 2.0, stroke = 0.55, colour = "grey15") +
  geom_text(data = det %>% filter(status == "retained"), aes(x, y, label = channel),
            size = 2.5, family = "Helvetica") +
  # the excluded channels keep their numbers, set radially just outside the cross
  # so the reader can see *which* channels the coverage constraint removes (6, 8, 9)
  geom_text(data = det %>% filter(status == "excluded") %>%
              mutate(r = sqrt(x^2 + y^2), lx = x * (r + 0.15) / r, ly = y * (r + 0.15) / r,
                     # channel 9's radial spot lands on the in-panel key's baseline;
                     # set it above-left of its cross instead
                     lx = ifelse(channel == 9, -0.41, lx),
                     ly = ifelse(channel == 9, -0.695, ly)),
            aes(lx, ly, label = channel), size = 2.4, family = "Helvetica",
            colour = "grey35") +
  geom_point(data = src, aes(x, y), shape = 22, size = 3.0, fill = "#7A1B45", colour = "black",
             stroke = 0.3) +
  geom_text(data = src, aes(x, y, label = sub("S", "", label)), size = 2.1,
            family = "Helvetica", colour = "white") +
  scale_fill_manual(values = REGION_FILL, guide = "none") +
  scale_alpha_manual(values = c(retained = 1, excluded = 0.30), guide = "none") +
  annotate("text", x = -1.22, y = 0.98, hjust = 0, size = 2.8, family = "Helvetica",
           colour = "grey25", label = "left") +
  annotate("text", x = 1.22, y = 0.98, hjust = 1, size = 2.8, family = "Helvetica",
           colour = "grey25", label = "right") +
  # in-panel key, drawn rather than typeset so the glyphs cannot go missing;
  # stacked into the corner the head leaves empty, so the head itself can be large
  annotate("point", x = -1.30, y = -0.72, shape = 22, size = 2.4, fill = "#7A1B45",
           colour = "black", stroke = 0.3) +
  annotate("text", x = -1.21, y = -0.72, hjust = 0, size = 2.7, family = "Helvetica",
           label = "source") +
  annotate("point", x = -1.30, y = -0.92, shape = 21, size = 3.2, fill = "grey75",
           colour = "black", stroke = 0.35) +
  annotate("text", x = -1.19, y = -0.92, hjust = 0, size = 2.7, family = "Helvetica",
           label = "detector (channel)") +
  annotate("point", x = -1.30, y = -1.12, shape = 21, size = 3.2, fill = "grey85",
           colour = "black", stroke = 0.35, alpha = 0.5) +
  annotate("point", x = -1.30, y = -1.12, shape = 4, size = 1.7, stroke = 0.5, colour = "grey15") +
  annotate("text", x = -1.19, y = -1.12, hjust = 0, size = 2.7, family = "Helvetica",
           label = "excluded for signal quality") +
  # extent matched to the cell, so the montage is drawn across panel e's whole
  # footprint (plot area plus its axis titles), not just e's plot area
  coord_fixed(xlim = c(-1.35, 1.35), ylim = c(-1.22, 1.18)) +
  theme_void(base_family = "Helvetica") +
  theme(plot.background = element_rect(fill = "white", colour = NA))

## ---- (e) acoustic constancy -------------------------------------------------
laeq_pairs <- samples %>%
  select(sample_id, side, laeq) %>%
  pivot_wider(names_from = side, values_from = laeq)
write_csv(laeq_pairs, file.path(LOCKDIR, "fig1e_laeq_pairs.csv"))
d  <- laeq_pairs$natural - laeq_pairs$composite
dz <- mean(d) / sd(d); pt <- t.test(d)$p.value
rng <- range(c(laeq_pairs$natural, laeq_pairs$composite))
p_laeq <- ggplot(laeq_pairs, aes(composite, natural)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.35, colour = "grey40") +
  geom_point(shape = 21, fill = "grey85", colour = "black", size = 1.9, stroke = 0.4) +
  annotate("text", x = rng[2] + 0.7, y = rng[1] + 1.10, hjust = 1, size = 3.1, family = "Helvetica",
           parse = TRUE,
           label = sprintf("Delta*italic(L)[plain('Aeq')] == '%+.2f dB(A)'", mean(d))) +
  annotate("text", x = rng[2] + 0.7, y = rng[1] + 0.10, hjust = 1, size = 3.1, family = "Helvetica",
           label = sprintf("italic(d)[z] == %.2f*',' ~italic(p) == %.2f", dz, pt), parse = TRUE) +
  coord_fixed(xlim = rng + c(-0.7, 0.7), ylim = rng + c(-0.7, 0.7)) +
  labs(x = expression("Composite-side "*italic(L)[plain("Aeq")]*" (dB(A))"),
       y = expression("Natural-side "*italic(L)[plain("Aeq")]*" (dB(A))")) +
  theme_pub()

# wrap_elements frees the montage from patchwork's panel-row alignment, so it is
# drawn across its whole cell — the same footprint as panel e with its axis titles.
# Row 2 raised 0.62 -> 0.74 and the canvas 178 -> 181 mm tall, so panels b and c
# gain ~17% height (Prof Zhang's round, 2026-08-20: both panels read cramped).
# 181 mm is the ceiling: at 185 mm the figure + caption overfilled the main-text
# page by 5.6 pt and the float escaped to its own page, leaving p8 blank.
# Row 2 raised again 0.74 -> 0.95 and panel b given 1.2x panel c's width, so the
# back-to-back seat drawing (the design's showpiece) is the largest element after
# the map (Prof Kang's round, 2026-08-28); canvas kept at the 181 mm ceiling and
# the map aspect re-tuned (3.0 -> 3.1) to the shorter row-1 cell.
fig <- p_map / ((p_scene | p_proto) + plot_layout(widths = c(1.2, 1))) /
  (wrap_elements(full = p_mont) | p_laeq) +
  plot_layout(heights = c(1.00, 0.95, 1.30)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 10, face = "bold", family = "Helvetica"),
        plot.tag.location = "plot", plot.tag.position = c(0, 1),
        plot.margin = margin(t = 11, r = 3, b = 2, l = 11))
save_fig(fig, file.path(FIGDIR, "fig1_design.png"), 178, 181)

write_csv(sites, file.path(LOCKDIR, "fig1a_sites.csv"))
write_csv(basemap, file.path(LOCKDIR, "fig1a_basemap_osm.csv"))
write_csv(montage, file.path(LOCKDIR, "fig1d_montage.csv"))
cat(sprintf("fig1 done; laeq d=%.3f dz=%.3f p=%.3f; %d retained channels\n",
            mean(d), dz, pt, sum(det$retained)))
