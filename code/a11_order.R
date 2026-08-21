# A11 — order and time structure.
# Families (declared): (i) 4 side x first_side interactions (p_iso,
# paq_monotonous, hbo_ch1, hbo_ch4), BH; (ii) 9 exposure_position trends
# (p_iso, e_iso, 7 HbO channels), BH. Session date/time descriptive.

source(file.path("code", "helpers.R"))
A_DIR <- analysis_output_dir("a11_order_and_time")

blocks <- load_p27_blocks()

# ---- (i) does starting side moderate the headline effects? -----------------
int_test <- function(v, fnirs = FALSE) {
  d <- if (fnirs) fnirs_channel_rows(blocks, sub("hbo_", "", v)) else blocks
  d <- d %>% mutate(yz = as.numeric(scale(.data[[v]])),
                    fs = factor(first_side, levels = c("composite", "natural")))
  fit <- p27_lmm(yz ~ side * fs, d,
                 re = c("(1 | participant)", "(1 | site)", "(1 | sample_id)"))
  co <- summary(fit)$coefficients
  tibble(var = v, n = nrow(d),
         side_main = co["sidenatural", "Estimate"],
         inter_est = co["sidenatural:fsnatural", "Estimate"],
         inter_se = co["sidenatural:fsnatural", "Std. Error"],
         inter_p = co["sidenatural:fsnatural", "Pr(>|t|)"],
         singular = isSingular(fit))
}
fam1 <- bind_rows(int_test("p_iso"), int_test("paq_monotonous"),
                  int_test("hbo_ch1", TRUE), int_test("hbo_ch4", TRUE)) %>%
  mutate(q = bh(inter_p))
write_outcome(fam1, A_DIR, "first_side_moderation.csv")

# ---- (ii) exposure-position trends -----------------------------------------
trend_test <- function(v, fnirs = FALSE) {
  d <- if (fnirs) fnirs_channel_rows(blocks, sub("hbo_", "", v)) else blocks
  d <- d %>% mutate(yz = as.numeric(scale(.data[[v]])),
                    pos = as.numeric(scale(exposure_position)))
  fit <- p27_lmm(yz ~ side + pos, d,
                 re = c("(1 | participant)", "(1 | site)", "(1 | sample_id)"))
  co <- summary(fit)$coefficients["pos", ]
  tibble(var = v, n = nrow(d), pos_beta = co[["Estimate"]], se = co[["Std. Error"]],
         p = co[["Pr(>|t|)"]], singular = isSingular(fit))
}
fam2 <- bind_rows(map_dfr(c("p_iso", "e_iso"), trend_test),
                  map_dfr(HBO_COLS, ~trend_test(.x, TRUE))) %>%
  mutate(q = bh(p))
write_outcome(fam2, A_DIR, "exposure_position_trends.csv")

# ---- descriptive: session date / time of day -------------------------------
# The deposit omits session clock-times (disclosure control); session_time is
# constant within sample, so dropping it from the distinct() changes nothing.
sess <- blocks %>% distinct(sample_id, participant, session_date) %>%
  left_join(load_p27_paired("p_iso") %>% select(sample_id, p_iso__diff), by = "sample_id")
by_date <- sess %>% group_by(session_date) %>%
  summarise(n = n(), mean_dpiso = mean(p_iso__diff, na.rm = TRUE))
write_outcome(by_date, A_DIR, "by_session_date.csv")

print(fam1); print(fam2)
cat("A11 done\n")
