# A04 — what drives P_ISO / E_ISO, level by level.
# Family (declared in ../README.md): per decomposition level, 16 tests =
# 8 predictors x 2 outcomes, BH within level.
#   between_sample:  sample means (n = 40; participant RE)
#   between_side:    within-sample side differences (n = 40 pairs; the causal level)
#   within_side:     block deviations (layer-filtered; LMM with full RE)
# Joint models = coefficient stability only.

source(file.path("code", "helpers.R"))
A_DIR <- analysis_output_dir("a04_perception_drivers")

blocks <- load_p27_blocks()
PREDICTORS <- c("naturalness", "complexity", "laeq", "l5_l95", SRC_COLS)
OUTCOMES <- c("p_iso", "e_iso")
layer_of <- function(x) {
  if (x %in% CONTEXT_COLS) "context" else if (x %in% ACOUSTIC_COLS) "acoustic" else "perception"
}
rows_for <- function(d, x) {
  switch(layer_of(x), context = context_rows(d), acoustic = acoustic_rows(d), d)
}

# ---- level 1: between-sample -----------------------------------------------
# sample means computed per layer on recorded rows; outcome from all 6 blocks
bs_tab <- blocks %>% group_by(sample_id, participant, site) %>%
  summarise(across(all_of(OUTCOMES), mean), .groups = "drop")
for (x in PREDICTORS) {
  xm <- rows_for(blocks, x) %>% group_by(sample_id) %>%
    summarise("{x}" := mean(.data[[x]], na.rm = TRUE), .groups = "drop")
  bs_tab <- left_join(bs_tab, xm, by = "sample_id")
}
bs_res <- map_dfr(OUTCOMES, function(y) map_dfr(PREDICTORS, function(x) {
  d <- bs_tab %>% filter(!is.na(.data[[x]])) %>%
    mutate(xz = as.numeric(scale(.data[[x]])), yz = as.numeric(scale(.data[[y]])))
  fit <- suppressMessages(suppressWarnings(lmer(yz ~ xz + (1 | participant), data = d)))
  co <- summary(fit)$coefficients["xz", ]
  tibble(level = "between_sample", outcome = y, predictor = x, n = nrow(d),
         beta_std = co[["Estimate"]], se = co[["Std. Error"]], p = co[["Pr(>|t|)"]],
         singular = isSingular(fit))
}))

# ---- level 2: within-sample-between-side (the causal level) ----------------
side_means <- blocks %>% group_by(sample_id, participant, side) %>%
  summarise(across(all_of(OUTCOMES), mean), .groups = "drop")
for (x in PREDICTORS) {
  xm <- rows_for(blocks, x) %>% group_by(sample_id, side) %>%
    summarise("{x}" := mean(.data[[x]], na.rm = TRUE), .groups = "drop")
  side_means <- left_join(side_means, xm, by = c("sample_id", "side"))
}
diffs <- side_means %>% pivot_wider(names_from = side,
    values_from = all_of(c(OUTCOMES, PREDICTORS)), names_glue = "{.value}__{side}")
for (v in c(OUTCOMES, PREDICTORS)) diffs[[paste0("d_", v)]] <-
  diffs[[paste0(v, "__natural")]] - diffs[[paste0(v, "__composite")]]
ws_res <- map_dfr(OUTCOMES, function(y) map_dfr(PREDICTORS, function(x) {
  d <- diffs %>% filter(!is.na(.data[[paste0("d_", x)]]), !is.na(.data[[paste0("d_", y)]])) %>%
    mutate(xz = as.numeric(scale(.data[[paste0("d_", x)]])),
           yz = as.numeric(scale(.data[[paste0("d_", y)]])))
  fit <- suppressMessages(suppressWarnings(lmer(yz ~ xz + (1 | participant), data = d)))
  co <- summary(fit)$coefficients["xz", ]
  tibble(level = "between_side", outcome = y, predictor = x, n = nrow(d),
         beta_std = co[["Estimate"]], se = co[["Std. Error"]], p = co[["Pr(>|t|)"]],
         singular = isSingular(fit))
}))

# ---- level 3: within-side block deviations ---------------------------------
wb_res <- map_dfr(OUTCOMES, function(y) map_dfr(PREDICTORS, function(x) {
  d <- rows_for(blocks, x) %>% decompose_bw(x) %>%
    filter(!is.na(.data[[paste0(x, "_wb")]])) %>%
    mutate(xz = as.numeric(scale(.data[[paste0(x, "_wb")]])),
           yz = as.numeric(scale(.data[[y]])))
  fit <- p27_lmm(yz ~ xz, d,
                 re = c("(1 | participant)", "(1 | site)", "(1 | sample_id)", "(1 | side_unit)"))
  co <- summary(fit)$coefficients["xz", ]
  tibble(level = "within_side", outcome = y, predictor = x, n = nrow(d),
         beta_std = co[["Estimate"]], se = co[["Std. Error"]], p = co[["Pr(>|t|)"]],
         singular = isSingular(fit))
}))

res <- bind_rows(bs_res, ws_res, wb_res) %>%
  group_by(level) %>% mutate(q = bh(p)) %>% ungroup() %>%
  mutate(sig_q05 = q < 0.05) %>% arrange(level, q)
write_outcome(res, A_DIR, "driver_screens.csv")

# ---- joint model at the causal level (stability, no family) ----------------
joint <- diffs %>% filter(!is.na(d_naturalness), !is.na(d_laeq)) %>%
  mutate(across(c(d_p_iso, d_naturalness, d_complexity, d_laeq), ~as.numeric(scale(.x))))
jfit <- suppressMessages(suppressWarnings(
  lmer(d_p_iso ~ d_naturalness + d_complexity + d_laeq + (1 | participant), data = joint)))
jco <- as_tibble(summary(jfit)$coefficients, rownames = "term") %>%
  mutate(singular = isSingular(jfit))
write_outcome(jco, A_DIR, "joint_model_between_side.csv")

print(res %>% filter(q < 0.10) %>% select(level, outcome, predictor, n, beta_std, p, q), n = 30)
print(jco)
cat("A04 done\n")
