# A08 — brain-context/acoustic coupling.
# Family (declared): per level, 21 tests = 3 predictors (naturalness,
# complexity, laeq) x 7 HbO channels, BH within level.

source(file.path("code", "helpers.R"))
A_DIR <- analysis_output_dir("a08_brain_context_acoustic_coupling")

blocks <- load_p27_blocks()
PREDS <- c("naturalness", "complexity", "laeq")
pred_layer_rows <- function(d, x) if (x == "laeq") acoustic_rows(d) else context_rows(d)

# ---- between_sample --------------------------------------------------------
samples <- load_p27_samples()
bs_tab <- samples %>% group_by(sample_id, participant) %>%
  summarise(across(all_of(c(PREDS, HBO_COLS)), ~mean(.x, na.rm = TRUE)), .groups = "drop")
bs <- map_dfr(PREDS, function(x) map_dfr(CHANNELS$ch, function(ch) {
  v <- paste0("hbo_", ch)
  d <- bs_tab %>% filter(is.finite(.data[[v]]), is.finite(.data[[x]])) %>%
    mutate(xz = as.numeric(scale(.data[[x]])), yz = as.numeric(scale(.data[[v]])))
  fit <- suppressMessages(suppressWarnings(lmer(yz ~ xz + (1 | participant), data = d)))
  co <- summary(fit)$coefficients["xz", ]
  tibble(level = "between_sample", predictor = x, channel = ch, n = nrow(d),
         beta_std = co[["Estimate"]], se = co[["Std. Error"]], p = co[["Pr(>|t|)"]],
         singular = isSingular(fit))
}))

# ---- between_side ----------------------------------------------------------
side_means <- blocks %>% group_by(sample_id, participant, side) %>%
  summarise(.groups = "drop") %>% distinct()
# per-layer side means
sm_pred <- map(PREDS, function(x) {
  pred_layer_rows(blocks, x) %>% group_by(sample_id, side) %>%
    summarise("{x}" := mean(.data[[x]], na.rm = TRUE), .groups = "drop")
}) %>% reduce(full_join, by = c("sample_id", "side"))
sm_hbo <- samples %>% select(sample_id, participant, side, all_of(HBO_COLS))
sm <- left_join(sm_hbo, sm_pred, by = c("sample_id", "side"))
wide <- sm %>% pivot_wider(names_from = side, values_from = all_of(c(PREDS, HBO_COLS)),
                           names_glue = "{.value}__{side}")
for (v in c(PREDS, HBO_COLS)) wide[[paste0("d_", v)]] <-
  wide[[paste0(v, "__natural")]] - wide[[paste0(v, "__composite")]]
ws <- map_dfr(PREDS, function(x) map_dfr(CHANNELS$ch, function(ch) {
  v <- paste0("hbo_", ch)
  d <- wide %>% filter(is.finite(.data[[paste0("d_", v)]]), is.finite(.data[[paste0("d_", x)]])) %>%
    mutate(xz = as.numeric(scale(.data[[paste0("d_", x)]])),
           yz = as.numeric(scale(.data[[paste0("d_", v)]])))
  fit <- suppressMessages(suppressWarnings(lmer(yz ~ xz + (1 | participant), data = d)))
  co <- summary(fit)$coefficients["xz", ]
  tibble(level = "between_side", predictor = x, channel = ch, n = nrow(d),
         beta_std = co[["Estimate"]], se = co[["Std. Error"]], p = co[["Pr(>|t|)"]],
         singular = isSingular(fit))
}))

# ---- within_side -----------------------------------------------------------
wb <- map_dfr(PREDS, function(x) map_dfr(CHANNELS$ch, function(ch) {
  v <- paste0("hbo_", ch)
  d <- fnirs_channel_rows(blocks, ch) %>%
    { if (x == "laeq") filter(., acoustic_recorded) else filter(., context_recorded) } %>%
    decompose_bw(x) %>% filter(is.finite(.data[[paste0(x, "_wb")]])) %>%
    mutate(xz = as.numeric(scale(.data[[paste0(x, "_wb")]])),
           yz = as.numeric(scale(.data[[v]])))
  fit <- p27_lmm(yz ~ xz, d,
                 re = c("(1 | participant)", "(1 | site)", "(1 | sample_id)", "(1 | side_unit)"))
  co <- summary(fit)$coefficients["xz", ]
  tibble(level = "within_side", predictor = x, channel = ch, n = nrow(d),
         beta_std = co[["Estimate"]], se = co[["Std. Error"]], p = co[["Pr(>|t|)"]],
         singular = isSingular(fit))
}))

res <- bind_rows(bs, ws, wb) %>% group_by(level) %>% mutate(q = bh(p)) %>%
  ungroup() %>% mutate(sig_q05 = q < 0.05) %>% arrange(level, q)
write_outcome(res, A_DIR, "context_acoustic_screens.csv")
print(res %>% filter(p < 0.10) %>% select(level, predictor, channel, n, beta_std, p, q), n = 30)
cat("A08 done\n")
