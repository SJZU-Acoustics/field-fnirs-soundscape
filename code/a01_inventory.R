# A01 — inventory and conventions.
# What data each layer genuinely holds, the variance structure the models must
# respect, and the reliability of the block level within each layer.

source(file.path("code", "helpers.R"))
A_DIR <- analysis_output_dir("a01_inventory_and_conventions")

blocks  <- load_p27_blocks()
samples <- load_p27_samples()

# ---- 1. Layer availability accounting -------------------------------------
avail <- bind_rows(
  blocks %>% count(side, acoustic_recorded) %>% mutate(layer = "acoustic", flag = acoustic_recorded) %>% select(layer, side, flag, n),
  blocks %>% count(side, context_recorded)  %>% mutate(layer = "context",  flag = context_recorded)  %>% select(layer, side, flag, n),
  blocks %>% count(side, fixation_recorded) %>% mutate(layer = "fixation", flag = fixation_recorded) %>% select(layer, side, flag, n),
  blocks %>% count(side, fnirs_recorded)    %>% mutate(layer = "fnirs",    flag = fnirs_recorded)    %>% select(layer, side, flag, n)
)
write_outcome(avail, A_DIR, "layer_availability.csv")

# usable block-rows per layer after the layer rule
usable <- tibble(
  layer = c("perception", "acoustic", "context", "fixation", "fnirs_any"),
  usable_blocks = c(nrow(blocks), nrow(acoustic_rows(blocks)), nrow(context_rows(blocks)),
                    nrow(fixation_rows(blocks)), nrow(fnirs_rows(blocks)))
)
ch_usable <- map_dfr(CHANNELS$ch, function(ch) {
  tibble(layer = paste0("fnirs_", ch), usable_blocks = nrow(fnirs_channel_rows(blocks, ch)))
})
write_outcome(bind_rows(usable, ch_usable), A_DIR, "usable_rows.csv")

# fNIRS sample-side units with a valid channel
val <- load_p27_validity()
write_outcome(val %>% count(channel, valid), A_DIR, "channel_validity_counts.csv")

# ---- 2. Descriptives by side (block level, layer-filtered) -----------------
desc_var <- function(d, v, layer) {
  d %>% group_by(side) %>%
    summarise(n = sum(!is.na(.data[[v]])), mean = mean(.data[[v]], na.rm = TRUE),
              sd = sd(.data[[v]], na.rm = TRUE), min = min(.data[[v]], na.rm = TRUE),
              max = max(.data[[v]], na.rm = TRUE), .groups = "drop") %>%
    mutate(var = v, layer = layer, .before = 1)
}
desc <- bind_rows(
  map_dfr(c(PAQ_COLS, "p_iso", "e_iso", SRC_COLS), ~desc_var(blocks, .x, "perception")),
  map_dfr(ACOUSTIC_COLS, ~desc_var(acoustic_rows(blocks), .x, "acoustic")),
  map_dfr(CONTEXT_COLS,  ~desc_var(context_rows(blocks), .x, "context")),
  desc_var(fixation_rows(blocks), "fixation_s", "fixation"),
  map_dfr(HBO_COLS, function(v) {
    ch <- sub("hbo_", "", v)
    desc_var(fnirs_channel_rows(blocks, ch), v, "fnirs")
  })
)
write_outcome(desc, A_DIR, "descriptives_by_side.csv")

# ---- 3. Variance structure (crossed random intercepts, null models) --------
vc_for <- function(d, v) {
  fit <- suppressMessages(suppressWarnings(
    lmer(reformulate("1 + (1 | participant) + (1 | site) + (1 | sample_id) + (1 | side_unit)",
                     response = v), data = d %>% filter(!is.na(.data[[v]])))))
  vc <- as.data.frame(VarCorr(fit))
  tot <- sum(vc$vcov)
  vc %>% transmute(var = v, component = grp, vcov, share = vcov / tot) %>%
    mutate(singular = isSingular(fit))
}
vc_targets <- list(
  list(blocks, "p_iso"), list(blocks, "e_iso"),
  list(acoustic_rows(blocks), "laeq"),
  list(context_rows(blocks), "naturalness"), list(context_rows(blocks), "complexity"),
  list(fixation_rows(blocks), "fixation_s")
)
vc_out <- map_dfr(vc_targets, ~vc_for(.x[[1]], .x[[2]]))
vc_hbo <- map_dfr(CHANNELS$ch, function(ch) vc_for(fnirs_channel_rows(blocks, ch), paste0("hbo_", ch)))
write_outcome(bind_rows(vc_out, vc_hbo), A_DIR, "variance_structure.csv")

# ---- 4. Block-level reliability within sample-side -------------------------
# ICC(1) of single blocks within their sample-side unit, per layer: how much
# of the block-level variance is the unit's stable signal. Also Spearman-Brown
# reliability of the 3-block mean.
icc_blocks <- function(d, v) {
  dd <- d %>% filter(!is.na(.data[[v]])) %>% group_by(side_unit) %>%
    filter(n() >= 2) %>% ungroup()
  if (nrow(dd) < 10) return(tibble(var = v, n_units = 0, icc1 = NA, icc_mean3 = NA))
  fit <- suppressMessages(suppressWarnings(
    lmer(reformulate("1 + (1 | side_unit)", response = v), data = dd)))
  vc <- as.data.frame(VarCorr(fit))
  vu <- vc$vcov[vc$grp == "side_unit"]; vr <- vc$vcov[vc$grp == "Residual"]
  icc1 <- vu / (vu + vr)
  k <- mean(table(dd$side_unit))
  tibble(var = v, n_units = n_distinct(dd$side_unit), mean_k = k, icc1 = icc1,
         icc_meank = (k * icc1) / (1 + (k - 1) * icc1), singular = isSingular(fit))
}
rel <- bind_rows(
  map_dfr(c("p_iso", "e_iso", PAQ_COLS, SRC_COLS), ~icc_blocks(blocks, .x)),
  map_dfr(ACOUSTIC_COLS, ~icc_blocks(acoustic_rows(blocks), .x)),
  map_dfr(CONTEXT_COLS,  ~icc_blocks(context_rows(blocks), .x)),
  icc_blocks(fixation_rows(blocks), "fixation_s"),
  map_dfr(CHANNELS$ch, function(ch) icc_blocks(fnirs_channel_rows(blocks, ch), paste0("hbo_", ch)))
)
write_outcome(rel, A_DIR, "block_reliability.csv")

# ---- 5. Design covariates sanity ------------------------------------------
pos <- blocks %>% count(first_side, side, block, exposure_position)
write_outcome(pos, A_DIR, "exposure_position_map.csv")

samp_acct <- samples %>% group_by(sample_id) %>% summarise(
  n_sides = n(), fnirs_blocks = sum(n_blocks_fnirs), acoustic_blocks = sum(n_blocks_acoustic),
  context_blocks = sum(n_blocks_context), fixation_blocks = sum(n_blocks_fixation), .groups = "drop")
write_outcome(samp_acct, A_DIR, "sample_accounting.csv")

cat("A01 done\n")
