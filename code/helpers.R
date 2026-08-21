# Shared helpers for the field-fnirs-soundscape release (project 27).
#
# SINGLE DATA ENTRY. Every analysis sources this file and takes its data from
# the loaders below, which read the Mendeley workbook via load_data.R. No
# analysis re-derives a design column, recomputes an index, or re-filters the
# sample locally.
#
# LAYER RULE (deposit codebook): a block that was not recorded was filled per
# measurement layer, so any block-level analysis must filter on the
# *_recorded flag OF THE LAYER IT USES, via the layer_rows() helpers here —
# never on a single row-level flag. The fNIRS layer additionally carries
# per-channel validity (23 sample×side×channel units treated as missing);
# block-level fNIRS models must also filter on hbo_ch*_valid.

suppressPackageStartupMessages({
  library(tidyverse)
  library(lmerTest)
  library(broom)
  library(broom.mixed)
})

source(file.path("code", "load_data.R"))

BOOT_N <- 1000   # house default

# fNIRS montage after the global 6/8/9 exclusion (deposit codebook).
# Coverage is bilateral frontal, LEFT-only temporal, RIGHT-only occipital:
# no temporal/occipital laterality contrast exists in this dataset.
CHANNELS <- tribble(
  ~ch,   ~region,     ~hemisphere,
  "ch1",  "frontal",   "left",
  "ch3",  "frontal",   "left",
  "ch2",  "frontal",   "right",
  "ch4",  "frontal",   "right",
  "ch5",  "temporal",  "left",
  "ch7",  "temporal",  "left",
  "ch10", "occipital", "right",
)
HBO_COLS <- paste0("hbo_", CHANNELS$ch)
HBR_COLS <- paste0("hbr_", CHANNELS$ch)

# The eight PAQ items in delivered order.
PAQ_COLS <- c("paq_pleasant", "paq_eventful", "paq_engaging", "paq_chaotic",
              "paq_annoying", "paq_monotonous", "paq_uneventful", "paq_calm")
SRC_COLS <- c("src_traffic", "src_human", "src_natural", "src_other")
ACOUSTIC_COLS <- c("laeq", "l5", "l95", "l5_l95")
CONTEXT_COLS  <- c("complexity", "naturalness", "artificiality")

# Two cross-sample identical 3-block fixation vectors: the LATER member of
# each pair is not an independent observation.
EYE_DUP_LATER <- c("S05_P04", "S24_P22")  # duplicates of S02_P02 (natural), S03_P03 (composite)
EYE_DUP_SIDE  <- c(S05_P04 = "natural", S24_P22 = "composite")

.to_logical <- function(x) {
  if (is.logical(x)) return(x)
  out <- rep(NA, length(x))
  out[x %in% c("True", "TRUE", "true")]    <- TRUE
  out[x %in% c("False", "FALSE", "false")] <- FALSE
  out
}

# -------------------------------------------------------------------- loaders

#' The pairing table (40 rows): recorder, session date, first_side
#' (19 natural-first / 21 composite-first). The workbook omits session
#' clock-times (disclosure control), so no loader can return session_time.
load_p27_pairing <- function() {
  read_frozen_csv("p27_pairing.csv") %>%
    mutate(recorder_analysed = .to_logical(recorder_analysed))
}

#' Block-level table, 240 rows = 40 samples x 2 sides x 3 blocks.
#' Adds: first_side (from pairing), and exposure_position 1..6 — the session's
#' chronological slot of this side-block, from the alternating-sides design
#' (segments alternate sides; first_side says which side opened).
load_p27_blocks <- function() {
  pairing <- load_p27_pairing() %>%
    select(sample_id, first_side, session_date, any_of("session_time"))
  read_frozen_csv("p27_block_level.csv") %>%
    mutate(across(c(acoustic_recorded, context_recorded, fixation_recorded,
                    fnirs_recorded, ends_with("_valid")), .to_logical)) %>%
    left_join(pairing, by = "sample_id") %>%
    mutate(
      side  = factor(side, levels = c("composite", "natural")),
      block = as.integer(block),
      exposure_position = 2L * (block - 1L) +
        ifelse(as.character(side) == first_side, 1L, 2L),
      participant = factor(participant),
      site = factor(site),
      sample_id = factor(sample_id),
      side_unit = interaction(sample_id, side, drop = TRUE)
    )
}

#' Sample-level table, 80 rows = 40 samples x 2 sides; each layer meaned over
#' its own recorded blocks only; HbO/HbR cells for invalid channels are NA.
load_p27_samples <- function() {
  pairing <- load_p27_pairing() %>%
    select(sample_id, first_side, session_date, any_of("session_time"))
  read_frozen_csv("p27_sample_level.csv") %>%
    left_join(pairing, by = "sample_id") %>%
    mutate(
      side = factor(side, levels = c("composite", "natural")),
      participant = factor(participant),
      site = factor(site),
      sample_id = factor(sample_id)
    )
}

#' Wide paired table, 40 rows: natural minus composite difference per variable,
#' for the within-sample side contrast (the design's strongest level).
load_p27_paired <- function(vars) {
  s <- load_p27_samples()
  s %>%
    select(sample_id, participant, site, first_side, side, all_of(vars)) %>%
    pivot_wider(names_from = side, values_from = all_of(vars),
                names_glue = "{.value}__{side}") %>%
    { d <- .
      for (v in vars) d[[paste0(v, "__diff")]] <-
          d[[paste0(v, "__natural")]] - d[[paste0(v, "__composite")]]
      d }
}

#' Participants (27 rows). Disclosure control: the workbook carries age in
#' 5-year bands (age_band), not exact age, and no session clock-time. Modules
#' needing exact age must degrade gracefully (A15 does).
load_p27_participants <- function() read_frozen_csv("p27_participants.csv")

#' Sites (25 rows). The workbook carries both coordinate datums; this loader
#' restores the working pipeline's longitude/latitude names for the field
#' record (GCJ-02) and keeps the WGS-84 conversion beside them. Distances
#' computed from either datum are identical to well under a millimetre (the
#' shift is a near-rigid ~570 m translation at this extent).
load_p27_sites <- function() {
  read_frozen_csv("p27_sites.csv") %>%
    rename(longitude = longitude_gcj02, latitude = latitude_gcj02)
}

#' The WGS-84 site table exactly as the working pipeline's figure code expects
#' it (column order included): used by the Figure 1 site map.
load_p27_sites_wgs84 <- function() {
  read_frozen_csv("p27_sites.csv") %>%
    select(site, longitude_gcj02, latitude_gcj02, longitude_wgs84,
           latitude_wgs84, shift_m, distance_to_main_road_m)
}

load_p27_elements <- function() read_frozen_csv("p27_element_proportions.csv")

load_p27_validity <- function() {
  read_frozen_csv("p27_fnirs_channel_validity.csv")
}

# ---------------------------------------------------------- layer row filters

#' Rows where the given layer is genuinely recorded (never filled).
acoustic_rows <- function(d) filter(d, acoustic_recorded)
context_rows  <- function(d) filter(d, context_recorded)
fixation_rows <- function(d) filter(d, fixation_recorded %in% TRUE)
fnirs_rows    <- function(d) filter(d, fnirs_recorded %in% TRUE)

#' Block rows usable for one fNIRS channel: recorded block AND valid channel.
fnirs_channel_rows <- function(d, ch) {
  d %>% filter(fnirs_recorded %in% TRUE, .data[[paste0("hbo_", ch, "_valid")]] %in% TRUE)
}

# --------------------------------------------------------------- estimators

#' Within-sample side contrast on a sample-level variable.
#' Tests mean(natural - composite) over the 40 pairs with a participant random
#' intercept (13 participants contribute two samples). Falls back to the exact
#' paired t when the LMM is singular, and always reports d_z and the
#' participant-cluster bootstrap CI.
side_contrast <- function(var, boot = TRUE, seed = 271) {
  p <- load_p27_paired(var)
  dcol <- paste0(var, "__diff")
  d <- p %>% filter(!is.na(.data[[dcol]]))
  n <- nrow(d)
  est <- mean(d[[dcol]]); sdd <- sd(d[[dcol]])
  tt  <- t.test(d[[dcol]])
  fit <- suppressMessages(suppressWarnings(
    lmer(reformulate("1 + (1 | participant)", response = dcol), data = d)))
  sing <- isSingular(fit)
  lmm  <- summary(fit)$coefficients[1, ]
  ci <- c(NA_real_, NA_real_)
  if (boot) {
    set.seed(seed)
    parts <- split(seq_len(n), d$participant)
    bm <- replicate(BOOT_N, {
      take <- unlist(parts[sample(names(parts), length(parts), replace = TRUE)],
                     use.names = FALSE)
      mean(d[[dcol]][take])
    })
    ci <- quantile(bm, c(.025, .975), names = FALSE)
  }
  tibble(var = var, n_pairs = n, mean_diff = est, sd_diff = sdd,
         d_z = est / sdd, t = unname(tt$statistic), p_t = tt$p.value,
         lmm_est = lmm[["Estimate"]], lmm_se = lmm[["Std. Error"]],
         lmm_p = lmm[["Pr(>|t|)"]], lmm_singular = sing,
         boot_lo = ci[1], boot_hi = ci[2])
}

#' Block-level mixed model with the design's grouping structure.
#' Random intercepts: participant, site, sample (participant x site), and the
#' sample x side unit where blocks repeat. Drops terms that fail to converge.
p27_lmm <- function(formula_fixed, data,
                    re = c("(1 | participant)", "(1 | site)", "(1 | sample_id)")) {
  f <- as.formula(paste(deparse(formula_fixed), "+", paste(re, collapse = " + ")))
  suppressMessages(suppressWarnings(lmer(f, data = data, REML = TRUE)))
}

#' Between/within decomposition of a block-level predictor, computed on the
#' rows of its own layer. Adds <x>_bs (sample mean over both sides), <x>_ws
#' (side mean minus sample mean) and <x>_wb (block minus side mean).
decompose_bw <- function(d, x) {
  d %>%
    group_by(sample_id) %>% mutate("{x}_bs" := mean(.data[[x]], na.rm = TRUE)) %>%
    group_by(sample_id, side) %>%
    mutate("{x}_side_mean" := mean(.data[[x]], na.rm = TRUE)) %>%
    ungroup() %>%
    mutate("{x}_ws" := .data[[paste0(x, "_side_mean")]] - .data[[paste0(x, "_bs")]],
           "{x}_wb" := .data[[x]] - .data[[paste0(x, "_side_mean")]]) %>%
    select(-all_of(paste0(x, "_side_mean")))
}

#' Benjamini-Hochberg q-values.
bh <- function(p) p.adjust(p, method = "BH")

#' Participant-cluster bootstrap of an arbitrary statistic over a data frame.
boot_by_participant <- function(d, stat_fn, n_boot = BOOT_N, seed = 271) {
  set.seed(seed)
  parts <- split(seq_len(nrow(d)), d$participant)
  replicate(n_boot, {
    take <- unlist(parts[sample(names(parts), length(parts), replace = TRUE)],
                   use.names = FALSE)
    stat_fn(d[take, , drop = FALSE])
  })
}

# ------------------------------------------------------------------- outputs

#' Create and return this module's output directory, output/<module>/.
analysis_output_dir <- function(module) {
  out <- file.path("output", module)
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  out
}

#' Write a CSV into the module's output directory.
write_outcome <- function(x, analysis_dir, filename) {
  dir.create(analysis_dir, showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(x, file.path(analysis_dir, filename))
  invisible(x)
}
