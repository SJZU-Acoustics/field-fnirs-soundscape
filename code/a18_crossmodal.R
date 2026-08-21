# A18 — cross-modal moderation: does facing the natural side change how
# strongly the HEARD environment drives appraisal? A04 ran main-effect drivers
# per level; A02 found source saliences flat across sides. Never tested: the
# side x heard-fluctuation interaction — the audio-visual interaction proper.
# Family (declared): 8 tests = side x x_wb for 4 predictors (src_natural,
# src_traffic, src_human, laeq) x 2 outcomes (p_iso, e_iso), BH.
# Block-level LMMs with the predictor's _bs/_ws main effects controlled.

source(file.path("code", "helpers.R"))
A_DIR <- analysis_output_dir("a18_crossmodal_moderation")

blocks <- load_p27_blocks()
RE4 <- c("(1 | participant)", "(1 | site)", "(1 | sample_id)", "(1 | side_unit)")

run_ix <- function(d, x, y) {
  d <- decompose_bw(d, x)
  xbs <- paste0(x, "_bs"); xws <- paste0(x, "_ws"); xwb <- paste0(x, "_wb")
  f <- as.formula(paste0(y, " ~ side * ", xwb, " + ", xbs, " + ", xws, " + ",
                         paste(RE4, collapse = " + ")))
  fit <- suppressMessages(suppressWarnings(lmer(f, data = d, REML = TRUE)))
  co <- summary(fit)$coefficients
  term <- paste0("sidenatural:", xwb)
  wbterm <- xwb
  take <- function(nm) if (nm %in% rownames(co)) unlist(co[nm, c("Estimate", "Std. Error", "Pr(>|t|)")]) else c(NA, NA, NA)
  ix <- take(term); mn <- take(wbterm)
  tibble(predictor = x, outcome = y, n_rows = nrow(d), singular = isSingular(fit),
         wb_main_est = mn[1], wb_main_p = mn[3],
         inter_est = ix[1], inter_se = ix[2], inter_p = ix[3])
}

preds <- c("src_natural", "src_traffic", "src_human", "laeq")
res <- bind_rows(lapply(preds, function(x) {
  d <- if (x == "laeq") acoustic_rows(blocks) else blocks
  bind_rows(run_ix(d, x, "p_iso"), run_ix(d, x, "e_iso"))
})) %>% mutate(q = bh(inter_p))
write_outcome(res, A_DIR, "crossmodal_moderation_screen.csv")

print(res, n = 20)
cat("A18 done\n")
