# =============================================================================
# Project 27 — generate all LaTeX table fragments from the exploration outcomes.
# Main-text fragments (tabular only, captions live in main.tex):
#   tables/table1_side_contrasts.tex, tables/table2_driver_screens.tex
# SI tables (complete environments incl. captions, \input in order):
#   supplementary/tabS_*.tex
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(readr)})

# Release paths: analysis modules write to output/aNN_<slug>/ (the working
# pipeline's exploration/analysis_NN_<slug>/outcome/), tables to output/.
TDIR  <- file.path("output", "tables")
SDIR  <- file.path("output", "supplementary")
dir.create(TDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(SDIR, showWarnings = FALSE, recursive = TRUE)
OUTC  <- function(a, f) read_csv(file.path("output", sub("^analysis_", "a", a), f),
                                 show_col_types = FALSE)

lab <- c(paq_pleasant = "pleasant", paq_eventful = "eventful", paq_engaging = "vibrant",
         paq_chaotic = "chaotic", paq_annoying = "annoying", paq_monotonous = "monotonous",
         paq_uneventful = "uneventful", paq_calm = "calm",
         p_iso = "ISO pleasantness", e_iso = "ISO eventfulness",
         src_traffic = "traffic sounds", src_human = "human sounds",
         src_natural = "natural sounds", src_other = "other sounds",
         laeq = "$L_{\\mathrm{Aeq}}$ (dB(A))", l5 = "$L_{5}$ (dB(A))",
         l95 = "$L_{95}$ (dB(A))", l5_l95 = "$L_{5}-L_{95}$ (dB(A))",
         complexity = "visual complexity", naturalness = "viewed naturalness",
         artificiality = "viewed artificiality", fixation_s = "total fixation duration (s)",
         hbo_ch1 = "HbO channel 1", hbo_ch2 = "HbO channel 2", hbo_ch3 = "HbO channel 3",
         hbo_ch4 = "HbO channel 4", hbo_ch5 = "HbO channel 5", hbo_ch7 = "HbO channel 7",
         hbo_ch10 = "HbO channel 10",
         hbr_ch1 = "HbR channel 1", hbr_ch2 = "HbR channel 2", hbr_ch3 = "HbR channel 3",
         hbr_ch4 = "HbR channel 4", hbr_ch5 = "HbR channel 5", hbr_ch7 = "HbR channel 7",
         hbr_ch10 = "HbR channel 10",
         hbt_ch1 = "HbT channel 1", hbt_ch3 = "HbT channel 3", hbt_ch4 = "HbT channel 4",
         hbo_frontal = "frontal HbO", d_piso = "$\\Delta$ISO pleasantness",
         d_eiso = "$\\Delta$ISO eventfulness", d_mono = "$\\Delta$monotonous",
         d_hbof = "$\\Delta$frontal HbO", d_nat = "$\\Delta$naturalness",
         mono = "monotonous", monotony = "monotonous", P_ISO = "ISO pleasantness")
L <- function(x) {
  x <- as.character(x)
  ifelse(x %in% names(lab), unname(lab[x]), gsub("_", "\\\\_", x))
}

fnum <- function(x, d = 3) ifelse(is.na(x), "--", formatC(x, format = "f", digits = d))
fp <- function(p) ifelse(is.na(p), "--",
                  ifelse(p < 0.001, "$<$0.001", formatC(p, format = "f", digits = 3)))
fci <- function(lo, hi, d = 2) sprintf("[%s, %s]", fnum(lo, d), fnum(hi, d))

rows_to_tabular <- function(mat, align, header) {
  body <- apply(mat, 1, paste, collapse = " & ")
  paste0("\\begin{tabular}{", align, "}\n\\toprule\n", header, " \\\\\n\\midrule\n",
         paste0(body, " \\\\", collapse = "\n"), "\n\\bottomrule\n\\end{tabular}\n")
}

write_frag <- function(txt, path) {cat(txt, file = path); message("wrote ", path)}

si_table <- function(mat, align, header, caption, label, path, long = FALSE,
                     size = "\\footnotesize") {
  if (!long) {
    txt <- paste0("\\begin{table}[htbp]\n\\centering\n", size, "\n\\caption{", caption,
                  "}\n\\label{", label, "}\n",
                  rows_to_tabular(mat, align, header), "\\end{table}\n")
  } else {
    body <- apply(mat, 1, paste, collapse = " & ")
    txt <- paste0("{", size, "\n\\begin{longtable}{", align, "}\n\\caption{", caption,
                  "}\\label{", label, "}\\\\\n\\toprule\n", header,
                  " \\\\\n\\midrule\n\\endfirsthead\n\\toprule\n", header,
                  " \\\\\n\\midrule\n\\endhead\n\\bottomrule\n\\endfoot\n",
                  paste0(body, " \\\\", collapse = "\n"), "\n\\end{longtable}\n}\n")
  }
  write_frag(txt, path)
}

## ============================ MAIN TABLE 1 ===================================
a02 <- OUTC("analysis_02_side_contrast", "side_contrast_families.csv")
a06 <- OUTC("analysis_06_fnirs_side_contrast", "fnirs_side_contrast.csv")

t1_block <- function(df, dq = 2) {
  cbind(L(df$var), df$n_pairs, fnum(df$mean_diff, dq), fci(df$boot_lo, df$boot_hi, dq),
        fnum(df$d_z, 2), fp(df$p_t), fp(df$q))
}
ord_perc <- c("paq_pleasant","paq_eventful","paq_engaging","paq_chaotic","paq_annoying",
              "paq_monotonous","paq_uneventful","paq_calm","p_iso","e_iso",
              "src_traffic","src_human","src_natural","src_other")
perc <- a02 %>% filter(family == "F1_perception") %>%
  mutate(var = factor(var, levels = ord_perc)) %>% arrange(var)
phys <- a02 %>% filter(family == "F2_physical_context") %>%
  mutate(var = factor(var, levels = c("laeq","l5","l95","l5_l95","complexity",
                                      "naturalness","artificiality"))) %>% arrange(var)
fixr <- a02 %>% filter(family == "F3_fixation")
hbo  <- a06 %>% filter(chromophore == "HbO") %>%
  mutate(var = factor(var, levels = paste0("hbo_ch", c(1,3,2,4,5,7,10)))) %>% arrange(var)

hdr1 <- "Measure & $n$ & Difference & 95\\% CI & $d_z$ & $p$ & $q$"
gh <- function(t) sprintf("\\multicolumn{7}{l}{\\textbf{%s}} \\\\", t)
b1 <- apply(t1_block(perc), 1, paste, collapse = " & ")
b2 <- apply(t1_block(phys), 1, paste, collapse = " & ")
b3 <- apply(t1_block(fixr, 2), 1, paste, collapse = " & ")
b4 <- apply(t1_block(hbo), 1, paste, collapse = " & ")
t1 <- paste0("\\begin{tabular}{lcccccc}\n\\toprule\n", hdr1, " \\\\\n\\midrule\n",
  gh("Perception family (14 tests, BH within family)"), "\n",
  paste0(b1, " \\\\", collapse = "\n"), "\n\\midrule\n",
  gh("Physical and scene family (7 tests, BH within family)"), "\n",
  paste0(b2, " \\\\", collapse = "\n"), "\n\\midrule\n",
  gh("Visual attention (single test, native $p$)"), "\n",
  paste0(b3, " \\\\", collapse = "\n"), "\n\\midrule\n",
  gh("fNIRS family (7 HbO tests, BH within family)"), "\n",
  paste0(b4, " \\\\", collapse = "\n"), "\n\\bottomrule\n\\end{tabular}\n")
write_frag(t1, file.path(TDIR, "table1_side_contrasts.tex"))

## ============================ MAIN TABLE 2 ===================================
a04 <- OUTC("analysis_04_perception_drivers", "driver_screens.csv")
bs <- a04 %>% filter(level == "between_sample") %>%
  arrange(desc(outcome), q) %>%
  transmute(o = L(outcome), p_ = L(predictor), b = fnum(beta_std, 2), s = fnum(se, 2),
            pp = fp(p), qq = fp(q))
t2 <- rows_to_tabular(as.matrix(bs), "llcccc",
  "Outcome & Predictor & $\\beta$ & SE & $p$ & $q$")
write_frag(t2, file.path(TDIR, "table2_driver_screens.tex"))

## ============================ SI TABLES ======================================
# S1 reliability
# ICC values from A01; unit counts and mean-k corrected by A25 (post-review):
# A01's n_units are side-units with >= 2 recorded blocks (ICC-estimable), and
# its mean_k divided those units' blocks by all 80 units. The corrected table
# shows both unit counts and the mean k within the ICC-estimable units.
a01 <- OUTC("analysis_01_inventory_and_conventions", "block_reliability.csv")
a25r <- OUTC("analysis_25_post_review_sensitivities", "reliability_units_corrected.csv")
a01 <- left_join(a01, a25r, by = c(var = "measure"))
stopifnot(identical(a01$n_units, a01$units_ge2))
si_table(cbind(L(a01$var), a01$units_ge1, a01$units_ge2,
               fnum(a01$mean_k_within_ge2, 2), fnum(a01$icc1, 2),
               fnum(a01$icc_meank, 2)), "lccccc",
  "Measure & $n_{\\geq 1}$ & $n_{\\geq 2}$ & Mean blocks & ICC(1) & ICC(mean $k$)",
  paste0("Between-block reliability of every measurement layer at the sample--side level. ",
         "$n_{\\geq 1}$ counts sample $\\times$ side cells with at least one recorded ",
         "block in that layer (the cells contributing to the paired contrasts in ",
         "Table~1); $n_{\\geq 2}$ counts cells with at least two recorded blocks, ",
         "over which the intraclass correlations are estimated, with the mean number of ",
         "recorded blocks within those cells. ICC(1) is the single-block intraclass ",
         "correlation from a one-way random-effects decomposition over side-units; ",
         "ICC(mean $k$) is the reliability of the mean over each unit's recorded blocks."),
  "tab:s-reliability", file.path(SDIR, "tabS01_reliability.tex"))

# S2 activation vs baseline
a05a <- OUTC("analysis_05_fnirs_signal_anatomy", "activation_vs_baseline.csv")
si_table(cbind(a05a$channel, a05a$chromophore, a05a$n_units, fnum(a05a$est, 3),
               fnum(a05a$se, 3), fp(a05a$p), ifelse(is.na(a05a$q), "--", fp(a05a$q))),
  "llccccc", "Channel & Chromophore & $n$ & Estimate & SE & $p$ & $q$",
  paste0("Net exposure-minus-baseline shift per channel (linear mixed model with participant ",
         "random intercept; $\\mu$mol/L). The HbO family (7 tests) is Benjamini--Hochberg ",
         "corrected; HbR is convergent evidence and uncorrected. No channel is significant at ",
         "$q<0.05$: exposure does not shift oxygenation from its masked baseline on average."),
  "tab:s-activation", file.path(SDIR, "tabS02_activation.tex"))

# S3 HbO-HbR coupling
a05c <- OUTC("analysis_05_fnirs_signal_anatomy", "hbo_hbr_coupling.csv")
si_table(cbind(a05c$channel, a05c$n, fnum(a05c$r_hbo_hbr, 2)), "lcc",
  "Channel & $n$ & $r$(HbO, HbR)",
  paste0("HbO--HbR coupling per channel across sample--side units. Canonical neurovascular ",
         "coupling implies a negative correlation; channel 5 is the only channel with a ",
         "clearly canonical value, and the strong positive coupling of channel 10 marks it ",
         "as systemic-suspect."),
  "tab:s-coupling-anatomy", file.path(SDIR, "tabS03_coupling_anatomy.tex"))

# S4 driver screens, all levels
a04f <- a04 %>%
  transmute(lv = gsub("_", " ", level), o = L(outcome), p_ = L(predictor),
            b = fnum(beta_std, 2), s = fnum(se, 2), pp = fp(p), qq = fp(q))
si_table(as.matrix(a04f), "lllcccc",
  "Level & Outcome & Predictor & $\\beta$ & SE & $p$ & $q$",
  paste0("The full three-level perception driver screen (16 tests per level: 8 predictors ",
         "$\\times$ 2 outcomes; BH within level). Predictors enter decomposed: between-sample ",
         "(cross-sectional), between-side (within-sample side difference) and within-side ",
         "(between-block fluctuation). Standardised slopes from mixed models with the design ",
         "random effects."),
  "tab:s-drivers", file.path(SDIR, "tabS04_driver_screens_full.tex"), long = TRUE, size = "\\scriptsize")

# S5 brain-perception coupling screens
a07 <- OUTC("analysis_07_brain_perception_coupling", "coupling_screens.csv")
a07f <- a07 %>% transmute(lv = gsub("_", " ", level), ax = L(axis), ch = channel,
                          n, b = fnum(beta_std, 2), s = fnum(se, 2), pp = fp(p), qq = fp(q))
si_table(as.matrix(a07f), "lllccccc",
  "Level & Axis & Channel & $n$ & $\\beta$ & SE & $p$ & $q$",
  paste0("Brain--perception coupling screens (per level, 14 tests: 2 ISO axes $\\times$ 7 HbO ",
         "channels; BH within level). No coupling is significant at any level: appraisal and ",
         "oxygenation shift in parallel without measurable per-pair covariation."),
  "tab:s-brain-perception", file.path(SDIR, "tabS05_brain_perception.tex"), long = TRUE)

# S6 brain-context/acoustic screens
a08 <- OUTC("analysis_08_brain_context_acoustic_coupling", "context_acoustic_screens.csv")
a08f <- a08 %>% transmute(lv = gsub("_", " ", level), pr = L(predictor), ch = channel,
                          n, b = fnum(beta_std, 2), s = fnum(se, 2), pp = fp(p), qq = fp(q))
si_table(as.matrix(a08f), "lllccccc",
  "Level & Predictor & Channel & $n$ & $\\beta$ & SE & $p$ & $q$",
  paste0("Brain--context and brain--acoustic coupling screens (per level, 21 tests: 3 ",
         "predictors $\\times$ 7 HbO channels; BH within level). Oxygenation is not ",
         "continuously related to measured naturalness, visual complexity or ",
         "$L_{\\mathrm{Aeq}}$ at any level."),
  "tab:s-brain-context", file.path(SDIR, "tabS06_brain_context.tex"), long = TRUE)

# S7 fixation strand
a10b <- OUTC("analysis_10_fixation_strand", "between_side_coupling.csv") %>% filter(set == "all")
a10w <- OUTC("analysis_10_fixation_strand", "within_side_coupling.csv")
m7 <- rbind(
  cbind("Between-side", L(a10b$target), a10b$n, fnum(a10b$r, 2), "--", fp(a10b$p), fp(a10b$q)),
  cbind("Within-side", L(a10w$target), a10w$n, fnum(a10w$beta_std, 2), fnum(a10w$se, 2),
        fp(a10w$p), fp(a10w$q)))
si_table(m7, "llccccc",
  "Level & Target & $n$ & Estimate & SE & $p$ & $q$",
  paste0("Fixation coupling screens (two BH families of four tests). Between-side rows are ",
         "correlations of $\\Delta$fixation with the paired shifts; within-side rows are ",
         "standardised block-level slopes. The estimate column holds $r$ (between-side) or ",
         "$\\beta$ (within-side). All null; conclusions unchanged when the later member of ",
         "each duplicated fixation vector is excluded."),
  "tab:s-fixation", file.path(SDIR, "tabS07_fixation.tex"))

# S8 first-side moderation
a11f <- OUTC("analysis_11_order_and_time", "first_side_moderation.csv")
si_table(cbind(L(a11f$var), a11f$n, fnum(a11f$side_main, 2), fnum(a11f$inter_est, 2),
               fnum(a11f$inter_se, 2), fp(a11f$inter_p), fp(a11f$q)), "lcccccc",
  "Outcome & $n$ rows & Side main effect & Side $\\times$ order & SE & $p$ & $q$",
  paste0("Exposure-order moderation of the four main side contrasts (side $\\times$ ",
         "first-side interaction, block-level mixed models; BH over 4 tests). No main contrast ",
         "depends on which side a session opened with."),
  "tab:s-first-side", file.path(SDIR, "tabS08_first_side.tex"))

# S9 exposure-position trends
a11p <- OUTC("analysis_11_order_and_time", "exposure_position_trends.csv")
si_table(cbind(L(a11p$var), a11p$n, fnum(a11p$pos_beta, 3), fnum(a11p$se, 3), fp(a11p$p),
               fp(a11p$q)), "lccccc",
  "Outcome & $n$ rows & Slope per position & SE & $p$ & $q$",
  paste0("Time-on-task trends across a session's six exposures (standardised slope on ",
         "exposure position 1--6, design random effects; BH over 9 tests). ISO pleasantness ",
         "declines across the session; the decline is orthogonal to the side contrast by the ",
         "alternating design."),
  "tab:s-position", file.path(SDIR, "tabS09_position_trends.tex"))

# S10 time course (side x block)
a16 <- OUTC("analysis_16_effect_time_course", "side_x_block_interactions.csv")
si_table(cbind(L(a16$outcome), a16$n_rows, fnum(a16$side_est, 3), fp(a16$side_p),
               fnum(a16$inter_est, 3), fnum(a16$inter_se, 3), fp(a16$inter_p),
               fp(a16$q_inter)), "lccccccc",
  "Outcome & $n$ & Side effect & Side $p$ & Side $\\times$ blk & SE & $p$ & $q$",
  paste0("Within-side time course of the side contrast (side $\\times$ block interaction, ",
         "block centred; BH over 4 tests). The side effect is present from the first 60-s ",
         "block and does not grow or fade across a side's three repetitions."),
  "tab:s-time-course", file.path(SDIR, "tabS10_time_course.tex"))

# S11 per-block profile (descriptive)
a16p <- OUTC("analysis_16_effect_time_course", "per_block_side_difference_profile.csv")
si_table(cbind(a16p$block, fnum(a16p$d_piso, 3), fnum(a16p$d_eiso, 3), fnum(a16p$d_mono, 3),
               fnum(a16p$d_hbof, 3)), "lcccc",
  paste0("Block & $\\Delta$ISO pleasantness & $\\Delta$ISO eventfulness & ",
         "$\\Delta$monotonous & $\\Delta$frontal HbO"),
  "Per-block natural $-$ composite differences (descriptive; recorded blocks only).",
  "tab:s-block-profile", file.path(SDIR, "tabS11_block_profile.tex"))

# S12 lateralisation
a17 <- OUTC("analysis_17_frontal_lateralisation", "frontal_lateralisation_tests.csv")
si_table(cbind(a17$test, a17$n_pairs, fnum(a17$mean, 2), fci(a17$boot_lo, a17$boot_hi, 2),
               fnum(a17$d_z, 2), fp(a17$p_t), fp(a17$q)), "lcccccc",
  "Contrast & $n$ & Mean & 95\\% CI & $d_z$ & $p$ & $q$",
  paste0("Frontal lateralisation tests (BH over 3 tests): side contrasts on the left- and ",
         "right-frontal HbO means and their difference-of-differences, the montage's only ",
         "possible hemisphere comparison. The drop is statistically symmetric."),
  "tab:s-lateralisation", file.path(SDIR, "tabS12_lateralisation.tex"))

# S13 cross-modal moderation
a18 <- OUTC("analysis_18_crossmodal_moderation", "crossmodal_moderation_screen.csv")
si_table(cbind(L(a18$predictor), L(a18$outcome), a18$n_rows, fnum(a18$wb_main_est, 3),
               fp(a18$wb_main_p), fnum(a18$inter_est, 3), fnum(a18$inter_se, 3),
               fp(a18$inter_p), fp(a18$q)), "llccccccc",
  paste0("Heard predictor & Outcome & $n$ & W-side slope & $p$ & Side $\\times$ ",
         "slope & SE & $p$ & $q$"),
  paste0("Audio--visual interaction screen (8 tests: 4 heard-environment predictors ",
         "$\\times$ 2 outcomes; BH). The heard$\\rightarrow$appraisal slopes are identical ",
         "on the two sides: visual context moves the appraisal intercept and leaves hearing ",
         "slopes untouched."),
  "tab:s-crossmodal", file.path(SDIR, "tabS13_crossmodal.tex"), size = "\\scriptsize")

# S14 element-class screen
a14 <- OUTC("analysis_14_element_class_screen", "element_class_screen.csv")
si_table(cbind(L(a14$outcome), gsub("_", " ", a14$class_group), a14$n, fnum(a14$r, 2),
               fp(a14$p), fp(a14$q)), "llcccc",
  "Outcome shift & Class group & $n$ & $r$ & $p$ & $q$",
  paste0("Element-class screen (16 tests: 8 aggregated semantic-segmentation class groups ",
         "$\\times$ 2 outcome shifts, between-side level; BH as one family). No class ",
         "group's between-side change predicts either outcome shift."),
  "tab:s-element", file.path(SDIR, "tabS14_element_screen.tex"))

# S15 scene-content moderators (A19 + A22)
a19 <- OUTC("analysis_19_scene_content_moderators", "content_moderator_screen.csv")
a22 <- OUTC("analysis_22_vegetation_structure", "vegetation_structure_tests.csv")
m15 <- rbind(
  cbind("Content screen (6 tests)",
        paste(gsub("_nat", " (natural side)", a19$moderator), "$\\times$", L(a19$outcome)),
        a19$n_pairs, fnum(a19$r, 2), fp(a19$p), fp(a19$q)),
  cbind("Vegetation structure (4 tests)", gsub("_nat", " (natural side)", a22$test),
        a22$n, fnum(a22$r, 2), fp(a22$p), fp(a22$q)))
m15[, 2] <- gsub("x dP_ISO", "$\\\\times$ $\\\\Delta$ISO pleasantness", m15[, 2])
m15[, 2] <- gsub("x dMonotony", "$\\\\times$ $\\\\Delta$monotonous", m15[, 2])
m15[, 2] <- gsub("d_piso", "$\\\\Delta$ISO pleasantness", m15[, 2])
m15[, 2] <- gsub("d_hbof", "$\\\\Delta$frontal HbO", m15[, 2])
m15[, 2] <- gsub("\\(SITE level\\)", "(site level)", m15[, 2])
si_table(m15, "p{2.9cm}p{4.6cm}cccc", "Family & Test & $n$ & $r$ & $p$ & $q$",
  paste0("Faced-scene content moderators of the appraisal gain. The content screen tests ",
         "the natural-side water, grass and sky proportions against both outcome shifts ",
         "(6 tests, BH); the vegetation-structure family separates ground-plane lawn from ",
         "trees and other vegetation and re-tests the lawn moderator at the site level ",
         "(4 tests, BH)."),
  "tab:s-scene-content", file.path(SDIR, "tabS15_scene_content.tex"))

# S16 person level
a15r <- OUTC("analysis_15_person_level", "per_person_effect_reliability.csv")
a15t <- OUTC("analysis_15_person_level", "trait_moderation.csv")
# Withheld input: exact age is banded in the deposit, so the two age-moderation
# rows cannot be recomputed. The fragment is not regenerated; the printed
# article version (reliability rows + sex and age moderation) stands.
if (!any(a15t$trait == "age")) {
  message("tabS16_person.tex not regenerated: age-moderation rows need exact ",
          "age, which the deposit bands (disclosure control).")
} else {
m16 <- rbind(
  cbind("Split-sample reliability", L(a15r$measure), a15r$n_participants, fnum(a15r$r, 2),
        "--", fp(a15r$p), "--"),
  cbind("Trait moderation", paste(L(a15t$outcome), "$\\times$", a15t$trait), a15t$n,
        fnum(a15t$est, 2), fnum(a15t$se, 2), fp(a15t$p), fp(a15t$q)))
si_table(m16, "p{2.3cm}p{3.5cm}ccccc", "Family & Test & $n$ & Estimate & SE & $p$ & $q$",
  paste0("Person-level reliability and moderation. Reliability rows: correlation between the two per-sample effects ",
         "of the 13 participants measured at two sites (the precondition for any responder ",
         "analysis; estimate is $r$). Moderation rows: sex and age moderation of the two main ",
         "shifts (4 tests, BH)."),
  "tab:s-person", file.path(SDIR, "tabS16_person.tex"))
}

# S17 site-level structure and dose
a20 <- OUTC("analysis_20_site_structure", "site_structure_tests.csv")
a12d <- OUTC("analysis_12_site_generality", "site_dose_response.csv")
a12h <- OUTC("analysis_12_site_generality", "site_heterogeneity_permutation.csv")
m17 <- rbind(
  cbind("Spatial and cross-outcome (BH over 3 tests)", a20$test, a20$n_sites, fnum(a20$r_obs, 2),
        fp(a20$p_perm), fp(a20$q)),
  cbind("Site dose-response (BH over 4 tests)",
        paste(L(a12d$outcome), "vs", ifelse(a12d$predictor == "d_nat",
              "$\\Delta$naturalness", "distance to main road")),
        a12d$n_sites, fnum(a12d$r, 2), fp(a12d$p), fp(a12d$q)),
  cbind("Site heterogeneity (permutation)", L(a12h$outcome), a12h$n_sites,
        fnum(a12h$stat_obs, 2), fp(a12h$p_perm), "--"))
m17[, 2] <- gsub("dP_ISO", "$\\\\Delta$ISO pleasantness", m17[, 2])
m17[, 2] <- gsub("dFrontal-HbO", "$\\\\Delta$frontal HbO", m17[, 2])
m17[, 2] <- gsub(" x ", " $\\\\times$ ", m17[, 2])
m17[, 2] <- gsub("spatial:", "spatial autocorrelation,", m17[, 2])
si_table(m17, "p{2.7cm}p{3.8cm}cccc", "Family & Test & Sites & Statistic & $p$ & $q$",
  paste0("Site-level structure. Spatial rows: Mantel-style correlation of pairwise effect ",
         "differences with geographic distance (1,000 site-label permutations). Dose rows: ",
         "site-mean effects against the naturalness gap and road distance. Heterogeneity ",
         "rows: random side-slope variance against a within-sample permutation null (the ",
         "statistic is the observed slope SD for ISO pleasantness and the likelihood-ratio ",
         "statistic for frontal HbO)."),
  "tab:s-site", file.path(SDIR, "tabS17_site_structure.tex"))

# S18 sensitivity battery
a13w <- OUTC("analysis_13_sensitivity_robustness", "wilcoxon_headlines.csv")
a13t <- OUTC("analysis_13_sensitivity_robustness", "hbt_convergence.csv")
a13s <- OUTC("analysis_13_sensitivity_robustness", "strict_acoustic_pairs.csv")
a13l <- OUTC("analysis_13_sensitivity_robustness", "loo_influence.csv")
a13p <- OUTC("analysis_13_sensitivity_robustness", "position_trend_sensitivity.csv")
m18 <- rbind(
  cbind("Wilcoxon signed-rank", L(a13w$var), a13w$n, fnum(a13w$median_diff, 2), "--",
        fp(a13w$p_wilcox)),
  cbind("HbT convergence", L(a13t$var), a13t$n, fnum(a13t$mean_diff, 2), fnum(a13t$d_z, 2),
        fp(a13t$p)),
  cbind("Strict acoustic pairs", L(a13s$var), a13s$n_pairs, fnum(a13s$mean_diff, 2),
        fnum(a13s$d_z, 2), fp(a13s$p)),
  cbind("Leave-one-participant-out", L(a13l$var), "--",
        sprintf("%s [%s, %s]", fnum(a13l$est_full, 2), fnum(a13l$est_min, 2),
                fnum(a13l$est_max, 2)), "--", fp(a13l$p_max)),
  cbind("Position trend without exposure 1", "ISO pleasantness", 200,
        fnum(a13p$beta, 3), "--", fp(a13p$p)))
a25t <- OUTC("analysis_25_post_review_sensitivities", "tost_acoustic_equivalence.csv")
a25s <- OUTC("analysis_25_post_review_sensitivities", "site_re_sensitivity.csv")
a25c <- OUTC("analysis_25_post_review_sensitivities", "matched_channel_frontal.csv")
a25c <- a25c[a25c$measure %in% c("frontal mean, matched complete-4 channels",
                                 "left minus right (diff-of-diffs), complete-4"), ]
m18 <- rbind(m18,
  cbind("TOST equivalence, $\\pm$1.5 dB(A)", L(a25t$var), a25t$n,
        sprintf("%s [%s, %s]", fnum(a25t$mean_diff, 2), fnum(a25t$ci90_lo, 2),
                fnum(a25t$ci90_hi, 2)), "--", fp(a25t$p_tost)),
  cbind("Site random intercept added", L(a25s$var), a25s$n_pairs,
        fnum(a25s$est_plus_site_re, 2), "--",
        paste0(fp(a25s$p_plus_site_re), ifelse(a25s$singular, "$^{\\dagger}$", ""))),
  cbind("Matched complete-4 frontal channels",
        c("frontal HbO", "left $-$ right frontal HbO"), a25c$n,
        fnum(a25c$mean_diff, 2), fnum(a25c$d_z, 2), fp(a25c$p)))
si_table(m18, "p{2.6cm}p{2.3cm}cccc", "Check & Measure & $n$ & Estimate & $d_z$ & $p$",
  paste0("Sensitivity analyses for every main estimate. Wilcoxon rows report median ",
         "paired differences; HbT rows the total-haemoglobin convergence on the frontal ",
         "channels; strict rows re-estimate on pairs whose acoustic layer is recorded on ",
         "both sides; leave-one-out rows give the full estimate with its range over ",
         "single-participant exclusions and the largest resulting $p$; the position-trend ",
         "row re-fits the session decline without the first exposure. The last three row ",
         "groups were added after an internal review of the complete draft (post-review ",
         "sensitivities, uncorrected): two one-sided tests of equivalence within ",
         "$\\pm$1.5 dB(A) (90\\% CI shown; the margin is half the $\\approx$3 dB(A) level ",
         "change conventionally taken as just noticeable for broadband environmental ",
         "sound); the primary contrasts re-tested with a site random intercept added to ",
         "the participant random intercept ($^{\\dagger}$site variance estimated at zero, ",
         "singular fit); and the frontal regional contrasts recomputed on pairs with all ",
         "four frontal channels valid on both sides."),
  "tab:s-sensitivity", file.path(SDIR, "tabS18_sensitivity.tex"))

# S19 decoding readouts
a09 <- OUTC("analysis_09_multivariate_pattern", "multivariate_readouts.csv")
a13d <- OUTC("analysis_13_sensitivity_robustness", "decoding_6ch.csv")
a23 <- OUTC("analysis_23_perception_decoding", "perception_decoding_readouts.csv")
m19 <- rbind(
  cbind("Side from 7-channel $\\Delta$HbO", a09$n_pairs[1], fnum(a09$statistic[1], 3),
        fnum(a09$null_q95[1], 3), fp(a09$p_perm[1])),
  cbind("$\\Delta$ISO pleasantness from $\\Delta$HbO ($r$)", a09$n_pairs[2],
        fnum(a09$statistic[2], 3), fnum(a09$null_q95[2], 3), fp(a09$p_perm[2])),
  cbind("Side from 6-channel $\\Delta$HbO (channel 2 dropped)", a13d$n_pairs,
        fnum(a13d$accuracy, 3), fnum(a13d$null_q95, 3), fp(a13d$p_perm)),
  cbind("Side from 8-item $\\Delta$PAQ profile", a23$n_pairs[1], fnum(a23$accuracy[1], 3),
        fnum(a23$null_q95[1], 3), fp(a23$p_perm[1])),
  cbind("Side from $\\Delta$PAQ $+$ $\\Delta$HbO", a23$n_pairs[2], fnum(a23$accuracy[2], 3),
        fnum(a23$null_q95[2], 3), fp(a23$p_perm[2])))
a25d <- OUTC("analysis_25_post_review_sensitivities", "matched_decoding.csv")
m19 <- rbind(m19,
  cbind("Side from $\\Delta$PAQ, matched 30 pairs", 30, fnum(a25d$statistic[2], 3),
        fnum(a25d$null_q95[2], 3), fp(a25d$p[2])),
  cbind("$\\Delta$HbO $-$ $\\Delta$PAQ accuracy, matched 30 pairs", 30,
        fnum(a25d$statistic[3], 3), "--", fp(a25d$p[3])))
si_table(m19, "p{5.1cm}cccc",
  "Analysis & $n$ pairs & Statistic & Null 95th pct. & Perm. $p$",
  paste0("Multivariate decoding analyses (leave-one-participant-out Fisher discriminant with ",
         "participant-level folds; 1,000 permutations each, sign-flipping the paired ",
         "differences so the permutation preserves the pairing). Accuracy for ",
         "the classification analyses; out-of-sample correlation for the regression analysis. ",
         "The last two rows were added post-review: the $\\Delta$PAQ classifier re-run on ",
         "the same 30 complete-channel pairs as the $\\Delta$HbO classifier, and the ",
         "paired accuracy difference on those pairs (McNemar exact test on discordant ",
         "held-out predictions, 11 vs 3)."),
  "tab:s-decoding", file.path(SDIR, "tabS19_decoding.tex"))

# S20 first-block between-group analysis
a21 <- OUTC("analysis_21_naive_first_exposure", "first_exposure_between_groups.csv")
a25f <- OUTC("analysis_25_post_review_sensitivities", "first_session_subset.csv")
m20 <- rbind(
  cbind(L(a21$outcome), a21$n, fnum(a21$mean_natural_first, 2),
        fnum(a21$mean_composite_first, 2), fnum(a21$lmm_est, 2),
        fnum(a21$lmm_se, 2), fp(a21$lmm_p), fp(a21$q)),
  cbind(paste0(L(a25f$outcome), ", first sessions only"),
        a25f$n_natural_first + a25f$n_composite_first,
        fnum(a25f$mean_natural_first, 2), fnum(a25f$mean_composite_first, 2),
        fnum(a25f$estimate, 2), "--", fp(a25f$p_welch), "--"))
si_table(m20, "p{3.2cm}ccccccc",
  "Outcome & $n$ & Nat.-first & Comp.-first & Est. & SE & $p$ & $q$",
  paste0("First-block between-group analysis: between-group contrasts at each session's ",
         "first exposure, by the side the session opened with (mixed model with ",
         "participant random intercept; BH over 3 tests). The opening side was varied ",
         "across sessions but not by a recorded randomisation procedure. The last three ",
         "rows were added post-review: the same contrasts restricted to each ",
         "participant's first session (27 sessions; Welch $t$, uncorrected), removing ",
         "all cross-session prior exposure. Within-visit role order (whether the ",
         "measured participant had first served as the eyes-closed binaural recorder) ",
         "was not recorded, so acoustic pre-exposure within a visit cannot be excluded ",
         "for the paired partner."),
  "tab:s-naive", file.path(SDIR, "tabS20_naive.tex"))

# S21 time of day
# Withheld input: session clock-times are not in the deposit, so A24 does not
# run and this fragment cannot be regenerated; the printed version stands.
a24_file <- file.path("output", "a24_time_of_day", "time_of_day_tests.csv")
if (!file.exists(a24_file)) {
  message("tabS21_time_of_day.tex not regenerated: session clock-times are ",
          "withheld in the deposit (disclosure control).")
} else {
a24 <- read_csv(a24_file, show_col_types = FALSE)
m21 <- cbind(gsub("start_hour x dP_ISO", "start hour $\\\\times$ $\\\\Delta$ISO pleasantness",
             gsub("start_hour x dFrontal-HbO",
                  "start hour $\\\\times$ $\\\\Delta$frontal HbO", a24$outcome)),
             a24$n_pairs, fnum(a24$r, 2), fp(a24$p), fp(a24$q))
si_table(m21, "lcccc", "Test & $n$ & $r$ & $p$ & $q$",
  paste0("Time-of-day check on the orientation confound (BH over 2 tests): neither ",
         "main shift varies with session start hour across the sampled 08:00--17:00 ",
         "range."),
  "tab:s-timeofday", file.path(SDIR, "tabS21_time_of_day.tex"))
}

# S22 fragment — post-review fixation–HbO coupling (cited as Supplementary Table S21)
a26f <- OUTC("analysis_26_fixation_brain_coupling", "fixation_hbo_family.csv")
a26d <- OUTC("analysis_26_fixation_brain_coupling", "per_channel_descriptive.csv")
m22 <- rbind(
  cbind("Family",
        c(paste0("pair level: $\\Delta$frontal HbO on $\\Delta$fixation, ",
                 "controlling $\\Delta L_{\\mathrm{Aeq}}$"),
          "block level: within-side fixation slope on frontal HbO"),
        a26f$n, fnum(a26f$estimate, 2), fnum(a26f$se, 2), fp(a26f$p), fp(a26f$q)),
  cbind("Descriptive", L(a26d$channel), a26d$n_pairs, fnum(a26d$beta_std, 2),
        fnum(a26d$se, 2), fp(a26d$p_descriptive), "--"))
si_table(m22, "lp{4.1cm}ccccc",
  "Set & Test & $n$ & Estimate & SE & $p$ & $q$",
  paste0("Post-review exploratory fixation--HbO coupling (a two-test BH family declared ",
         "before the analysis ran). The pair-level row is the standardised slope of the ",
         "natural-minus-composite frontal HbO difference on the fixation difference, ",
         "controlling the $L_{\\mathrm{Aeq}}$ difference, with a participant random ",
         "intercept over the 30 eye-tracked samples; the block-level row is the ",
         "standardised within-side fixation slope on frontal HbO under the full mixed ",
         "model. Descriptive rows repeat the pair-level model per channel and for frontal ",
         "HbR. The pair-level association is rank-robust (Spearman $\\rho = -0.47$, ",
         "$p = 0.010$, holding after removal of the three largest $\\Delta$fixation pairs) ",
         "and unchanged under duplicate-vector exclusion, a first-side covariate, a site ",
         "random intercept and omission of the acoustic control; its parametric ",
         "leave-one-out $p$ reaches 0.19 when the most extreme pair is dropped, while the ",
         "rank-based estimate remains significant after the same drop ($p = 0.029$)."),
  "tab:s-fixation-hbo", file.path(SDIR, "tabS22_fixation_hbo.tex"))

message("All tables written.")
