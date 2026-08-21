#!/usr/bin/env Rscript
# =============================================================================
# run_all.R — reproduce every analysis, figure and table of the article from
# the deposited workbook.
#
#   Rscript run_all.R        (from the repository root)
#
# Data input:  data/In_situ_riverside_park_soundscape_fNIRS_eyetracking_data.xlsx
#              (Mendeley Data, CC BY 4.0 — see README.md; hash-verified on read)
# Output:      output/aNN_*/           analysis outcome CSVs (one per module)
#              output/figures/         fig1..fig4 PNG (600 dpi) + data_lock CSVs
#              output/tables/          main-text LaTeX table fragments
#              output/supplementary/   SI LaTeX table environments
#
# Two items are deliberately not regenerated (withheld inputs, disclosure
# control — see README.md): the two age-moderation rows of Supplementary
# Table S17 and Supplementary Table S22 (time of day). The affected modules
# (A15, A24) and table fragments degrade gracefully with a printed message.
# =============================================================================

if (!file.exists(file.path("code", "helpers.R"))) {
  stop("Run from the repository root: Rscript run_all.R", call. = FALSE)
}

t0 <- Sys.time()

# ---- analysis modules (a01 .. a26b, alphabetical) ---------------------------
modules <- sort(list.files("code", pattern = "^a[0-9].*[.]R$", full.names = TRUE))
for (m in modules) {
  message("== ", m)
  source(m)
}

# ---- display layer ----------------------------------------------------------
for (s in c("fig1_design.R", "fig2_perception.R", "fig3_physiology.R",
            "fig4_generality.R", "make_tables.R")) {
  message("== code/", s)
  source(file.path("code", s))
}

message(sprintf("run_all finished in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
