# field-fnirs-soundscape

R code for reproducing the statistical analyses, figures and tables for the manuscript "Visual natural context in an urban park improves soundscape appraisal and lowers frontal oxygenation in situ".

## Requirements

- R 4.5+ (developed and verified on R 4.5.3)
- CRAN packages: `tidyverse`, `xml2`, `lme4`, `lmerTest`, `broom`, `broom.mixed`, `patchwork`, `ggrepel`, `ggtext`, `ragg`
- Install:

  ```r
  install.packages(c("tidyverse", "xml2", "lme4", "lmerTest", "broom",
                     "broom.mixed", "patchwork", "ggrepel", "ggtext", "ragg"))
  ```

- No non-standard hardware is required. A complete run takes a few minutes on a normal desktop; the permutation and bootstrap modules account for most of it.

## Data

The analysis reads a single deposited input: the Mendeley Data workbook
`In_situ_riverside_park_soundscape_fNIRS_eyetracking_data.xlsx` (CC BY 4.0).

1. Download it from Mendeley Data — see the manuscript's data-availability statement for the DOI.
2. Place the `.xlsx` file in the `data/` folder (see `data/README.md` for the sheet map).

The workbook holds the analysis-ready record of the field study: block-level PAQ ratings, ISO coordinates, source salience, in-situ acoustics, viewed-scene annotations, eye-tracking fixation durations and per-channel fNIRS HbO/HbR with the per-layer recorded flags, the channel-validity table, participants, sites and the pairing design.

Two small Figure 1 inputs ship with this repository in `data/figure_inputs/` because they are transcriptions or third-party geometry rather than measurements: the OpenStreetMap base geometry behind the site map (© OpenStreetMap contributors, ODbL) and the digitised fNIRS montage layout. See `data/figure_inputs/README.md`.

### Withheld inputs and the two items that are not regenerated

The deposit applies disclosure control: participant age is released in 5-year bands, and session clock-times are dropped (dates remain). Two display items therefore cannot be recomputed from the deposit and **stand as printed in the article**:

- **Supplementary Table S17** — its two age-moderation rows (module A15 skips them with a message; the sex-moderation and reliability rows still run).
- **Supplementary Table S22** — time-of-day bounding of the orientation confound (module A24 needs session start hour and is skipped with a message).

Every other figure, table and analysis regenerates from the deposit alone.

## File structure

- `run_all.R` — master script: runs every analysis module, then every display script.
- `code/load_data.R` — single data entry: verifies the workbook's SHA-256 and maps each canonical table to its workbook sheet. Sheets are read by a small direct XLSX reader (not readxl, whose cell conversion loses the last ulp against the frozen working tables), so every value is bit-identical to the working pipeline.
- `code/helpers.R` — the loaders, the layer-filter rule (a block not recorded was filled per layer; analyses filter on the flag of the layer they use), the paired side-contrast estimator, the mixed-model wrapper, the between/within decomposition and the bootstrap.
- `code/style.R` — shared publication plot styling.
- `code/a*.R` — the 27 analysis modules (A01–A26, with A26b the robustness appendix): inventory and reliability (A01), side contrasts (A02), perception structure (A03), drivers (A04), fNIRS signal anatomy (A05–A06), brain–perception and brain–context coupling (A07–A08), multivariate pattern (A09), fixation strand (A10), order and time (A11), site generality (A12), sensitivity battery (A13), element classes (A14), person level (A15), time course (A16), frontal lateralisation (A17), cross-modal moderation (A18), scene-content moderators (A19, A22), site structure (A20), naive first exposure (A21), perception decoding (A23), time of day (A24, skipped — see above), post-review sensitivities (A25), fixation–brain coupling (A26, A26b).
- `code/fig1_design.R` … `code/fig4_generality.R` — Figures 1–4; `code/make_tables.R` — Table 1, Table 2 and the Supplementary Tables.

## Usage

From the repository root:

```bash
Rscript run_all.R
```

Outputs are written to:

- `output/figures/` — Figures 1–4 (PNG, 600 dpi, 178 mm) and `output/figures/data_lock/`, the plotted values behind every figure panel (CSV)
- `output/tables/` — main-text LaTeX fragments; `output/supplementary/` — the SI table environments
- `output/aNN_*/` — the full regenerated result tables of each analysis module

To keep the console log alongside the outputs:

```bash
Rscript run_all.R 2>&1 | tee output/run_log.txt
```

`output/` is produced at run time and is safe to delete.

## Verification

Run against the deposited workbook, this pipeline reproduces the manuscript's display items: Figures 1, 2 and 4 and all 22 regenerable LaTeX table fragments (2 main-text, 20 supplementary) are **byte-identical** to the paper's, and every figure's underlying data-lock CSV is byte-identical. Figure 3 differs only in panel c: its points are jittered, the working figure's jitter was an unseeded draw, and the release seeds it (`set.seed(271)`) so the render reproduces exactly run to run — the points sit at different positions within their group strips, nothing else; the figure's panels a and b are pixel-identical and its data-lock CSVs are byte-identical. All bootstraps and permutations are seeded, so interval endpoints reproduce exactly.

The two withheld-input items named above (Supplementary Table S17's age rows; Supplementary Table S22) are not regenerated.

## Notes

- The design is fully paired: 40 paired sessions, each facing the natural (river) and composite (city) sides of one of 25 riverside sites in alternation, 19 natural-first / 21 composite-first. 13 participants each contributed two sessions at different sites.
- The fNIRS montage is bilateral frontal plus left-only temporal and right-only occipital after the global channel exclusion, so no temporal/occipital laterality contrast exists in the dataset.
- Two cross-sample identical 3-block fixation vectors are handled per the manuscript: the later member of each pair is excluded where the duplication matters.

## License

Code in this repository is released under the MIT License (see `LICENSE`). The input data are archived separately under CC BY 4.0 at Mendeley Data. The OpenStreetMap geometry in `data/figure_inputs/` is © OpenStreetMap contributors, ODbL.
