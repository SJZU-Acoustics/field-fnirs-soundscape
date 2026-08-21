# Data folder

## The Mendeley Data workbook (required)

Download `In_situ_riverside_park_soundscape_fNIRS_eyetracking_data.xlsx` from
Mendeley Data (CC BY 4.0; the DOI is given in the manuscript's data-availability
statement) and place it in **this folder**:

```
data/In_situ_riverside_park_soundscape_fNIRS_eyetracking_data.xlsx
```

`run_all.R` reads nothing else from the deposit. The workbook is not tracked in
this repository — it is archived at Mendeley Data, which is its citable home.
The pipeline verifies the workbook's SHA-256 on first read, so a substituted or
corrupted file stops the run rather than silently analysing the wrong data.

Sheets used by the analysis (a codebook is the workbook's own `dictionary` and
`README` sheets):

| Sheet | What it holds |
|---|---|
| `block_level` | 240 rows = 40 samples × 2 sides × 3 blocks: PAQ items, ISO pleasantness/eventfulness, source salience, acoustics (LAeq, L5, L95), viewed-scene complexity/naturalness/artificiality, fixation durations, per-channel HbO/HbR, and the per-layer recorded flags |
| `sample_level` | 80 rows: each layer meaned over its own recorded blocks only; per-channel HbO/HbR cells for invalid channels are NA |
| `element_proportions` | semantic-segmentation class proportions per sample × side (natural scene) |
| `fnirs_channel_validity` | per sample × side × channel signal-quality calls (23 sample×side×channel units treated as missing) |
| `participants` | 27 participants: sex, 5-year **age bands** (exact age is withheld for disclosure control) |
| `sites` | the 25 riverside sites: GCJ-02 field coordinates, WGS-84 conversion, datum shift, distance to main road, faced-side bearings |
| `pairing` | the 40 recorder pairs: session **date** (clock-times withheld), opening side, recorder-analysed flag |

Two disclosure-control omissions are deliberate: exact participant age (banded
instead) and session clock-times (dates only). They make two supplementary items
not regenerable from the deposit — see the repository README.

## `figure_inputs/` (shipped here)

Two figure inputs that are transcriptions/third-party geometry rather than
measurements, so they travel with the code — see `figure_inputs/README.md` for
provenance and licence notes.
