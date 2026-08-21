# =============================================================================
# load_data.R — single data entry for the field-fnirs-soundscape release.
#
# Data source: the Mendeley Data workbook
#   data/In_situ_riverside_park_soundscape_fNIRS_eyetracking_data.xlsx  (CC BY 4.0)
#
# The workbook's SHA-256 is verified on first read: if the file is not the
# deposited copy, the run stops rather than analysing the wrong data.
#
# Sheets are read by a small direct XLSX reader below — NOT via readxl. Both
# readxl's text and numeric conversions round cell values with only 16
# significant digits, one short of IEEE-754 round-trip, which lost the last
# ulp against the frozen working tables (verification 2026-08-21). This
# reader takes the cell's own decimal string from the sheet XML and parses it
# with base R (correctly rounded), so every value is bit-identical to the
# working pipeline's frozen CSVs. Types follow the working pipeline's
# conventions: booleans -> logical, numbers -> numeric, text -> character.
#
# Every analysis module reads its input through read_frozen_csv(); the mapping
# below resolves each original processed-CSV name to its workbook sheet.
# Nothing else touches the workbook.
#
# WITHHELD INPUTS (disclosure control, deliberate): exact participant age
# (the workbook carries 5-year bands) and session clock-times (the workbook
# carries dates). Two analysis modules therefore degrade gracefully:
#   - A15 person level: the two age-moderation tests are skipped (sex rows run);
#   - A24 time of day: skipped entirely (session start hour is not derivable).
# The display items that depend on them (Supplementary Tables S17 age rows and
# S22) are named in the repository README.
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(tibble)
})

XLSX_PATH <- file.path("data", "In_situ_riverside_park_soundscape_fNIRS_eyetracking_data.xlsx")

# SHA-256 of the deposited workbook (DEPOSIT_MANIFEST.json, publication copy
# built 2026-08-21). Mendeley Data serves the uploaded file byte-identically.
EXPECTED_XLSX_SHA256 <- "24acdbc465edcf96f9892571bb93de6664d45b7adc52cb922c198edf85f6e0fa"

# processed-CSV name -> workbook sheet
SHEET_FOR_FILE <- c(
  "p27_block_level.csv"            = "block_level",
  "p27_sample_level.csv"           = "sample_level",
  "p27_element_proportions.csv"    = "element_proportions",
  "p27_fnirs_channel_validity.csv" = "fnirs_channel_validity",
  "p27_participants.csv"           = "participants",
  "p27_sites.csv"                  = "sites",
  "p27_pairing.csv"                = "pairing"
)

.verify_workbook <- function() {
  if (!file.exists(XLSX_PATH)) {
    stop("Workbook not found at ", XLSX_PATH, "\n",
         "Download it from Mendeley Data and place it in data/ — see README.",
         call. = FALSE)
  }
  got <- tools::sha256sum(XLSX_PATH)
  if (!identical(unname(got), EXPECTED_XLSX_SHA256)) {
    stop("Workbook hash mismatch: this is not the deposited copy.\n",
         "  expected ", EXPECTED_XLSX_SHA256, "\n  got      ", got,
         call. = FALSE)
  }
  invisible(TRUE)
}

# ------------------------------------------------------------ direct XLSX read
.col_index <- function(letters) {
  n <- 0L
  for (ch in strsplit(letters, "", fixed = TRUE)[[1]]) {
    n <- n * 26L + (match(ch, LETTERS))
  }
  n
}

.read_xlsx_sheet <- function(xlsx, sheet_name) {
  if (!requireNamespace("xml2", quietly = TRUE)) {
    stop("Package 'xml2' is required (it ships with the tidyverse).", call. = FALSE)
  }
  tmp <- tempfile(pattern = "p27xlsx_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  utils::unzip(xlsx, exdir = tmp)

  wbx   <- xml2::read_xml(file.path(tmp, "xl", "workbook.xml"))
  rels  <- xml2::read_xml(file.path(tmp, "xl", "_rels", "workbook.xml.rels"))
  sh    <- xml2::xml_find_all(wbx, "//*[local-name()='sheet']")
  rids  <- xml2::xml_attr(sh, "id")          # r:id shows up as plain "id"
  names(rids) <- xml2::xml_attr(sh, "name")
  if (!sheet_name %in% names(rids)) {
    stop("Sheet '", sheet_name, "' not found in the workbook.", call. = FALSE)
  }
  rel_nodes <- xml2::xml_find_all(rels, "//*[local-name()='Relationship']")
  targets <- xml2::xml_attr(rel_nodes, "Target")
  names(targets) <- xml2::xml_attr(rel_nodes, "Id")
  tgt <- targets[[rids[[sheet_name]]]]
  sheet_file <- if (startsWith(tgt, "/")) file.path(tmp, sub("^/", "", tgt)) else file.path(tmp, "xl", tgt)

  shared <- character()
  ss_path <- file.path(tmp, "xl", "sharedStrings.xml")
  if (file.exists(ss_path)) {
    ss <- xml2::read_xml(ss_path)
    shared <- vapply(xml2::xml_find_all(ss, "//*[local-name()='si']"),
                     function(si) paste(xml2::xml_text(xml2::xml_find_all(si, ".//*[local-name()='t']")),
                                        collapse = ""),
                     character(1))
  }

  sx <- xml2::read_xml(sheet_file)
  rows <- xml2::xml_find_all(sx, "//*[local-name()='sheetData']/*[local-name()='row']")
  cells <- xml2::xml_find_all(rows, "./*[local-name()='c']")
  refs  <- xml2::xml_attr(cells, "r")
  types <- xml2::xml_attr(cells, "t")
  # numeric/boolean/shared-string cells carry <v>; openpyxl writes text as
  # inline strings (<is><t>), so take whichever is present
  vtxt  <- vapply(xml2::xml_find_first(cells, "./*[local-name()='v']"),
                  function(v) if (inherits(v, "xml_missing")) NA_character_ else xml2::xml_text(v),
                  character(1))
  itxt  <- vapply(xml2::xml_find_first(cells, "./*[local-name()='is']/*[local-name()='t']"),
                  function(v) if (inherits(v, "xml_missing")) NA_character_ else xml2::xml_text(v),
                  character(1))
  vals  <- ifelse(is.na(vtxt), itxt, vtxt)

  row_no  <- as.integer(sub("^[A-Z]+", "", refs))
  col_no  <- vapply(sub("[0-9]+$", "", refs), .col_index, integer(1))
  hrow    <- row_no == 1L
  header  <- vals[hrow]
  hcols   <- col_no[hrow]
  htypes  <- types[hrow]
  htypes[is.na(htypes)] <- "n"
  resolved <- ifelse(htypes == "s", shared[as.integer(header) + 1L], header)
  ncols <- max(hcols)
  col_names <- character(ncols)
  col_names[hcols] <- resolved

  out <- vector("list", ncols)
  names(out) <- col_names
  for (j in seq_len(ncols)) out[[j]] <- rep(NA, length(rows) - 1L)

  body <- row_no > 1L
  bcol <- col_no[body]; brow <- row_no[body] - 1L
  btyp <- types[body]; bval <- vals[body]
  btyp[is.na(btyp)] <- "n"
  # `[[<-` on an atomic vector coerces the vector to the assigned type, so a
  # homogeneous column ends up logical / numeric / character exactly as the
  # working pipeline's read_csv guessed it.
  for (k in which(!is.na(bval))) {
    out[[bcol[k]]][[brow[k]]] <- switch(btyp[k],
      s = shared[as.integer(bval[k]) + 1L],
      b = bval[k] == "1",
      n = readr::parse_double(bval[k]),  # same parser as the working pipeline
      bval[k])
  }
  tibble::as_tibble(out)
}

# Cache of parsed sheets within the session (the workbook is verified once).
.sheet_cache <- new.env(parent = emptyenv())

read_workbook_sheet <- function(sheet) {
  if (exists(sheet, envir = .sheet_cache)) return(get(sheet, envir = .sheet_cache))
  .verify_workbook()
  df <- .read_xlsx_sheet(XLSX_PATH, sheet)
  assign(sheet, df, envir = .sheet_cache)
  df
}

#' Read a canonical processed table by its original file name.
read_frozen_csv <- function(filename, ...) {
  sheet <- SHEET_FOR_FILE[[filename]]
  if (is.null(sheet)) {
    stop("No workbook sheet is registered for ", filename, call. = FALSE)
  }
  read_workbook_sheet(sheet)
}
