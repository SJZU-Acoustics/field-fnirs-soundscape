# A16 — time course of the side effect within a side: does the natural-side
# benefit appear from the first 60 s, and does it build or fade across the
# three repeated exposures? Never tested by A11 (which ran session-position
# main-effect trends and first_side moderation, not side x block).
# Family (declared): 4 tests = side x block_c interaction (block centred at 2)
# for P_ISO, E_ISO, monotony, frontal HbO; BH. Per-block profile descriptive.

source(file.path("code", "helpers.R"))
A_DIR <- analysis_output_dir("a16_effect_time_course")

blocks <- load_p27_blocks() %>% mutate(block_c = block - 2L)

# block-level frontal HbO: mean over VALID frontal channels, recorded blocks only
hbo_long <- blocks %>%
  filter(fnirs_recorded %in% TRUE) %>%
  select(sample_id, participant, site, side, block, block_c, side_unit,
         starts_with("hbo_ch") & matches("hbo_ch[1-4]($|_)"))
frontal_blk <- hbo_long %>%
  rowwise() %>%
  mutate(hbo_frontal = mean(c_across(all_of(paste0("hbo_ch", 1:4)))[
    c_across(all_of(paste0("hbo_ch", 1:4, "_valid")))], na.rm = TRUE)) %>%
  ungroup() %>%
  filter(is.finite(hbo_frontal)) %>%
  select(sample_id, side, block, block_c, side_unit, participant, site, hbo_frontal)

RE4 <- c("(1 | participant)", "(1 | site)", "(1 | sample_id)", "(1 | side_unit)")

fit_one <- function(d, y) {
  f <- as.formula(paste0(y, " ~ side * block_c + ",
                         paste(RE4, collapse = " + ")))
  fit <- suppressMessages(suppressWarnings(lmer(f, data = d, REML = TRUE)))
  co <- summary(fit)$coefficients
  ix <- rownames(co)
  take <- function(nm) if (nm %in% ix) unlist(co[nm, c("Estimate", "Std. Error", "Pr(>|t|)")]) else c(NA, NA, NA)
  main <- take("sidenatural"); inter <- take("sidenatural:block_c")
  tibble(outcome = y, n_rows = nrow(d), singular = isSingular(fit),
         side_est = main[1], side_se = main[2], side_p = main[3],
         inter_est = inter[1], inter_se = inter[2], inter_p = inter[3])
}

res <- bind_rows(
  fit_one(blocks, "p_iso"),
  fit_one(blocks, "e_iso"),
  fit_one(blocks, "paq_monotonous"),
  frontal_blk %>% rename(hbo_frontal_y = hbo_frontal) %>%
    { fit_one(rename(., hbo_frontal = hbo_frontal_y), "hbo_frontal") }
) %>% mutate(q_inter = bh(inter_p))
write_outcome(res, A_DIR, "side_x_block_interactions.csv")

# descriptive: per-block side-difference profile (natural - composite)
prof <- bind_rows(
  blocks %>% group_by(side, block) %>%
    summarise(p_iso = mean(p_iso), e_iso = mean(e_iso),
              monotony = mean(paq_monotonous), .groups = "drop") %>%
    pivot_wider(names_from = side, values_from = c(p_iso, e_iso, monotony)) %>%
    mutate(d_piso = p_iso_natural - p_iso_composite,
           d_eiso = e_iso_natural - e_iso_composite,
           d_mono = monotony_natural - monotony_composite) %>%
    select(block, d_piso, d_eiso, d_mono),
  frontal_blk %>% group_by(side, block) %>%
    summarise(hbof = mean(hbo_frontal), .groups = "drop") %>%
    pivot_wider(names_from = side, values_from = hbof) %>%
    mutate(d_hbof = natural - composite) %>% select(block, d_hbof)
) %>% arrange(block) %>%
  group_by(block) %>% summarise(across(everything(), ~ first(.x[!is.na(.x)])), .groups = "drop")
write_outcome(prof, A_DIR, "per_block_side_difference_profile.csv")

print(res); print(prof)
cat("A16 done\n")
