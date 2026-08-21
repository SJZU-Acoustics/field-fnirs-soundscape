# A24 — time-of-day bounding of the orientation confound. Natural sides face
# north (river), composite south, so sun/thermal asymmetry cannot be separated
# from scene category within a site — but session time varies sun angle across
# the day, and an effect stable across times of day is harder to read as glare
# or thermal. Family (declared): 2 tests — session start hour x dP_ISO and
# x dFrontal-HbO (pair level), BH. By-slot table descriptive.

source(file.path("code", "helpers.R"))
A_DIR <- analysis_output_dir("a24_time_of_day")

samples <- load_p27_samples()

# Disclosure control: the deposit omits session clock-times, so session start
# hour cannot be derived and this analysis cannot be recomputed from it.
# Skipped with a message; SI Table S22 stands as printed in the article.
if (!"session_time" %in% names(samples)) {
  message("A24 skipped: session clock-times are withheld in the deposit ",
          "(disclosure control), so start hour is not derivable. ",
          "SI Table S22 is not regenerated.")
} else {

pw <- samples %>%
  rowwise() %>%
  mutate(hbo_frontal = mean(c_across(all_of(paste0("hbo_ch", 1:4))), na.rm = TRUE)) %>%
  ungroup() %>%
  select(sample_id, participant, site, side, p_iso, hbo_frontal, session_time) %>%
  pivot_wider(names_from = side, values_from = c(p_iso, hbo_frontal),
              names_sep = "__") %>%
  mutate(d_piso = p_iso__natural - p_iso__composite,
         d_hbof = hbo_frontal__natural - hbo_frontal__composite,
         start_hour = as.numeric(sub(":.*", "", session_time)))

ct_row <- function(y, label, d) {
  dd <- d %>% filter(is.finite(.data[[y]]), is.finite(start_hour))
  ct <- cor.test(dd$start_hour, dd[[y]])
  tibble(outcome = label, n_pairs = nrow(dd),
         r = unname(ct$estimate), p = ct$p.value)
}

res <- bind_rows(ct_row("d_piso", "start_hour x dP_ISO", pw),
                 ct_row("d_hbof", "start_hour x dFrontal-HbO", pw)) %>%
  mutate(q = bh(p))
write_outcome(res, A_DIR, "time_of_day_tests.csv")

slot <- pw %>% group_by(session_time) %>%
  summarise(n = n(), start_hour = first(start_hour),
            d_piso_mean = mean(d_piso, na.rm = TRUE),
            d_hbof_mean = mean(d_hbof, na.rm = TRUE), .groups = "drop") %>%
  arrange(start_hour)
write_outcome(slot, A_DIR, "by_session_slot_descriptive.csv")

print(res); print(slot)
cat("A24 done\n")

}
