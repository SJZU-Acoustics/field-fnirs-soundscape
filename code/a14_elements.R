# A14 — element-class screen (exhaustion-audit gap-fill: the 44-class layer
# below the summary indices had never been consumed).
# Family (declared): 16 tests = 8 aggregated class groups x 2 outcomes
# (dP_ISO, dfrontal-HbO), between-side level, BH as one family.
# Class order differs between sides in the delivered sheets; the frozen
# element table is long-format by name, so the join is by name (house trap noted).

source(file.path("code", "helpers.R"))
A_DIR <- analysis_output_dir("a14_element_class_screen")

blocks <- load_p27_blocks()
elem   <- load_p27_elements() %>%
  mutate(context_recorded = .to_logical(context_recorded)) %>%
  filter(context_recorded)

GROUPS <- list(
  water = c("水体", "河流", "海洋"),
  trees = c("树木", "棕榈树"),
  grass = c("草地"),
  other_vegetation = c("植物", "花卉"),
  sky = c("天空"),
  buildings = c("建筑", "房屋"),
  paving = c("人行道", "小路", "道路", "地板", "跑道"),
  people = c("行人")
)

g <- elem %>% select(sample_id, side, block)
for (nm in names(GROUPS)) {
  cols <- intersect(GROUPS[[nm]], names(elem))
  g[[nm]] <- rowSums(elem[, cols, drop = FALSE], na.rm = TRUE)
}
gsm <- g %>% group_by(sample_id, side) %>%
  summarise(across(all_of(names(GROUPS)), mean), .groups = "drop")

# outcomes at sample-side level
frontal_cols <- paste0("hbo_ch", 1:4)
out <- load_p27_samples() %>%
  mutate(hbo_frontal = rowMeans(across(all_of(frontal_cols)), na.rm = TRUE)) %>%
  select(sample_id, participant, side, p_iso, hbo_frontal)

dd <- left_join(out, gsm, by = c("sample_id", "side")) %>%
  pivot_wider(names_from = side,
              values_from = all_of(c("p_iso", "hbo_frontal", names(GROUPS))),
              names_glue = "{.value}__{side}")
for (v in c("p_iso", "hbo_frontal", names(GROUPS)))
  dd[[paste0("d_", v)]] <- dd[[paste0(v, "__natural")]] - dd[[paste0(v, "__composite")]]

res <- map_dfr(c("p_iso", "hbo_frontal"), function(y) map_dfr(names(GROUPS), function(x) {
  d <- dd %>% filter(is.finite(.data[[paste0("d_", y)]]), is.finite(.data[[paste0("d_", x)]]))
  if (sd(d[[paste0("d_", x)]]) == 0) return(tibble(outcome = y, class_group = x,
    n = nrow(d), r = NA_real_, p = NA_real_))
  ct <- cor.test(d[[paste0("d_", x)]], d[[paste0("d_", y)]])
  tibble(outcome = y, class_group = x, n = nrow(d), r = unname(ct$estimate), p = ct$p.value)
})) %>% mutate(q = bh(p), sig_q05 = q < 0.05) %>% arrange(q)
write_outcome(res, A_DIR, "element_class_screen.csv")

# descriptive: what actually differs between sides at class-group level
side_diff <- map_dfr(names(GROUPS), function(x) {
  x_ <- dd[[paste0("d_", x)]]; x_ <- x_[is.finite(x_)]
  tibble(class_group = x, mean_diff = mean(x_), d_z = mean(x_)/sd(x_),
         p = t.test(x_)$p.value)
})
write_outcome(side_diff, A_DIR, "class_group_side_differences.csv")

print(res, n = 20); print(side_diff)
cat("A14 done\n")
