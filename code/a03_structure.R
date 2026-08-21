# A03 — perception structure in situ.
# (i) Does the 8-item PAQ behave as a two-dimensional circumplex in the field?
# (ii) Are the delivered P_ISO/E_ISO adequate summaries (vs the first PCs)?
# (iii) Is the A02 side-shift profile across the 8 items consistent with a pure
#       pleasantness-axis rotation, item shift ~ a*cos(theta) + b*sin(theta)?
# Descriptive/structural analysis: no new inferential family beyond the stated
# profile fit (2 coefficients, native p).

source(file.path("code", "helpers.R"))
A_DIR <- analysis_output_dir("a03_perception_structure")

blocks  <- load_p27_blocks()
samples <- load_p27_samples()

# Chinese-corrected item angles (codebook section 2), degrees
ANGLES <- c(paq_pleasant = 0, paq_eventful = 18, paq_engaging = 38,
            paq_chaotic = 154, paq_annoying = 171, paq_monotonous = 196,
            paq_uneventful = 217, paq_calm = 318)

# ---- 1. Correlation structure at three levels ------------------------------
# between-sample (sample means), within-sample-between-side (side means minus
# sample mean), within-side (block deviations)
d <- blocks %>% select(sample_id, side, block, side_unit, all_of(PAQ_COLS))
bs <- d %>% group_by(sample_id) %>% summarise(across(all_of(PAQ_COLS), mean))
ws <- d %>% group_by(sample_id, side) %>% summarise(across(all_of(PAQ_COLS), mean), .groups = "drop") %>%
  group_by(sample_id) %>% mutate(across(all_of(PAQ_COLS), ~.x - mean(.x))) %>% ungroup()
wb <- d %>% group_by(side_unit) %>% mutate(across(all_of(PAQ_COLS), ~.x - mean(.x))) %>% ungroup()

corr_long <- function(m, level) {
  cc <- cor(m)
  as_tibble(cc, rownames = "item_i") %>%
    pivot_longer(-item_i, names_to = "item_j", values_to = "r") %>%
    filter(item_i < item_j) %>% mutate(level = level)
}
cors <- bind_rows(
  corr_long(bs %>% select(all_of(PAQ_COLS)), "between_sample"),
  corr_long(ws %>% select(all_of(PAQ_COLS)), "within_sample_between_side"),
  corr_long(wb %>% select(all_of(PAQ_COLS)), "within_side_block")
)
write_outcome(cors, A_DIR, "item_correlations_by_level.csv")

# ---- 2. PCA on sample-level side means; congruence with delivered axes -----
sm <- samples %>% select(sample_id, side, all_of(PAQ_COLS), p_iso, e_iso)
pca <- prcomp(sm %>% select(all_of(PAQ_COLS)), scale. = TRUE)
ev <- tibble(pc = seq_along(pca$sdev), var_share = pca$sdev^2 / sum(pca$sdev^2))
write_outcome(ev, A_DIR, "pca_variance.csv")
loadings <- as_tibble(pca$rotation[, 1:3], rownames = "item") %>%
  mutate(angle_deg = ANGLES[item])
write_outcome(loadings, A_DIR, "pca_loadings.csv")

scores <- bind_cols(sm %>% select(sample_id, side, p_iso, e_iso),
                    as_tibble(pca$x[, 1:2]))
cong <- tibble(
  pair = c("p_iso~PC1", "p_iso~PC2", "e_iso~PC1", "e_iso~PC2"),
  r = c(cor(scores$p_iso, scores$PC1), cor(scores$p_iso, scores$PC2),
        cor(scores$e_iso, scores$PC1), cor(scores$e_iso, scores$PC2))
)
# best linear recovery of each ISO axis from the two PCs
cong <- bind_rows(cong,
  tibble(pair = c("p_iso~PC1+PC2 (R)", "e_iso~PC1+PC2 (R)"),
         r = c(sqrt(summary(lm(p_iso ~ PC1 + PC2, scores))$r.squared),
               sqrt(summary(lm(e_iso ~ PC1 + PC2, scores))$r.squared))))
write_outcome(cong, A_DIR, "iso_axis_congruence.csv")

# ---- 3. The side-shift profile against the circumplex ----------------------
paired <- load_p27_paired(PAQ_COLS)
shift <- map_dfr(PAQ_COLS, function(v) {
  x <- paired[[paste0(v, "__diff")]]
  tibble(item = v, mean_shift = mean(x), se = sd(x)/sqrt(length(x)),
         d_z = mean(x)/sd(x), angle_deg = ANGLES[v])
}) %>%
  mutate(cos_t = cos(angle_deg * pi/180), sin_t = sin(angle_deg * pi/180))
# weighted least squares of shift on (cos, sin): a = P-axis displacement,
# b = E-axis displacement implied by the item profile
fit <- lm(mean_shift ~ 0 + cos_t + sin_t, data = shift, weights = 1/shift$se^2)
co <- summary(fit)$coefficients
profile_fit <- tibble(
  term = rownames(co), est = co[, 1], se = co[, 2], t = co[, 3], p = co[, 4],
  r2 = summary(fit)$r.squared)
write_outcome(profile_fit, A_DIR, "shift_profile_fit.csv")
shift <- shift %>% mutate(fitted = predict(fit), resid = mean_shift - fitted)
write_outcome(shift, A_DIR, "item_shift_profile.csv")

# internal consistency of the two delivered axes (Cronbach alpha on the
# pole items, reversed): P: pleasant, calm vs chaotic, annoying;
# E: eventful, engaging vs monotonous, uneventful (per the ISO angles)
alpha_of <- function(m) {
  k <- ncol(m); cm <- cov(m)
  (k/(k-1)) * (1 - sum(diag(cm))/sum(cm))
}
pa <- sm %>% transmute(paq_pleasant, paq_calm, chaotic_r = 6 - paq_chaotic, annoying_r = 6 - paq_annoying)
ea <- sm %>% transmute(paq_eventful, paq_engaging, monotonous_r = 6 - paq_monotonous, uneventful_r = 6 - paq_uneventful)
write_outcome(tibble(axis = c("P_items", "E_items"),
                     alpha = c(alpha_of(pa), alpha_of(ea))), A_DIR, "axis_alpha.csv")

print(ev, n = 3); print(cong); print(profile_fit); print(shift %>% select(item, mean_shift, fitted, resid, d_z))
cat("A03 done\n")
