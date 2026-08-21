# A20 — site-level structure: is the site texture (A12 heterogeneity,
# p = 0.002) SPATIALLY organised along the linear park, and do the two
# response systems co-vary at the site grain (the one coupling level A07 did
# not cover)?
# Family (declared): 3 tests, BH — (i) spatial autocorrelation of site-mean
# dP_ISO: Mantel-style cor of pairwise effect differences with pairwise
# geographic distances, 1,000 site-label permutations, one-sided (similarity
# decays with distance => negative r); (ii) same for dFrontal-HbO;
# (iii) cross-outcome coupling: cor(site-mean dP_ISO, site-mean dFrontal-HbO).

source(file.path("code", "helpers.R"))
A_DIR <- analysis_output_dir("a20_site_structure")

samples <- load_p27_samples()
sites   <- load_p27_sites()

pw <- samples %>%
  rowwise() %>%
  mutate(hbo_frontal = mean(c_across(all_of(paste0("hbo_ch", 1:4))), na.rm = TRUE)) %>%
  ungroup() %>%
  select(sample_id, site, side, p_iso, hbo_frontal) %>%
  pivot_wider(names_from = side, values_from = c(p_iso, hbo_frontal)) %>%
  mutate(d_piso = p_iso_natural - p_iso_composite,
         d_hbof = hbo_frontal_natural - hbo_frontal_composite)

site_eff <- pw %>% group_by(site) %>%
  summarise(d_piso = mean(d_piso, na.rm = TRUE),
            d_hbof = mean(d_hbof, na.rm = TRUE), n = n(), .groups = "drop") %>%
  left_join(sites %>% mutate(site = factor(site)) %>%
              select(site, longitude, latitude), by = "site") %>%
  arrange(site)
write_outcome(site_eff, A_DIR, "site_effects_with_coordinates.csv")

# geographic distance matrix (equirectangular approximation, metres)
lon <- site_eff$longitude; lat <- site_eff$latitude
lat0 <- mean(lat) * pi / 180
dx <- outer(lon, lon, `-`) * 111320 * cos(lat0)
dy <- outer(lat, lat, `-`) * 111320
geo <- sqrt(dx^2 + dy^2)

mantel_perm <- function(eff) {
  ok <- is.finite(eff)
  e <- eff[ok]; g <- geo[ok, ok]
  iu <- upper.tri(g)
  obs <- cor(abs(outer(e, e, `-`))[iu], g[iu])
  set.seed(271)
  nulls <- replicate(BOOT_N, {
    ep <- sample(e)
    cor(abs(outer(ep, ep, `-`))[iu], g[iu])
  })
  tibble(r_obs = obs, null_median = median(nulls),
         null_q05 = quantile(nulls, .05),   # one-sided: observed MORE negative
         p_perm = mean(nulls <= obs))
}

m_piso <- mantel_perm(site_eff$d_piso) %>% mutate(test = "spatial: dP_ISO")
m_hbof <- mantel_perm(site_eff$d_hbof) %>% mutate(test = "spatial: dFrontal-HbO")

ct <- cor.test(site_eff$d_piso, site_eff$d_hbof)
xout <- tibble(test = "cross-outcome coupling: dP_ISO x dFrontal-HbO",
               r_obs = unname(ct$estimate), null_median = NA_real_,
               null_q05 = NA_real_, p_perm = ct$p.value)

res <- bind_rows(m_piso, m_hbof, xout) %>%
  mutate(n_sites = nrow(site_eff), q = bh(p_perm)) %>%
  relocate(test, n_sites)
write_outcome(res, A_DIR, "site_structure_tests.csv")

print(res)
cat("A20 done\n")
