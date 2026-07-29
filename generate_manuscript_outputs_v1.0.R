cat("Loading packages...\n")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lubridate)
  library(viridis)
  library(scales)
  library(gridExtra)
  library(grid)
  library(ragg)
  library(maps)
  library(broom)
})

# Figures containing the Unicode omega (ω) in axis labels are rendered with the
# ragg AGG device; the base R png() device fails on U+03C9 ("conversion failure
# on 'ω'"). save_grob() writes any assembled grob to a PNG via ragg.
save_grob <- function(path, grob, width, height) {
  agg_png(path, width = width, height = height, units = "in", res = 300)
  grid.newpage(); grid.draw(grob); dev.off()
}

setwd("/Users/hagen/Desktop/neonSoilHealth/projects/descriptor/manuscript")
dir.create("manuscript_outputs", showWarnings = FALSE)
dir.create("manuscript_outputs/tables", showWarnings = FALSE)
dir.create("manuscript_outputs/figures", showWarnings = FALSE)
dir.create("manuscript_outputs/data", showWarnings = FALSE)

cat("Output directories created.\n")

cat("Loading dataset...\n")
neon_plfa <- read.csv("../../../archive/superseded_data/neon_soil_health_2026-05-18.csv", stringsAsFactors = FALSE)
neon_plfa$collectDate <- as.Date(neon_plfa$collectDate)

cat("Loaded", nrow(neon_plfa), "samples from", length(unique(neon_plfa$siteID)), "sites\n")

cat("\n=== Enriching with pH, sample depth, and biomass condition (DP1.10086.001) ===\n")

# Reuses the already-downloaded DP1.10086.001 cache from the mane project (same NEON
# release, same 47 sites) rather than re-pulling from NEON. sls_soilpH and
# sls_soilCoreCollection share the PLFA sampleID scheme (see match-rate check), unlike
# sls_soilChemistry, which only overlaps ~16% of PLFA sampleIDs.
soil_periodic_cache <- readRDS("../../mane/cache/soilperiodic_DP1_10086.rds")

ph_add <- soil_periodic_cache$sls_soilpH %>%
  select(sampleID, soilInWaterpH, soilInCaClpH) %>%
  distinct(sampleID, .keep_all = TRUE)

core_add <- soil_periodic_cache$sls_soilCoreCollection %>%
  select(sampleID, sampleTopDepth, sampleBottomDepth, biomassSampleCondition) %>%
  distinct(sampleID, .keep_all = TRUE)

neon_plfa <- neon_plfa %>%
  left_join(ph_add, by = "sampleID") %>%
  left_join(core_add, by = "sampleID")

cat("soilInWaterpH non-NA:", sum(!is.na(neon_plfa$soilInWaterpH)), "/", nrow(neon_plfa), "\n")
cat("sampleTopDepth non-NA:", sum(!is.na(neon_plfa$sampleTopDepth)), "/", nrow(neon_plfa), "\n")
cat("biomassSampleCondition non-NA:", sum(!is.na(neon_plfa$biomassSampleCondition)), "/", nrow(neon_plfa), "\n")

cat("\n=== Calculating Gram-positive:Gram-negative ratio ===\n")

neon_plfa$gram_positive <- rowSums(neon_plfa[, c(
  "i14To0ScaledConcentration",
  "i15To0ScaledConcentration",
  "aC15To0ScaledConcentration",
  "i16To0ScaledConcentration",
  "i17To0ScaledConcentration",
  "c17To0AnteisoScaledConcentration"
)], na.rm = TRUE)

neon_plfa$gram_negative <- rowSums(neon_plfa[, c(
  "c16To1n7ScaledConcentration",
  "c18To1n11ScaledConcentration"
)], na.rm = TRUE)

neon_plfa$GP_GN_ratio <- ifelse(neon_plfa$gram_negative > 0,
                                 neon_plfa$gram_positive / neon_plfa$gram_negative,
                                 NA)

cat("GP:GN ratio calculated for", sum(!is.na(neon_plfa$GP_GN_ratio)), "samples\n")

cat("\n=== Calculating cyclopropyl:precursor stress indices ===\n")

# Physiological stress indices based on the accumulation of cyclopropyl fatty
# acids relative to their monoenoic precursors (higher = greater stress).
# NEON column mapping: cyclo19To0 = cy19:0; cyclo17To0 = cy17:0;
# c18To1n11 = 18:1w7c (cis-vaccenic); c16To1n7 = 16:1w7c.
# Primary index (cy19:0/18:1w7c) spans the full record; the cy17:0-based and
# combined indices depend on cyclo17To0, which NEON reports from ~2021 onward,
# so they are populated only where that precursor pair is measured.

neon_plfa$stress_index_cy19 <- ifelse(
  !is.na(neon_plfa$c18To1n11ScaledConcentration) &
    neon_plfa$c18To1n11ScaledConcentration > 0,
  neon_plfa$cyclo19To0ScaledConcentration / neon_plfa$c18To1n11ScaledConcentration,
  NA)

neon_plfa$stress_index_cy17 <- ifelse(
  !is.na(neon_plfa$cyclo17To0ScaledConcentration) &
    !is.na(neon_plfa$c16To1n7ScaledConcentration) &
    neon_plfa$c16To1n7ScaledConcentration > 0,
  neon_plfa$cyclo17To0ScaledConcentration / neon_plfa$c16To1n7ScaledConcentration,
  NA)

stress_num <- neon_plfa$cyclo17To0ScaledConcentration + neon_plfa$cyclo19To0ScaledConcentration
stress_den <- neon_plfa$c16To1n7ScaledConcentration + neon_plfa$c18To1n11ScaledConcentration
neon_plfa$stress_index_combined <- ifelse(
  !is.na(neon_plfa$cyclo17To0ScaledConcentration) &
    !is.na(neon_plfa$cyclo19To0ScaledConcentration) &
    !is.na(stress_den) & stress_den > 0,
  stress_num / stress_den,
  NA)

cat("Stress index (cy19:0/18:1w7c) calculated for",
    sum(!is.na(neon_plfa$stress_index_cy19)), "samples;",
    "cy17:0-based for", sum(!is.na(neon_plfa$stress_index_cy17)), "samples;",
    "combined for", sum(!is.na(neon_plfa$stress_index_combined)), "samples\n")

cat("\n=== Generating clean dataset ===\n")

# Redundant join-artifact columns: byte-identical to their unsuffixed counterpart
# because core/moisture/chemistry all carry their own copies of sampleID's parent
# fields (sampleID, horizon, plotID, collectDate, year, month).
duplicate_join_cols <- c(
  "horizon.core", "plotID.core", "collectDateTime.core", "collectDate.core", "year.core", "month.core",
  "horizon.moist", "plotID.moist", "collectDateTime.moist", "collectDate.moist", "year.moist", "month.moist",
  "plotID.chem", "collectDateTime.chem", "collectDate.chem", "year.chem", "month.chem"
)

# total_microbial_biomass is a literal copy of microbial_biomass_nmol_g (see
# compute_soil_health_metrics() in R/neon_plfa.R); only the latter is published.
# The superseded legacy single-column "stress_index" is dropped in favor of the
# three newly derived stress indices (stress_index_cy19 / cy17 / combined).
drop_cols <- c(duplicate_join_cols, "total_microbial_biomass", "stress_index")

neon_plfa_clean <- neon_plfa %>%
  select(-any_of(drop_cols))

write.csv(neon_plfa_clean,
          "manuscript_outputs/data/neon_plfa_synthesis_v1.2.csv",
          row.names = FALSE)

cat("✓ Saved: neon_plfa_synthesis_v1.2.csv\n")

cat("\n=== Generating supplementary tables ===\n")

site_metadata <- neon_plfa_clean %>%
  group_by(siteID) %>%
  summarize(
    domain = first(domainID),
    latitude = round(first(decimalLatitude), 4),
    longitude = round(first(decimalLongitude), 4),
    elevation_m = round(first(elevation), 1),
    soil_order = first(soilOrder),
    soil_suborder = first(soilSuborder),
    n_samples = n(),
    n_years = n_distinct(year),
    years_sampled = paste(sort(unique(year)), collapse=", "),
    first_sample = min(collectDate, na.rm = TRUE),
    last_sample = max(collectDate, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(domain, siteID)

write.csv(site_metadata, "manuscript_outputs/tables/Table1_site_metadata.csv", row.names = FALSE)

calc_summary <- function(x) {
  x <- x[!is.na(x)]
  data.frame(
    n = length(x),
    mean = round(mean(x), 2),
    sd = round(sd(x), 2),
    min = round(min(x), 2),
    q25 = round(quantile(x, 0.25), 2),
    median = round(median(x), 2),
    q75 = round(quantile(x, 0.75), 2),
    max = round(max(x), 2)
  )
}

summary_stats <- bind_rows(
  calc_summary(neon_plfa_clean$microbial_biomass_nmol_g) %>%
    mutate(metric = "Microbial biomass (nmol/g)", .before = 1),
  calc_summary(neon_plfa_clean$fungal_plfa_nmol_g) %>%
    mutate(metric = "Fungal PLFA (nmol/g)", .before = 1),
  calc_summary(neon_plfa_clean$bacterial_plfa_nmol_g) %>%
    mutate(metric = "Bacterial PLFA (nmol/g)", .before = 1),
  calc_summary(neon_plfa_clean$FB_ratio) %>%
    mutate(metric = "F:B ratio", .before = 1),
  calc_summary(neon_plfa_clean$stress_index_cy19) %>%
    mutate(metric = "Cyclopropyl:precursor stress index [cy19:0/18:1ω7c]", .before = 1),
  calc_summary(neon_plfa_clean$stress_index_cy17) %>%
    mutate(metric = "Cyclopropyl:precursor stress index [cy17:0/16:1ω7c; 2021+]", .before = 1),
  calc_summary(neon_plfa_clean$stress_index_combined) %>%
    mutate(metric = "Cyclopropyl:precursor stress index [combined; 2021+]", .before = 1),
  calc_summary(neon_plfa_clean$GP_GN_ratio) %>%
    mutate(metric = "GP:GN ratio", .before = 1)
)

write.csv(summary_stats, "manuscript_outputs/tables/Table2_summary_statistics.csv", row.names = FALSE)

temporal_matrix <- neon_plfa_clean %>%
  group_by(siteID, year) %>%
  summarize(n_samples = n(), .groups = "drop") %>%
  pivot_wider(names_from = year, values_from = n_samples, values_fill = 0)
write.csv(temporal_matrix, "manuscript_outputs/tables/temporal_coverage_matrix.csv", row.names = FALSE)

completeness_by_site <- neon_plfa_clean %>%
  group_by(siteID, domainID) %>%
  summarize(
    n_total = n(),
    pct_biomass = round(100 * sum(!is.na(microbial_biomass_nmol_g)) / n(), 1),
    pct_fungal = round(100 * sum(!is.na(fungal_plfa_nmol_g)) / n(), 1),
    pct_bacterial = round(100 * sum(!is.na(bacterial_plfa_nmol_g)) / n(), 1),
    pct_FB = round(100 * sum(!is.na(FB_ratio)) / n(), 1),
    pct_stress = round(100 * sum(!is.na(stress_index_cy19)) / n(), 1),
    pct_GPGN = round(100 * sum(!is.na(GP_GN_ratio)) / n(), 1),
    .groups = "drop"
  )
write.csv(completeness_by_site, "manuscript_outputs/tables/completeness_by_site.csv", row.names = FALSE)

horizon_summary <- neon_plfa_clean %>%
  group_by(horizon) %>%
  summarize(
    n_samples = n(),
    pct_total = round(100 * n() / nrow(neon_plfa_clean), 1),
    n_sites = n_distinct(siteID),
    mean_biomass = round(mean(microbial_biomass_nmol_g, na.rm = TRUE), 2),
    sd_biomass = round(sd(microbial_biomass_nmol_g, na.rm = TRUE), 2),
    mean_FB = round(mean(FB_ratio, na.rm = TRUE), 3),
    sd_FB = round(sd(FB_ratio, na.rm = TRUE), 3),
    mean_stress = round(mean(stress_index_cy19, na.rm = TRUE), 2),
    sd_stress = round(sd(stress_index_cy19, na.rm = TRUE), 2),
    mean_GPGN = round(mean(GP_GN_ratio, na.rm = TRUE), 2),
    sd_GPGN = round(sd(GP_GN_ratio, na.rm = TRUE), 2),
    .groups = "drop"
  )
write.csv(horizon_summary, "manuscript_outputs/tables/horizon_summary.csv", row.names = FALSE)

domain_summary <- neon_plfa_clean %>%
  group_by(domainID) %>%
  summarize(
    n_samples = n(),
    n_sites = n_distinct(siteID),
    mean_biomass = round(mean(microbial_biomass_nmol_g, na.rm = TRUE), 2),
    sd_biomass = round(sd(microbial_biomass_nmol_g, na.rm = TRUE), 2),
    mean_FB = round(mean(FB_ratio, na.rm = TRUE), 3),
    sd_FB = round(sd(FB_ratio, na.rm = TRUE), 3),
    mean_stress = round(mean(stress_index_cy19, na.rm = TRUE), 2),
    sd_stress = round(sd(stress_index_cy19, na.rm = TRUE), 2),
    mean_GPGN = round(mean(GP_GN_ratio, na.rm = TRUE), 2),
    sd_GPGN = round(sd(GP_GN_ratio, na.rm = TRUE), 2),
    .groups = "drop"
  )
write.csv(domain_summary, "manuscript_outputs/tables/domain_summary.csv", row.names = FALSE)

cat("✓ All tables generated\n")


cat("\n=== Generating Figures ===\n")
theme_set(theme_bw(base_size = 11))

cat("Generating Figure 1a: Site Map...\n")

north_america <- map_data("world", region = c("USA", "Canada", "Mexico", "Puerto Rico"))

fig1a <- ggplot() +
  geom_polygon(data = north_america,
               aes(x = long, y = lat, group = group),
               fill = "gray95", color = "gray70", linewidth = 0.3) +
  geom_point(data = site_metadata,
             aes(x = longitude, y = latitude, size = n_samples, fill = domain),
             shape = 21, alpha = 0.8, stroke = 0.5) +
  scale_size_continuous(range = c(2, 10), name = "N samples",
                        breaks = c(100, 200, 300, 400, 500)) +
  scale_fill_viridis_d(option = "turbo", name = "Domain") +
  coord_fixed(1.3, xlim = c(-180, -50), ylim = c(14, 72)) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(legend.position = "right",
        legend.box = "vertical",
        panel.background = element_rect(fill = "aliceblue"),
        plot.margin = margin(10, 40, 10, 10))  # Extra right margin for legend

ggsave("manuscript_outputs/figures/Figure1a_site_map.png", fig1a,
       width = 12, height = 7, dpi = 300)
ggsave("manuscript_outputs/figures/Figure1a_site_map.pdf", fig1a,
       width = 12, height = 7)

cat("✓ Figure 1a saved\n")


cat("Generating Figure 1b: Temporal Coverage...\n")

temporal_long <- temporal_matrix %>%
  pivot_longer(-siteID, names_to = "year", values_to = "n_samples") %>%
  mutate(year = as.numeric(year))

fig1b <- ggplot(temporal_long, aes(x = year, y = siteID, fill = n_samples)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_viridis_c(option = "plasma", name = "N samples",
                       breaks = c(0, 50, 100, 200, 300, 400, 500)) +
  scale_x_continuous(breaks = 2017:2024) +
  labs(x = "Year", y = "Site") +
  theme_minimal(base_size = 9) +
  theme(axis.text.y = element_text(size = 6),
        panel.grid = element_blank())

ggsave("manuscript_outputs/figures/Figure1b_temporal_heatmap.png", fig1b,
       width = 8, height = 10, dpi = 300)
ggsave("manuscript_outputs/figures/Figure1b_temporal_heatmap.pdf", fig1b,
       width = 8, height = 10)

cat("✓ Figure 1b saved\n")

cat("Generating Figure 1c: Horizon Distribution...\n")

fig1c <- ggplot(horizon_summary, aes(x = "", y = n_samples, fill = horizon)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y") +
  scale_fill_manual(values = c("M" = "#8B4513", "O" = "#2E8B57"),
                    labels = c("M" = paste0("Mineral (", horizon_summary$pct_total[horizon_summary$horizon == "M"], "%)"),
                               "O" = paste0("Organic (", horizon_summary$pct_total[horizon_summary$horizon == "O"], "%)"))) +
  labs(fill = "Horizon") +
  theme_void()

ggsave("manuscript_outputs/figures/Figure1c_horizon_distribution.png", fig1c,
       width = 6, height = 5, dpi = 300)

cat("✓ Figure 1c saved\n")

cat("Generating Figure 2: Metric Distributions...\n")

# Shared styling for the cyclopropyl:precursor stress index (primary, cy19:0/18:1ω7c)
ST_COLOR <- "darkorange3"
ST_FILL  <- "navajowhite"
STRESS_LAB_FULL  <- "Cyclopropyl:precursor stress ratio (cy19:0/18:1ω7c)"
STRESS_LAB_MEAN  <- "Mean stress ratio (cy19:0/18:1ω7c)"
STRESS_LAB_SHORT <- "Stress ratio (cy19:0/18:1ω7c)"
bold_title <- theme(plot.title = element_text(face = "bold", size = 14, hjust = 0))

fig2a <- ggplot(neon_plfa_clean, aes(x = log10(microbial_biomass_nmol_g))) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
  labs(x = "log10(Biomass) [nmol/g]", y = "Count") +
  ggtitle("a") + theme_minimal() + bold_title

fig2b <- ggplot(neon_plfa_clean, aes(x = log10(FB_ratio + 0.001))) +
  geom_histogram(bins = 50, fill = "darkgreen", alpha = 0.7) +
  labs(x = "log10(Fungal:Bacterial Ratio + 0.001)", y = "Count") +
  ggtitle("b") + theme_minimal() + bold_title

# Stress index truncated at 5 for display (long right tail; see Table 2 for full range)
fig2c <- ggplot(neon_plfa_clean %>% filter(!is.na(stress_index_cy19) & stress_index_cy19 < 5),
                aes(x = stress_index_cy19)) +
  geom_histogram(bins = 50, fill = ST_COLOR, alpha = 0.7) +
  labs(x = STRESS_LAB_FULL, y = "Count") +
  ggtitle("c") + theme_minimal() + bold_title

fig2d <- ggplot(neon_plfa_clean %>% filter(!is.na(GP_GN_ratio) & GP_GN_ratio < 10),
                aes(x = log10(GP_GN_ratio))) +
  geom_histogram(bins = 50, fill = "purple", alpha = 0.7) +
  labs(x = "log10(Gram+:Gram- Ratio)", y = "Count") +
  ggtitle("d") + theme_minimal() + bold_title

save_grob("manuscript_outputs/figures/Figure2_distributions.png",
          arrangeGrob(fig2a, fig2b, fig2c, fig2d, ncol = 1), 8, 13)

cat("✓ Figure 2 saved\n")

cat("Generating Figure 3: Latitudinal Gradients...\n")

site_averages <- neon_plfa_clean %>%
  group_by(siteID) %>%
  summarize(
    latitude = first(decimalLatitude),
    mean_biomass = mean(microbial_biomass_nmol_g, na.rm = TRUE),
    mean_FB = mean(FB_ratio, na.rm = TRUE),
    mean_stress = mean(stress_index_cy19, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  ) %>%
  mutate(log_mean_biomass = log10(mean_biomass))

add_lm_stats <- function(data, x_col, y_col) {
  df <- data.frame(x = data[[x_col]], y = data[[y_col]]) %>%
    filter(!is.na(x) & !is.na(y))

  if (nrow(df) < 3) return(NULL)

  fit <- lm(y ~ x, data = df)
  fit_summary <- summary(fit)

  r2 <- round(fit_summary$r.squared, 3)
  pval <- fit_summary$coefficients["x", "Pr(>|t|)"]
  pval_text <- if (pval < 0.001) "p < 0.001" else paste0("p = ", round(pval, 3))
  n <- nrow(df)

  list(
    label = paste0("N = ", n, "\nR² = ", r2, "\n", pval_text),
    r2 = r2,
    pval = pval,
    n = n
  )
}

# Fit on log10(mean_biomass) so the reported R²/p-value match the log10-transformed
# axis plotted below and the "log10(Biomass)" label in Table S1 (previously the lm()
# was fit on raw mean_biomass while the figure/label implied a log10 fit).
stats_biomass <- add_lm_stats(site_averages, "latitude", "log_mean_biomass")

fig3a <- ggplot(site_averages, aes(x = latitude, y = log_mean_biomass)) +
  geom_point(size = 3, alpha = 0.7, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "darkblue", fill = "lightblue") +
  annotate("text", x = 20, y = 2.9,
           label = stats_biomass$label, hjust = 0, vjust = 1, size = 3.5,
           fontface = "bold") +
  labs(x = "Latitude (°N)", y = "log10(Mean Biomass) [nmol/g]") +
  ggtitle("a") +
  theme_minimal() +
  theme(plot.margin = margin(10, 10, 10, 10),
        plot.title = element_text(face = "bold", size = 14, hjust = 0))

stats_FB <- add_lm_stats(site_averages, "latitude", "mean_FB")

fig3b <- ggplot(site_averages, aes(x = latitude, y = mean_FB)) +
  geom_point(size = 3, alpha = 0.7, color = "darkgreen") +
  geom_smooth(method = "lm", se = TRUE, color = "darkgreen", fill = "lightgreen") +
  annotate("text", x = 20, y = 0.29,
           label = stats_FB$label, hjust = 0, vjust = 1, size = 3.5,
           fontface = "bold") +
  labs(x = "Latitude (°N)", y = "Mean F:B Ratio") +
  ggtitle("b") +
  theme_minimal() +
  theme(plot.margin = margin(10, 10, 10, 10),
        plot.title = element_text(face = "bold", size = 14, hjust = 0))

stats_stress <- add_lm_stats(site_averages, "latitude", "mean_stress")

fig3c <- ggplot(site_averages, aes(x = latitude, y = mean_stress)) +
  geom_point(size = 3, alpha = 0.7, color = ST_COLOR) +
  geom_smooth(method = "lm", se = TRUE, color = ST_COLOR, fill = ST_FILL) +
  annotate("text", x = 45, y = max(site_averages$mean_stress, na.rm = TRUE) * 0.98,
           label = stats_stress$label, hjust = 0, vjust = 1, size = 3.5,
           fontface = "bold") +
  labs(x = "Latitude (°N)", y = STRESS_LAB_MEAN) +
  ggtitle("c") +
  theme_minimal() +
  theme(plot.margin = margin(10, 10, 10, 10),
        plot.title = element_text(face = "bold", size = 14, hjust = 0))

save_grob("manuscript_outputs/figures/Figure3_latitudinal_gradients.png",
          arrangeGrob(fig3a, fig3b, fig3c, ncol = 1), 8, 14)

lat_stats <- data.frame(
  metric = c("log10(Biomass)", "F:B Ratio", "Cyclopropyl:precursor stress index"),
  n_sites = c(stats_biomass$n, stats_FB$n, stats_stress$n),
  r_squared = c(stats_biomass$r2, stats_FB$r2, stats_stress$r2),
  p_value = c(stats_biomass$pval, stats_FB$pval, stats_stress$pval)
)
write.csv(lat_stats, "manuscript_outputs/tables/latitudinal_gradient_stats.csv", row.names = FALSE)

cat("✓ Figure 3 saved\n")


cat("Generating Figure 4: Domain patterns...\n")

fig4a <- ggplot(neon_plfa_clean, aes(x = domainID, y = log10(microbial_biomass_nmol_g))) +
  geom_boxplot(fill = "steelblue", alpha = 0.6, outlier.size = 0.5) +
  labs(x = "NEON Domain", y = "log10(Biomass) [nmol/g]") +
  ggtitle("a") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 14, hjust = 0))

fig4b <- ggplot(neon_plfa_clean, aes(x = domainID, y = FB_ratio)) +
  geom_boxplot(fill = "darkgreen", alpha = 0.6, outlier.size = 0.5) +
  labs(x = "NEON Domain", y = "F:B Ratio") +
  ggtitle("b") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 14, hjust = 0))

fig4c <- ggplot(neon_plfa_clean %>% filter(!is.na(stress_index_cy19)),
                aes(x = domainID, y = stress_index_cy19)) +
  geom_boxplot(fill = ST_COLOR, alpha = 0.6, outlier.size = 0.5) +
  coord_cartesian(ylim = c(0, quantile(neon_plfa_clean$stress_index_cy19, 0.99, na.rm = TRUE))) +
  labs(x = "NEON Domain", y = STRESS_LAB_SHORT) +
  ggtitle("c") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 14, hjust = 0))

save_grob("manuscript_outputs/figures/Figure4_domain_patterns.png",
          arrangeGrob(fig4a, fig4b, fig4c, ncol = 1), 10, 15)

cat("✓ Figure 4 saved\n")


cat("Generating Figure 5: Data Quality...\n")

plfa_completeness <- data.frame(
  metric = c("Biomass", "Fungal PLFA", "Bacterial PLFA", "F:B Ratio", "Stress Index", "GP:GN Ratio"),
  pct_complete = c(
    100 * sum(!is.na(neon_plfa_clean$microbial_biomass_nmol_g)) / nrow(neon_plfa_clean),
    100 * sum(!is.na(neon_plfa_clean$fungal_plfa_nmol_g)) / nrow(neon_plfa_clean),
    100 * sum(!is.na(neon_plfa_clean$bacterial_plfa_nmol_g)) / nrow(neon_plfa_clean),
    100 * sum(!is.na(neon_plfa_clean$FB_ratio)) / nrow(neon_plfa_clean),
    100 * sum(!is.na(neon_plfa_clean$stress_index_cy19)) / nrow(neon_plfa_clean),
    100 * sum(!is.na(neon_plfa_clean$GP_GN_ratio)) / nrow(neon_plfa_clean)
  )
)

fig5a <- ggplot(plfa_completeness, aes(x = reorder(metric, pct_complete), y = pct_complete)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
  geom_text(aes(label = paste0(round(pct_complete, 1), "%")),
            hjust = -0.1, size = 3.8) +
  coord_flip() +
  ylim(0, 105) +
  labs(x = "", y = "% Complete") +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 13, face = "bold"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 11),
        panel.grid.major.y = element_blank(),
        plot.margin = margin(10, 20, 10, 10))

site_year_means <- neon_plfa_clean %>%
  group_by(siteID, year) %>%
  summarize(
    mean_biomass = mean(microbial_biomass_nmol_g, na.rm = TRUE),
    mean_FB = mean(FB_ratio, na.rm = TRUE),
    mean_stress = mean(stress_index_cy19, na.rm = TRUE),
    mean_GPGN = mean(GP_GN_ratio, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  )

temporal_stability <- site_year_means %>%
  group_by(siteID) %>%
  summarize(
    n_years = n(),
    cv_biomass = 100 * sd(mean_biomass, na.rm = TRUE) / mean(mean_biomass, na.rm = TRUE),
    cv_FB = 100 * sd(mean_FB, na.rm = TRUE) / mean(mean_FB, na.rm = TRUE),
    cv_stress = 100 * sd(mean_stress, na.rm = TRUE) / mean(mean_stress, na.rm = TRUE),
    cv_GPGN = 100 * sd(mean_GPGN, na.rm = TRUE) / mean(mean_GPGN, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_years >= 3)

write.csv(temporal_stability,
          "manuscript_outputs/tables/temporal_stability_cv.csv",
          row.names = FALSE)

cv_summary <- data.frame(
  metric = c("Biomass", "F:B Ratio", "GP:GN Ratio", "Cyclopropyl:precursor stress index"),
  n_sites = c(
    sum(!is.na(temporal_stability$cv_biomass)),
    sum(!is.na(temporal_stability$cv_FB)),
    sum(!is.na(temporal_stability$cv_GPGN)),
    sum(!is.na(temporal_stability$cv_stress))
  ),
  median_cv = c(
    round(median(temporal_stability$cv_biomass, na.rm = TRUE), 1),
    round(median(temporal_stability$cv_FB, na.rm = TRUE), 1),
    round(median(temporal_stability$cv_GPGN, na.rm = TRUE), 1),
    round(median(temporal_stability$cv_stress, na.rm = TRUE), 1)
  ),
  mean_cv = c(
    round(mean(temporal_stability$cv_biomass, na.rm = TRUE), 1),
    round(mean(temporal_stability$cv_FB, na.rm = TRUE), 1),
    round(mean(temporal_stability$cv_GPGN, na.rm = TRUE), 1),
    round(mean(temporal_stability$cv_stress, na.rm = TRUE), 1)
  )
)

write.csv(cv_summary,
          "manuscript_outputs/tables/cv_summary_stats.csv",
          row.names = FALSE)

cat("Temporal stability: n =", nrow(temporal_stability), "sites with ≥3 years\n")
cat("Median CVs: Biomass =", cv_summary$median_cv[1], "%, F:B =", cv_summary$median_cv[2],
    "%, GP:GN =", cv_summary$median_cv[3], "%, Stress =", cv_summary$median_cv[4], "%\n")

calc_icc <- function(data, value_col) {
  model_data <- data %>%
    select(siteID, year, value = all_of(value_col)) %>%
    filter(!is.na(value))

  if(nrow(model_data) < 10) return(NA)

  site_means <- model_data %>%
    group_by(siteID) %>%
    summarize(site_mean = mean(value, na.rm = TRUE), .groups = "drop")

  grand_mean <- mean(model_data$value, na.rm = TRUE)
  between_var <- var(site_means$site_mean)

  within_var <- model_data %>%
    left_join(site_means, by = "siteID") %>%
    mutate(dev = (value - site_mean)^2) %>%
    pull(dev) %>%
    mean(na.rm = TRUE)

  icc <- between_var / (between_var + within_var)
  return(icc)
}

icc_biomass <- calc_icc(site_year_means, "mean_biomass")
icc_FB <- calc_icc(site_year_means, "mean_FB")
icc_GPGN <- calc_icc(site_year_means, "mean_GPGN")
icc_stress <- calc_icc(site_year_means, "mean_stress")

icc_results <- data.frame(
  metric = c("Biomass", "F:B Ratio", "GP:GN Ratio", "Cyclopropyl:precursor stress index"),
  ICC = round(c(icc_biomass, icc_FB, icc_GPGN, icc_stress), 3)
)

write.csv(icc_results,
          "manuscript_outputs/tables/icc_results.csv",
          row.names = FALSE)

cat("ICCs: Biomass =", round(icc_biomass, 3), ", F:B =", round(icc_FB, 3),
    ", GP:GN =", round(icc_GPGN, 3), ", Stress =", round(icc_stress, 3), "\n")

cv_long <- temporal_stability %>%
  select(siteID, n_years, cv_biomass, cv_FB, cv_stress, cv_GPGN) %>%
  pivot_longer(cols = starts_with("cv_"),
               names_to = "metric",
               values_to = "cv") %>%
  mutate(metric = recode(metric,
                         "cv_biomass" = "Biomass",
                         "cv_FB" = "F:B Ratio",
                         "cv_stress" = "Stress Index",
                         "cv_GPGN" = "GP:GN Ratio"),
         metric = factor(metric, levels = c("Biomass", "F:B Ratio", "Stress Index", "GP:GN Ratio")))

fig5b <- ggplot(cv_long, aes(x = metric, y = cv)) +
  geom_boxplot(fill = "lightblue", alpha = 0.7, outlier.alpha = 0.3) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1.5) +
  geom_hline(yintercept = c(20, 30, 40), linetype = "dashed",
             color = "gray60", linewidth = 0.3) +
  labs(x = "", y = "Coefficient of Variation (%)") +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 13, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        axis.text.y = element_text(size = 11),
        plot.margin = margin(10, 10, 10, 20))

# Combine completeness (a) and temporal stability (b) into a single two-panel figure
fig5a_lab <- arrangeGrob(fig5a, top = textGrob("a", x = 0.02, hjust = 0,
                                               gp = gpar(fontface = "bold", fontsize = 16)))
fig5b_lab <- arrangeGrob(fig5b, top = textGrob("b", x = 0.02, hjust = 0,
                                               gp = gpar(fontface = "bold", fontsize = 16)))
save_grob("manuscript_outputs/figures/Figure5_combined.png",
          arrangeGrob(fig5a_lab, fig5b_lab, ncol = 2, widths = c(1, 1)), 14, 6)

cat("✓ Figure 5 saved\n")


cat("Generating Figure 6: PLFA Metric Relationships with statistics...\n")

get_stats_label <- function(data, x_var, y_var) {
  df <- data.frame(x = data[[x_var]], y = data[[y_var]]) %>%
    filter(!is.na(x) & !is.na(y) & is.finite(x) & is.finite(y))

  if (nrow(df) < 3) return(list(label = "Insufficient data", n = 0))

  fit <- lm(y ~ x, data = df)
  fit_summary <- summary(fit)

  r2 <- round(fit_summary$r.squared, 3)
  pval <- fit_summary$coefficients["x", "Pr(>|t|)"]
  pval_text <- if (pval < 0.001) "p < 0.001" else paste0("p = ", round(pval, 3))

  list(
    label = paste0("N = ", nrow(df), "\nR² = ", r2, "\n", pval_text),
    n = nrow(df),
    r2 = r2,
    pval = pval
  )
}

data_6a <- neon_plfa_clean %>%
  mutate(log_biomass = log10(microbial_biomass_nmol_g)) %>%
  filter(!is.na(log_biomass) & !is.na(FB_ratio) & is.finite(log_biomass))

stats_6a <- get_stats_label(data_6a, "log_biomass", "FB_ratio")

fig6a <- ggplot(data_6a, aes(x = log_biomass, y = FB_ratio)) +
  geom_point(alpha = 0.2, size = 0.8, color = "steelblue") +
  geom_smooth(method = "lm", color = "darkblue", se = TRUE, fill = "lightblue") +
  annotate("text", x = 0.8, y = 2.1,
           label = stats_6a$label, hjust = 0, vjust = 1, size = 3.5,
           fontface = "bold") +
  labs(x = "log10(Biomass) [nmol/g]", y = "F:B Ratio") +
  ggtitle("a") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0))

data_6b <- neon_plfa_clean %>%
  mutate(log_bacterial = log10(bacterial_plfa_nmol_g),
         log_fungal = log10(fungal_plfa_nmol_g + 0.1)) %>%
  filter(!is.na(log_bacterial) & !is.na(log_fungal) &
         is.finite(log_bacterial) & is.finite(log_fungal))

stats_6b <- get_stats_label(data_6b, "log_bacterial", "log_fungal")

fig6b <- ggplot(data_6b, aes(x = log_bacterial, y = log_fungal)) +
  geom_point(alpha = 0.2, size = 0.8, color = "darkgreen") +
  geom_smooth(method = "lm", color = "darkgreen", se = TRUE, fill = "lightgreen") +
  annotate("text", x = 0.2, y = 2.7,
           label = stats_6b$label, hjust = 0, vjust = 1, size = 3.5,
           fontface = "bold") +
  labs(x = "log10(Bacterial PLFA) [nmol/g]",
       y = "log10(Fungal PLFA + 0.1) [nmol/g]") +
  ggtitle("b") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0))

# Panel (c): stress index vs F:B ratio. Statistics use all paired observations;
# display truncates the stress axis at the 99th percentile for readability.
data_6c_full <- neon_plfa_clean %>%
  filter(is.finite(stress_index_cy19) & !is.na(FB_ratio))
stats_6c <- get_stats_label(data_6c_full, "FB_ratio", "stress_index_cy19")
data_6c <- data_6c_full %>%
  filter(stress_index_cy19 < quantile(stress_index_cy19, 0.99, na.rm = TRUE))

fig6c <- ggplot(data_6c, aes(x = FB_ratio, y = stress_index_cy19)) +
  geom_point(alpha = 0.2, size = 0.8, color = ST_COLOR) +
  geom_smooth(method = "lm", color = ST_COLOR, se = TRUE, fill = ST_FILL) +
  annotate("text", x = max(data_6c$FB_ratio, na.rm = TRUE) * 0.6,
           y = max(data_6c$stress_index_cy19, na.rm = TRUE) * 0.95,
           label = stats_6c$label, hjust = 0, vjust = 1, size = 3.5,
           fontface = "bold") +
  labs(x = "F:B Ratio", y = STRESS_LAB_SHORT) +
  ggtitle("c") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0))

save_grob("manuscript_outputs/figures/Figure6_plfa_relationships.png",
          arrangeGrob(fig6a, fig6b, fig6c, ncol = 1), 8, 14)

fig6_stats <- data.frame(
  panel = c("(a) F:B vs Biomass", "(b) Fungal vs Bacterial", "(c) Stress index vs F:B"),
  n_samples = c(stats_6a$n, stats_6b$n, stats_6c$n),
  r_squared = c(stats_6a$r2, stats_6b$r2, stats_6c$r2),
  p_value = c(stats_6a$pval, stats_6b$pval, stats_6c$pval)
)
write.csv(fig6_stats, "manuscript_outputs/tables/figure6_relationship_stats.csv", row.names = FALSE)

cat("✓ Figure 6 saved\n")


cat("\n=== Saving R workspace ===\n")
save.image("manuscript_outputs/neon_plfa_analysis_workspace.RData")
cat("✓ Saved workspace\n")