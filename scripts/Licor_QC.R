# =============================================================================
# LI-6800 QC: leaf-area correction + flagging of problematic values
#
# Input:   licor2       one row per logged observation
# Output:  licor_qc     all rows, area-corrected values (*_c) + flag columns
#          licor_clean  "exclude" rows removed, iWUE masked where unreliable
#          review_table rows needing a look, with reasons
#
# Corrected columns end in _c (A_c, E_c, Emm_c, gsw_c, gtc_c, Ci_c, Pci_c,
# CiCa_c, WUE_c). Original columns are left untouched.
# =============================================================================

# This script corrects A, E, gsw, gtc, Ci and iWUE for the actual leaf area in licor_qc.R. It then sorts each observation into “ok,” “review,” or “exclude” and lists the reason for each flag.

# How the correction works: A, E, and total conductance scale exactly by 6 / actual area, so those are rescaled directly. gsw, gtc, and Ci are then recomputed with the LI-6800 equations. The boundary-layer terms are back-calculated from each row’s own output, so the correction keeps whatever settings the machine used. Section 3 checks that the equations reproduce the machine’s logged Ci and WUE. Both differences should be near zero; if they aren’t, stop and let me know.

# Things to know before running it:
  
#   Check what area was entered on each machine. The script assumes both machines used 6 cm². If the second machine had a different area entered, the correction factor for those rows will be wrong.
# Leaves larger than 6 cm² are capped at 6, not scaled. The chamber can only hold 6 cm² of leaf, so the script treats an area over 6 as a full chamber and flags it for review.
# Leaves under 1 cm² are flagged. The correction multiplies their values more than 6×, which also magnifies measurement noise.
# iWUE will barely change after correction. A and gsw both scale with area, so the area mostly cancels out of the ratio. A and gsw themselves will change a lot.
# Thresholds are set for C3 plants. They’re grouped at the top of the script and are reasonable for Bromus. Adjust them if you have C4 species.

# The flag_summary table shows whether flags cluster on the second machine. Rows flagged “review” stay in licor_clean until you’ve looked at them in review_table.

library(dplyr)
library(tidyr)
library(ggplot2)

# ---- 1. Settings (edit these) -----------------------------------------------

# Leaf area (cm^2) the LI-6800 used in its calculations. If the second machine
# had a different area entered, make this a per-row column instead.
area_entered <- 6

thr <- list(
  area_min    = 1,     # cm^2; below this the correction (>6x) amplifies noise
  A_max       = 40,    # µmol m-2 s-1; implausibly high for most C3 leaves
  gsw_max     = 1.5,   # mol m-2 s-1; implausibly high
  gsw_min_wue = 0.02,  # mol m-2 s-1; below this A/gsw is unreliable
  wue_max     = 250,   # µmol mol-1; implausible for C3
  CiCa_min    = 0.2,   # very low Ci/Ca with positive A suggests leak/error (C3)
  CiCa_max    = 1.0,   # Ci > Ca with positive A is physically impossible
  rep_dev     = 0.30,  # relative deviation of an obs from its plant's median A
  mad_z       = 3.5    # robust z-score cutoff for species-level outliers
)

# ---- 2. Leaf-area correction ------------------------------------------------
# A, E and gtw are all divided by leaf area in the LI-6800 equations, so they
# scale exactly by (area_entered / actual_area). gsw, gtc and Ci are then
# recomputed from the corrected values. The boundary-layer terms are
# back-calculated from each row's own output, so the correction keeps whatever
# boundary-layer conductance and stomatal ratio the machine used.

licor2 <- read.csv(here::here("data/licor2.csv"))

licor_qc <- licor2 %>%
  mutate(
    WUE_source = if_else(is.na(WUE), "manual (2nd machine)", "machine"),
    
    leaf_area = case_when(
      is.na(total.leaf.area) | total.leaf.area <= 0 ~ NA_real_,
      total.leaf.area > area_entered ~ area_entered,  # leaf larger than aperture
      TRUE ~ total.leaf.area
    ),
    area_factor = area_entered / leaf_area,
    
    # Back-calculated terms from the machine's own output
    r_bw  = 1 / gtw - 1 / gsw,        # boundary-layer term, water
    r_bc  = 1 / gtc - 1.6 / gsw,      # boundary-layer term, CO2
    P_kPa = 1000 * Pca / Ca,          # chamber pressure
    
    # Fluxes scale directly with 1/area
    A_c   = A   * area_factor,
    E_c   = E   * area_factor,
    Emm_c = E_c * 1000,
    gtw_c = gtw * area_factor,
    
    # Recompute conductances and Ci
    gsw_den = 1 / gtw_c - r_bw,
    gsw_c   = if_else(gsw_den > 0, 1 / gsw_den, NA_real_),
    gtc_c   = 1 / (1.6 / gsw_c + r_bc),
    Ci_c    = ((gtc_c - E_c / 2) * Ca - A_c) / (gtc_c + E_c / 2),
    Pci_c   = Ci_c * P_kPa / 1000,
    CiCa_c  = Ci_c / Ca,
    
    # Intrinsic WUE (A/gsw) for every row, including the missing ones
    WUE_c   = A_c / gsw_c
  )

# ---- 3. Sanity checks: do the equations reproduce the machine's output? -----
# Both differences should be tiny (rounding only). If not, stop and check.

ci_diff  <- with(licor_qc, ((gtc - E / 2) * Ca - A) / (gtc + E / 2) - Ci)
wue_diff <- with(licor_qc, A / gsw - WUE)
cat("Max |Ci reconstructed - Ci logged|: ",
    round(max(abs(ci_diff), na.rm = TRUE), 3), "ppm\n")
cat("Max |A/gsw - WUE logged|:           ",
    round(max(abs(wue_diff), na.rm = TRUE), 3), "\n")

# ---- 4. Replicate stability and species-level outliers ----------------------

safe_z <- function(x) {
  m <- mad(x, na.rm = TRUE)
  if (is.na(m) || m == 0) return(rep(NA_real_, length(x)))
  (x - median(x, na.rm = TRUE)) / m
}

licor_qc <- licor_qc %>%
  group_by(Plant_ID, Date) %>%
  mutate(
    n_obs   = n(),
    A_med   = median(A_c, na.rm = TRUE),
    rep_dev = abs(A_c - A_med) / pmax(abs(A_med), 1)  # floor avoids blow-up near 0
  ) %>%
  group_by(Sp_Code) %>%
  mutate(
    z_A   = safe_z(A_c),
    z_WUE = safe_z(WUE_c)
  ) %>%
  ungroup()

# ---- 5. Flags ---------------------------------------------------------------
# h_ = hard flags: technical failure or physically impossible -> exclude
# s_ = soft flags: unusual but possibly real biology          -> review

licor_qc <- licor_qc %>%
  mutate(
    # Hard
    h_area_missing      = is.na(leaf_area),
    h_correction_failed = !is.na(leaf_area) & is.na(gsw_c),
    h_Ci_negative       = !is.na(Ci_c) & Ci_c < 0,
    h_Ci_gt_Ca          = !is.na(CiCa_c) & CiCa_c > thr$CiCa_max & A_c > 0,
    
    # Soft
    s_area_capped  = !is.na(total.leaf.area) & total.leaf.area > area_entered,
    s_area_small   = !is.na(leaf_area) & leaf_area < thr$area_min,
    s_A_negative   = !is.na(A_c) & A_c < 0,
    s_A_high       = !is.na(A_c) & A_c > thr$A_max,
    s_E_negative   = !is.na(E_c) & E_c < 0,
    s_gsw_low      = !is.na(gsw_c) & gsw_c < thr$gsw_min_wue,
    s_gsw_high     = !is.na(gsw_c) & gsw_c > thr$gsw_max,
    s_WUE_high     = !is.na(WUE_c) & WUE_c > thr$wue_max,
    s_CiCa_low     = !is.na(CiCa_c) & CiCa_c < thr$CiCa_min & A_c > 0,
    s_unstable_rep = n_obs > 1 & !is.na(rep_dev) & rep_dev > thr$rep_dev,
    s_outlier_A    = !is.na(z_A)   & abs(z_A)   > thr$mad_z,
    s_outlier_WUE  = !is.na(z_WUE) & abs(z_WUE) > thr$mad_z
  )

hard_flags <- grep("^h_", names(licor_qc), value = TRUE)
soft_flags <- grep("^s_", names(licor_qc), value = TRUE)

list_flags <- function(df, cols) {
  m <- as.matrix(df[cols])
  apply(m, 1, function(r) paste(cols[r], collapse = "; "))
}

licor_qc$hard_reasons <- list_flags(licor_qc, hard_flags)
licor_qc$soft_reasons <- list_flags(licor_qc, soft_flags)

licor_qc <- licor_qc %>%
  mutate(qc_status = case_when(
    hard_reasons != "" ~ "exclude",
    soft_reasons != "" ~ "review",
    TRUE               ~ "ok"
  ))

# ---- 6. Summaries -----------------------------------------------------------

# How many rows are ok / review / exclude, by data source
licor_qc %>% count(WUE_source, qc_status) %>% print()

# Which flags fire, and whether they cluster on one machine
flag_summary <- licor_qc %>%
  group_by(WUE_source) %>%
  summarise(across(all_of(c(hard_flags, soft_flags)), sum), .groups = "drop") %>%
  pivot_longer(-WUE_source, names_to = "flag", values_to = "n") %>%
  pivot_wider(names_from = WUE_source, values_from = n) %>%
  arrange(flag)
print(flag_summary, n = Inf)

# Rows to inspect by hand
review_table <- licor_qc %>%
  filter(qc_status != "ok") %>%
  select(Plant_ID, Sp_Code, Plot, Date, obs, WUE_source, qc_status,
         total.leaf.area, A, A_c, gsw_c, Ci_c, CiCa_c, WUE_c,
         hard_reasons, soft_reasons) %>%
  arrange(desc(qc_status), Plant_ID, obs)
print(review_table, n = 50)
write.csv(review_table, here::here("data/licor_review_table.csv"), row.names = FALSE)

# ---- 7. Clean dataset -------------------------------------------------------
# Drops hard-flagged rows and masks iWUE where gsw is too low for the ratio
# to be meaningful (A and gsw themselves are kept). Soft-flagged rows stay in;
# decide on those after reviewing, and run analyses with and without them.

licor_clean <- licor_qc %>%
  filter(qc_status != "exclude") %>%
  mutate(WUE_clean = if_else(gsw_c < thr$gsw_min_wue, NA_real_, WUE_c))

write.csv(licor_clean,here::here("data/licor_clean.csv"))

# ---- 8. Quick visual check --------------------------------------------------

ggplot(licor_qc, aes(gsw_c, A_c, colour = qc_status, shape = WUE_source)) +
  geom_point(alpha = 0.7) +
  scale_colour_manual(values = c(ok = "grey50", review = "orange",
                                 exclude = "red")) +
  labs(x = expression(g[sw]~"(mol m"^-2~"s"^-1*", area-corrected)"),
       y = expression(A~"(µmol m"^-2~"s"^-1*", area-corrected)"),
       colour = "QC status", shape = "Source") +
  theme_bw()
