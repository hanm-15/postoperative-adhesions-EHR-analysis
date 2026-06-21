library(tidyverse)
source("fix10.R")

# --- 1. THE CHRONIC DERIVATION LOOKUP (The Anchor) ---
# Manually mapping the three worlds to ensure 100% naming robustness.
chronic_der_lookup <- tibble(
  generic_name = c("ph", "be", "bicarb", "lactate", "neutro", "plate", "mono", "lympho", "wbc", 
                   "crp", "fib", "mpv", "rdw", "hb", "alb", "bun_creat_ratio", "hr", "map", "spo2"),
  preop_col    = c("pre_ph", "pre_be", "pre_bicarb", "pre_lactate", "pre_neutro", "pre_plate", "pre_mono", "pre_lympho", "pre_wbc", "pre_crp", "pre_fib", "pre_mpv", "pre_rdw", "pre_hb", "pre_alb", "pre_bun_creat_ratio", "pre_hr", "pre_map", "pre_spo2"),
  readmit_col  = c("readmit_ph", "readmit_be", "readmit_bicarb", "readmit_lactate", "readmit_neutro", "readmit_plate", "readmit_mono", "readmit_lympho", "readmit_wbc", "readmit_crp", "readmit_fib", "readmit_mpv", "readmit_rdw", "readmit_hb", "readmit_alb", "readmit_bun_creat_ratio", "readmit_hr", "readmit_map", "readmit_spo2"),
  bayesian_col = paste0(generic_name, "_mean")
)


# --- A. EXTRACT REOPERATED PRE-OP BASELINES ---
reop_chronic_baselines <- analysis_data %>%
  filter(patient_group == "reoperated" & lysis_confirm == TRUE) %>%
  inner_join(final_assignments_table %>% select(!!sym(PATIENT_ID_VAR), Class), by = PATIENT_ID_VAR) %>%
  group_by(!!sym(PATIENT_ID_VAR), !!sym(ENCOUNTER_ID_VAR)) %>%
  do({
    # THE FIX: Ensure these are single values (scalars), not vectors
    row_enc   <- .[[ENCOUNTER_ID_VAR]][1]
    row_pat   <- .[[PATIENT_ID_VAR]][1] 
    t0_stable <- .[[SURG_START_VAR]][1] 
    
    # Labs
    this_preop_labs <- lab_results %>% 
      filter(!!sym(ENCOUNTER_ID_VAR) == row_enc) %>%
      filter(!!sym(CHART_TIME_VAR) < t0_stable) %>%
      mutate(period = "pre_op", !!sym(PATIENT_ID_VAR) := row_pat) %>%
      group_by(variable) %>%
      slice_min(!!sym(CHART_TIME_VAR), n = 1, with_ties = FALSE) %>%
      ungroup()
    
    # Vitals
    this_preop_vitals <- vitals_results %>% 
      filter(!!sym(ENCOUNTER_ID_VAR) == row_enc) %>%
      filter(!!sym(CHART_TIME_VAR) < t0_stable) %>%
      mutate(period = "pre_op", !!sym(PATIENT_ID_VAR) := row_pat) %>%
      rename(!!sym(MAP_VAR) := any_of("map_val"), 
             !!sym(HR_VAR)  := any_of("hr_val"), 
             !!sym(SPO2_VAR) := any_of("spo2_val")) %>%
      slice_min(!!sym(CHART_TIME_VAR), n = 1, with_ties = FALSE) %>%
      ungroup()
    
    get_baseline_features(this_preop_labs, this_preop_vitals, PATIENT_ID_VAR, CHART_TIME_VAR)
  }) %>%
  ungroup() %>%
  # Add Class back
  left_join(final_assignments_table %>% select(!!sym(PATIENT_ID_VAR), Class), by = PATIENT_ID_VAR)

if (nrow(reop_chronic_baselines) == 0) {
  stop("ERROR: No baseline features were extracted. Check if 'period == pre_op' exists in labs_labeled for these encounters.")
}


# --- 2. GENERATE COMPARISON MATRIX (SCALED) ---
safe_scale <- function(x) {
  if(all(is.na(x)) || sd(x, na.rm = TRUE) == 0) {
    return(rep(0, length(x)))
  }
  as.vector(scale(x))
}

# --- 2. GENERATE COMPARISON MATRIX (FIXED: SHARED YARDSTICK) ---

# A. Capture the "Chronic Yardstick" from Reoperated Patients
# This is our gold standard for what "Chronic" looks like.
chronic_ref_stats <- reop_chronic_baselines %>%
  summarise(across(all_of(chronic_der_lookup$preop_col), 
                   list(m = ~mean(.x, na.rm = TRUE), s = ~sd(.x, na.rm = TRUE))))

# B. Standardize Reoperated Medians (Against their own population)
reop_stats <- reop_chronic_baselines %>%
  pivot_longer(cols = all_of(chronic_der_lookup$preop_col), names_to = "preop_col", values_to = "raw_reop") %>%
  group_by(Class, preop_col) %>%
  summarise(val_reop = median(raw_reop, na.rm = TRUE), .groups = "drop") %>%
  rowwise() %>%
  mutate(
    m_ref = chronic_ref_stats[[paste0(preop_col, "_m")]],
    s_ref = chronic_ref_stats[[paste0(preop_col, "_s")]],
    z_reop = if_else(s_ref > 0, (val_reop - m_ref) / s_ref, 0)
  ) %>%
  left_join(chronic_der_lookup, by = "preop_col")

# C. Standardize Non-Surgical Medians (AGAINST THE CHRONIC YARDSTICK)
# This is the critical change: measuring Readmissions by the Chronic scale.
nonsurg_stats <- peak_readmission_validation %>%
  rename(Class = predicted_class) %>%
  mutate(Class = as.character(Class)) %>%
  pivot_longer(cols = all_of(chronic_der_lookup$readmit_col), names_to = "readmit_col", values_to = "val_ns") %>%
  group_by(Class, readmit_col) %>%
  summarise(val_ns = median(val_ns, na.rm = TRUE), .groups = "drop") %>%
  left_join(chronic_der_lookup, by = "readmit_col") %>%
  rowwise() %>%
  mutate(
    # We use the same m_ref and s_ref from the Chronic population!
    m_ref = chronic_ref_stats[[paste0(preop_col, "_m")]], 
    s_ref = chronic_ref_stats[[paste0(preop_col, "_s")]],
    z_ns = if_else(s_ref > 0, (val_ns - m_ref) / s_ref, 0)
  ) %>%
  ungroup()

# --- 3. THE "NUDGE" BALANCING ACT ---
comparison_results <- inner_join(
  reop_stats, # This has z_reop from your previous block
  nonsurg_stats %>% select(Class, generic_name, z_ns, val_ns), # <--- NOW val_ns DEFINITELY EXISTS
  by = c("Class", "generic_name")
) %>%
  mutate(
    delta_z = abs(z_reop - z_ns),
    consistency_factor = 1 / (1 + delta_z)
  )

# --- 4. FINAL WEIGHT CALCULATION ---
# Get the RAW Reop means to apply the nudge
reop_raw_means <- reop_chronic_baselines %>%
  group_by(Class) %>%
  summarise(across(all_of(chronic_der_lookup$preop_col), ~mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  pivot_longer(-Class, names_to = "preop_col", values_to = "raw_reop_mean")

chronic_weight_set <- reop_raw_means %>%
  left_join(comparison_results, by = c("Class", "preop_col")) %>%
  mutate(
    # Use the 'preop_col' name (e.g., "pre_ph") directly to ensure prefixes are kept
    chronic_val = val_ns + (raw_reop_mean - val_ns) * consistency_factor,
    final_col_name = paste0(preop_col, "_mean") # Result: "pre_ph_mean"
  ) %>%
  select(Class, final_col_name, chronic_val) %>%
  pivot_wider(names_from = final_col_name, values_from = chronic_val)

# --- 4. BIOLOGICAL REALITY CHECK (Fixed Column Name) ---
message("\n--- TOP 5 MARKERS WITH HIGHEST CHRONIC PENALTY ---")
# Use 'delta_z' instead of 'delta'
print(comparison_results %>% 
        arrange(desc(delta_z)) %>% 
        select(Class, generic_name, delta_z, consistency_factor) %>% 
        head(5))

# --- PEEK AT THE FINAL WEIGHTS (Corrected for _mean suffix) ---
message("\n--- PREVIEW OF BALANCED WEIGHTS (Human Units) ---")

# Dynamically select the available 'pre_' mean columns
available_peek <- intersect(
  c("pre_ph_mean", "pre_lactate_mean", "pre_map_mean"), 
  names(chronic_weight_set)
)

if(length(available_peek) > 0) {
  print(chronic_weight_set %>% select(Class, all_of(available_peek)))
} else {
  message("Note: Standard peek columns not found. Printing first 4 columns instead:")
  print(chronic_weight_set[, 1:min(4, ncol(chronic_weight_set))])
}


# --- FINAL CALIBRATION: SCALE WEIGHTS AGAINST THE DERIVED POPULATION ---
message("\n--- FINAL CALIBRATION: STANDARDIZING AGAINST DERIVED CHRONIC POPULATION ---")

# 1. Use the data you JUST built as the absolute baseline
chronic_baseline_stats <- reop_chronic_baselines %>%
  summarise(across(all_of(chronic_der_lookup$preop_col), 
                   list(m = ~mean(.x, na.rm=TRUE), s = ~sd(.x, na.rm=TRUE))))

# 2. Map the balanced 'chronic_val' to this baseline's Z-scores
chronic_weight_set_final <- chronic_weight_set %>%
  pivot_longer(-Class, names_to = "final_col_name", values_to = "raw_val") %>%
  mutate(generic_preop = str_remove(final_col_name, "_mean$")) %>%
  rowwise() %>%
  mutate(
    m_val = chronic_baseline_stats[[paste0(generic_preop, "_m")]],
    s_val = chronic_baseline_stats[[paste0(generic_preop, "_s")]],
    # Calculate Z-score: (Raw - Population Mean) / Population SD
    z_val = if_else(!is.na(s_val) && s_val > 0, (raw_val - m_val) / s_val, 0)
  ) %>%
  ungroup() %>%
  select(Class, final_col_name, z_val) %>%
  pivot_wider(names_from = final_col_name, values_from = z_val)

# 3. Overwrite for LTA readiness
chronic_weight_set <- chronic_weight_set_final

# Check if both weight sets have the same number of features
n_acute   <- length(acute_vars)
n_chronic <- ncol(chronic_weight_set) - 1 # subtracting 'Class' column

message("Acute Features: ", n_acute, " | Chronic Features: ", n_chronic)

# Compare the Z-score range (Should generally be between -3 and 3)
summary(as.vector(as.matrix(chronic_weight_set[,-1])))

# Do the classes actually look different in their Chronic Baselines?
chronic_weight_set %>%
  select(Class, any_of(c("pre_lactate_mean", "pre_map_mean", "pre_ph_mean"))) %>%
  print()

# Verify that markers with HIGH delta_z have LOW consistency factors
comparison_results %>%
  select(generic_name, delta_z, consistency_factor) %>%
  arrange(desc(delta_z)) %>%
  head(10)

message("SUCCESS: Chronic weights are now standardized against the pre-op cohort.")