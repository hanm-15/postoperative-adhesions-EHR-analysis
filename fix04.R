library(tidyverse)
library(lubridate)
source("fix02.R")

# check if the "first lysis" of everyone in the reoperated subset actually works as intended
# reconsider resection/excision and repair codes
# map variables later
# check line 112

# --- DEFINE THIS AT THE TOP OF fix04.R ---
required_raw_vars <- c(
  PCS_VAR, DX_VAR, DURATION_VAR, WEIGHT_VAR, EBL_VAR, 
  CRYSTALLOID_VAR, COLLOID_VAR, BLOOD_PROD_VAR, 
  IRRIGATION_VAR, IRRIGATION_IN_VAR, URINE_VAR,
  DOPE_VAR, DOBU_VAR, EPI_VAR, NORE_VAR, MIL_VAR, VASO_VAR
)

drug_vars <- c(DOPE_VAR, DOBU_VAR, EPI_VAR, NORE_VAR, MIL_VAR, VASO_VAR)


# --- 0. LCA Specific Prefixes and Other Processing ---
# These should eventually live in 01.R, but we'll define them here for now
resection_prefix <- "^0[DFU7][BT].*"
repair_prefix    <- "^0[DFU7]Q.*"

# 1. Summarize Meds (Find the MAX dose given during the encounter)
vaso_summary <- vaso_results %>%
  mutate(enc_id = as.character(enc_id)) %>%
  group_by(enc_id, variable) %>%
  summarize(max_dose = max(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = variable, values_from = max_dose)

# 2. Summarize Fluids (Find the TOTAL volume given/out)
# We use the mapping from your globals (URINE_VAR, etc.)
fluid_summary <- fluid_results %>%
  mutate(enc_id = as.character(enc_id)) %>%
  group_by(enc_id, variable) %>%
  summarize(total_vol = sum(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = variable, values_from = total_vol)

# vitals_data is the raw high-frequency table (1 or 5 min pings)
vitals_summary <- vitals_results %>%
  # --- THE FIX: Force ID to be Character ---
  mutate(across(all_of(ENCOUNTER_ID_VAR), as.character)) %>% 
  group_by(!!sym(ENCOUNTER_ID_VAR)) %>%
  summarize(
    calc_map_cv = sd(!!sym(MAP_VAR), na.rm = TRUE) / mean(!!sym(MAP_VAR), na.rm = TRUE),
    calc_hr_cv   = sd(!!sym(HR_VAR), na.rm = TRUE) / mean(!!sym(HR_VAR), na.rm = TRUE)
  ) %>%
  mutate(across(c(calc_map_cv, calc_hr_cv), ~replace_na(.x, 0)))

# --- 1. Isolate the LCA Population ---
lca_sample <- analysis_data %>%
  filter(reoperated_subset == TRUE) %>% 
  filter(lysis_confirm == TRUE) %>%
  group_by(!!sym(PATIENT_ID_VAR)) %>%
  arrange(!!sym(DATE_VAR)) %>%
  slice(1) %>% 
  ungroup()

# --- 2. Feature Extraction & Calculation (Safe Version) ---

# UPDATED HELPER: More robust column injection
add_missing_cols <- function(df, cols) {
  for (col in cols) {
    if (!col %in% names(df)) {
      df[[col]] <- NA_real_ # Force it to be a numeric NA
    }
  }
  return(df)
}

# 1. Prepare the Canvas: Join all clinical summaries
lca_prepared <- lca_sample %>%
  mutate(across(all_of(ENCOUNTER_ID_VAR), as.character)) %>%
  left_join(vitals_summary, by = ENCOUNTER_ID_VAR) %>%
  left_join(vaso_summary,   by = ENCOUNTER_ID_VAR) %>%
  left_join(fluid_summary,  by = ENCOUNTER_ID_VAR) %>%
  # Pass 1: Add all required columns (defaults to NA)
  add_missing_cols(required_raw_vars) %>%
  add_missing_cols(c("calc_map_cv", "calc_hr_cv")) %>%
  # Pass 2: Only zero-out drugs for the VIS score
  mutate(across(any_of(drug_vars), ~replace_na(.x, 0)))

# 2. Now run the calculation
lca_features <- lca_prepared %>%
  mutate(
    feat_cryst_kg = if_else(coalesce(.data[[WEIGHT_VAR]], 0) > 0, 
                            coalesce(.data[[CRYSTALLOID_VAR]], 0) / .data[[WEIGHT_VAR]], 0),
    feat_resection = if_else(grepl(resection_prefix, .data[[PCS_VAR]]), 1, 0),
    feat_iatrogenic = if_else(grepl("T812", .data[[DX_VAR]]) | grepl(repair_prefix, .data[[PCS_VAR]]), 1, 0),
    
    feat_duration   = as.numeric(coalesce(.data[[DURATION_VAR]], 0)),
    feat_ebl_kg     = if_else(coalesce(.data[[WEIGHT_VAR]], 0) > 0, 
                              coalesce(.data[[EBL_VAR]], 0) / .data[[WEIGHT_VAR]], 0),
    
    total_in  = rowSums(across(all_of(c(CRYSTALLOID_VAR, COLLOID_VAR, BLOOD_PROD_VAR))), na.rm = TRUE),
    total_out = rowSums(across(all_of(c(URINE_VAR, EBL_VAR))), na.rm = TRUE),
    
    feat_net_fluid_kg = if_else(coalesce(.data[[WEIGHT_VAR]], 0) > 0, 
                                (total_in - total_out) / .data[[WEIGHT_VAR]], 0),
    
    feat_uop_kg_hr = if_else(coalesce(.data[[WEIGHT_VAR]], 0) > 0 & feat_duration > 0,
                             (coalesce(.data[[URINE_VAR]], 0) / .data[[WEIGHT_VAR]]) / (feat_duration / 60), 0),
    
    feat_vis_score = coalesce(.data[[DOPE_VAR]], 0) + coalesce(.data[[DOBU_VAR]], 0) + 
      (100 * coalesce(.data[[EPI_VAR]], 0)) + (100 * coalesce(.data[[NORE_VAR]], 0)) + 
      (10 * coalesce(.data[[MIL_VAR]], 0)) + (10000 * coalesce(.data[[VASO_VAR]], 0)),
    
    # Use .data[[]] here to be 100% explicit to the R compiler
    feat_map_cv = coalesce(.data[["calc_map_cv"]], 0),
    feat_hr_cv  = coalesce(.data[["calc_hr_cv"]], 0)
    )
  # --- 3. Final Feature Scrubbing (Dynamic Filter) ---
  
  # 1. First, identify "Viable" columns: 
  # They must have more than 1 unique value (LCA crashes on "flat" columns)
  # and they must be in your list of features to scale.
candidate_cols <- c("feat_duration", "feat_uop_kg_hr", "feat_resection", "feat_iatrogenic",
                    "feat_ebl_kg", "feat_cryst_kg", "feat_net_fluid_kg", "feat_vis_score", 
                    "feat_map_cv", "feat_hr_cv")

# Find which of those actually have data and variance in YOUR current lca_features
viable_cols <- lca_features %>%
  select(any_of(candidate_cols)) %>%
  summarise(across(everything(), ~ n_distinct(.x, na.rm = TRUE))) %>%
  pivot_longer(everything()) %>%
  filter(value > 1) %>%
  pull(name)

# 2. Finalize the Dataset: Keep ID + only the viable features
lca_final <- lca_features %>%
  select(all_of(PATIENT_ID_VAR), all_of(viable_cols)) %>%
  # Fill any sparse NAs with the median so mclust doesn't delete the whole row
  mutate(across(all_of(viable_cols), ~replace_na(.x, median(.x, na.rm = TRUE))))

# 3. Scale only the columns we actually kept
lca_final <- lca_final %>%
  mutate(across(all_of(viable_cols), ~as.vector(scale(.x))))

# 4. Prepare the matrix for the model
lca_matrix <- lca_final %>% select(all_of(viable_cols))

message("!!! SUCCESS: LCA Feature Matrix ready with ", nrow(lca_final), " patients.")
message("Kept ", length(viable_cols), " clinical features with sufficient variance: ", 
        paste(viable_cols, collapse=", "))