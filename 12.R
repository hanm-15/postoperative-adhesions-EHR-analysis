library(tidyverse)
library(lubridate)
source("11.R")

message("Weight Sets Loaded:")
message("- Chronic Weights: ", nrow(chronic_weight_set), " classes.")
message("- Acute Weights:   ", nrow(acute_weights), " classes.")
message("- Anchor Weights:  ", nrow(anchor_weights), " classes.")

# --- LONGITUDINAL TRANSITION NODES ---

# --- Encounter Definition ---
lta_skeleton <- raw_data %>%
  group_by(!!sym(ENCOUNTER_ID_VAR)) %>%
  mutate(
    intraperitoneal_confirm = grepl(intraperitoneal_surg, !!sym(PCS_VAR)),
    
    lysis_confirm = grepl(lysis_proc, !!sym(PCS_VAR)) & 
      grepl(adhesion_diag, !!sym(DX_VAR)),
    
    readmission_confirm = grepl(adhesion_diag, !!sym(DX_VAR)) & 
      (!grepl(surgical_intervention, !!sym(PCS_VAR)) | is.na(!!sym(PCS_VAR)))
  ) %>%
  ungroup() %>%
  mutate(across(ends_with("_confirm"), ~replace_na(.x, FALSE))) %>%
  mutate(
    node_label = case_when(
      lysis_confirm           ~ "lysis_confirm",
      intraperitoneal_confirm ~ "intraperitoneal_confirm",
      readmission_confirm      ~ "readmission_confirm",
      TRUE                    ~ "excl_node"
    )
  ) %>%
  filter(node_label != "excl_node") %>%
  group_by(!!sym(PATIENT_ID_VAR)) %>%
  # CRITICAL FIX: Use .data[[ ]] for both arrange and difftime
  arrange(.data[[ADMIT_TIME_VAR]], .by_group = TRUE) %>% 
  mutate(
    node_sequence = paste0("n", row_number()),
    days_between_nodes = as.numeric(difftime(
      .data[[ADMIT_TIME_VAR]], 
      lag(.data[[ADMIT_TIME_VAR]]), 
      units = "days"
    )),
    days_between_nodes = replace_na(days_between_nodes, 0)
  ) %>%
  ungroup()

# --- 6. GENERATE COVARIANCE LISTS (THE NAMING FIX) ---

# A. CHRONIC COV (Uses 'pre_' prefixes from Appendix 4.C)
chronic_cov_list <- map(as.character(unique(chronic_weight_set$Class)), ~{
  calculate_class_covariance(
    data = reop_chronic_baselines, 
    vars = chronic_der_lookup$preop_col, 
    class_label = .x
  )
}) %>% set_names(as.character(unique(chronic_weight_set$Class)))

# B. ANCHOR COV (Uses 'feat_' prefixes from Appendix 3.A)
anchor_data_prepared <- lca_features %>%
  left_join(final_assignments_table %>% select(!!sym(PATIENT_ID_VAR), Class), by = PATIENT_ID_VAR)

anchor_cov_list <- map(as.character(unique(anchor_weights$Class)), ~{
  calculate_class_covariance(anchor_data_prepared, viable_cols, .x)
}) %>% set_names(as.character(unique(anchor_weights$Class)))

patient_journeys_wide <- lta_skeleton %>%
  filter(node_label != "excl_node") %>%
  select(!!sym(PATIENT_ID_VAR), node_sequence, node_label) %>%
  pivot_wider(names_from = node_sequence, values_from = node_label)

message("--- LTA NODES DEFINED: N1, N2 Sequence established ---")
print(head(patient_journeys_wide))

# --- 7. REFINED WEIGHT-SET RULES (THE DUAL-ROLE MAPPING) ---
lta_weight_mapping <- lta_skeleton %>%
  mutate(
    weight_sets_to_use = case_when(
      node_label == "lysis_confirm"           ~ list(c("chronic_weight_set", "anchor_weights", "acute_weights")),
      node_label == "intraperitoneal_confirm" ~ list(c("chronic_weight_set", "acute_weights")),
      node_label == "readmission_confirm"     ~ list(c("chronic_weight_set")),
      TRUE                                    ~ list(NULL)
    )
  )


# --- 8. VECTORIZED FEATURE ENGINE (LOOP-FREE) ---

# A. Create a Period-Labeled long table for all high-freq data
label_all_periods <- function(raw_df, skeleton) {
  # 1. Prepare a clean skeleton with ONLY the timing info we need
  skel_timing <- skeleton %>% 
    select(all_of(c(ENCOUNTER_ID_VAR, ADMIT_TIME_VAR, "node_label")), 
           any_of(c(SURG_START_VAR, SURG_END_VAR))) %>%
    distinct()
  
  # 2. Join, but remove any existing timing columns from the raw_df first
  # This prevents the "surg_start.x" vs "surg_start.y" nightmare
  raw_df %>%
    select(-any_of(c(SURG_START_VAR, SURG_END_VAR, ADMIT_TIME_VAR, "node_label"))) %>% 
    inner_join(skel_timing, by = ENCOUNTER_ID_VAR) %>%
    mutate(
      t0 = if_else(!is.na(.data[[SURG_START_VAR]]), .data[[SURG_START_VAR]], .data[[ADMIT_TIME_VAR]]),
      t1 = if_else(!is.na(.data[[SURG_END_VAR]]), .data[[SURG_END_VAR]], t0),
      obs_time = .data[[CHART_TIME_VAR]]
    ) %>%
    mutate(
      period = case_when(
        obs_time < t0 ~ "pre_op",
        obs_time >= t0 & obs_time <= t1 ~ "intra_op",
        as.numeric(difftime(obs_time, t1, units="hours")) > 0 & 
          as.numeric(difftime(obs_time, t1, units="hours")) <= 24 ~ "post_24",
        as.numeric(difftime(obs_time, t1, units="hours")) > 24 & 
          as.numeric(difftime(obs_time, t1, units="hours")) <= 48 ~ "post_48",
        TRUE ~ "other"
      )
    ) %>%
    filter(period != "other")
}

message("Labeling data periods...")
labs_labeled   <- label_all_periods(lab_results, lta_skeleton)
vitals_labeled <- label_all_periods(vitals_results, lta_skeleton)
fluids_labeled <- label_all_periods(fluid_results, lta_skeleton)
vaso_labeled   <- label_all_periods(vaso_results, lta_skeleton)


# --- 8B. SHAPE ALIGNMENT (CRITICAL) ---
# Pivot labs WIDE once here so all lenses can see the columns
labs_wide <- labs_labeled %>%
  pivot_wider(names_from = variable, values_from = value)


# --- LENS 1: CHRONIC ---

# Path A: Surgical Chronic (NOW USING labs_wide)
chronic_surg <- get_baseline_features(
  labs_wide %>% filter(node_label != "readmission_confirm"), 
  vitals_labeled %>% filter(node_label != "readmission_confirm"), 
  ENCOUNTER_ID_VAR, CHART_TIME_VAR
)

# Path B: Readmission Chronic (Uses your specific readmission function)
chronic_readmit <- get_readmission_features(
  lab_results %>% filter(!!sym(ENCOUNTER_ID_VAR) %in% (lta_skeleton %>% filter(node_label == "readmission_confirm") %>% pull(!!sym(ENCOUNTER_ID_VAR)))),
  vitals_results %>% filter(!!sym(ENCOUNTER_ID_VAR) %in% (lta_skeleton %>% filter(node_label == "readmission_confirm") %>% pull(!!sym(ENCOUNTER_ID_VAR)))),
  ENCOUNTER_ID_VAR, CHART_TIME_VAR
) %>% 
  rename_with(~str_replace(., "readmit_", "pre_"), starts_with("readmit_"))

chronic_features <- bind_rows(chronic_surg, chronic_readmit)



# --- LENS 2: ANCHOR (Fortified with NA Placeholders) ---

# Helper to sum a column ONLY if it exists, otherwise return NA
safe_sum_na <- function(df, col_name) {
  if (col_name %in% colnames(df)) {
    return(sum(df[[col_name]], na.rm = TRUE))
  } else {
    return(NA_real_) # Keeps it as a numeric NA for matrix safety
  }
}

anchor_features <- vitals_labeled %>%
  filter(period == "intra_op", node_label == "lysis_confirm") %>%
  group_by(!!sym(ENCOUNTER_ID_VAR)) %>%
  summarize(
    calc_map_cv = if(n()>1) sd(!!sym(MAP_VAR), na.rm=T)/mean(!!sym(MAP_VAR), na.rm=T) else NA_real_,
    calc_hr_cv  = if(n()>1) sd(!!sym(HR_VAR), na.rm=T)/mean(!!sym(HR_VAR), na.rm=T) else NA_real_
  ) %>%
  left_join(
    # --- FIXED PIVOT HERE ---
    fluids_labeled %>% 
      filter(period == "intra_op") %>%
      # Pivot so variables like 'fluid_cryst_ml' become actual column headers
      pivot_wider(names_from = variable, values_from = value, values_fn = sum) %>% 
      group_by(!!sym(ENCOUNTER_ID_VAR)) %>%
      summarize(
        tmp_cryst = safe_sum_na(pick(everything()), "fluid_cryst_ml"),
        tmp_coll  = safe_sum_na(pick(everything()), "intra_colloid"),
        tmp_blood = safe_sum_na(pick(everything()), "intra_blood"),
        tmp_urine = safe_sum_na(pick(everything()), "urine_ml"),
        tmp_ebl   = safe_sum_na(pick(everything()), "intra_ebl"),
        .groups = "drop"
      ), by = ENCOUNTER_ID_VAR
  ) %>%
  left_join(
    vaso_labeled %>% 
      filter(period == "intra_op") %>%
      # 1. CRITICAL STEP: Move drug names from rows to columns
      pivot_wider(names_from = variable, 
                  values_from = value, 
                  values_fn = max) %>% 
      group_by(!!sym(ENCOUNTER_ID_VAR)) %>%
      # 2. Now safe_sum_na can find NORE_VAR and VASO_VAR
      summarize(
        tmp_vis = (1 * safe_sum_na(pick(everything()), DOPE_VAR)) + 
          (100 * safe_sum_na(pick(everything()), NORE_VAR)) + 
          (1 * safe_sum_na(pick(everything()), DOBU_VAR)) + 
          (100 * safe_sum_na(pick(everything()), EPI_VAR)) + 
          (10000 * safe_sum_na(pick(everything()), VASO_VAR)) + 
          (10 * safe_sum_na(pick(everything()), MIL_VAR)),
        .groups = "drop"
      ), by = ENCOUNTER_ID_VAR
  ) %>%
  left_join(lta_skeleton, by = ENCOUNTER_ID_VAR) %>%
  mutate(
    feat_resection  = if_else(grepl(resection_prefix, !!sym(PCS_VAR)), 1, 0),
    feat_iatrogenic = if_else(grepl("T812", !!sym(DX_VAR)) | grepl(repair_prefix, !!sym(PCS_VAR)), 1, 0),
    feat_duration   = as.numeric(coalesce(!!sym(DURATION_VAR), 0)),
    
    # Feature calculation: if tmp_ebl is NA, feat_ebl_kg becomes NA automatically
    feat_ebl_kg     = tmp_ebl / !!sym(WEIGHT_VAR),
    feat_cryst_kg   = tmp_cryst / !!sym(WEIGHT_VAR),
    
    # We use rowSums with na.rm=FALSE so if any component is missing, the total is NA
    # This signals to the LTA "Do not use this combined feature"
    total_in        = tmp_cryst + tmp_coll + tmp_blood,
    total_out       = tmp_urine + tmp_ebl,
    
    feat_net_fluid_kg = (total_in - total_out) / !!sym(WEIGHT_VAR),
    feat_uop_kg_hr    = (tmp_urine / !!sym(WEIGHT_VAR)) / (feat_duration / 60),
    feat_vis_score    = tmp_vis,
    feat_map_cv       = calc_map_cv,
    feat_hr_cv        = calc_hr_cv
  ) %>%
  select(!!sym(ENCOUNTER_ID_VAR), starts_with("feat_"))

# --- LENS 3: ACUTE (With Leak Protection) ---

# --- 1. Get Intraoperative Features ---
acute_intra <- get_intraop_features(
  vitals_df   = vitals_labeled %>% filter(node_label != "readmission_confirm"),
  vaso_df     = vaso_labeled,  # ensure_wide() inside 07.R handles this
  fluid_df    = fluids_labeled, # ensure_wide() inside 07.R handles this
  metadata_df = lta_skeleton,
  id_col      = ENCOUNTER_ID_VAR,
  time_col    = CHART_TIME_VAR
)

acute_post <- get_postop_features(
  lab_df      = labs_labeled %>% filter(node_label != "readmission_confirm"),
  vitals_df   = vitals_labeled %>% filter(node_label != "readmission_confirm"),
  fluid_df    = fluids_labeled,
  vaso_df     = vaso_labeled,
  metadata_df = lta_skeleton,
  id_col      = ENCOUNTER_ID_VAR,
  time_col    = CHART_TIME_VAR
)

acute_features <- acute_intra %>%
  left_join(acute_post, by = ENCOUNTER_ID_VAR) %>%
  mutate(across(where(is.numeric), ~replace_na(.x, 0)))



# --- 9. FINAL ASSEMBLY ---

lta_input_data <- lta_weight_mapping %>%
  left_join(chronic_features, by = ENCOUNTER_ID_VAR) %>%
  left_join(anchor_features, by = ENCOUNTER_ID_VAR) %>%
  left_join(acute_features, by = ENCOUNTER_ID_VAR)



message("--- LTA INPUT READY ---")
# Print Check
print(lta_input_data %>% 
        select(any_of(c(PATIENT_ID_VAR, "node_sequence", "node_label", "weight_sets_to_use"))) %>% 
        head())


# --- SECTION 9.B: PRESERVE RAW CLINICAL VALUES ---
# This creates a "Source of Truth" for audits and clinical tables
lta_input_copy <- lta_input_data 

message("--- RAW DATA PRESERVED AS 'lta_input_copy' ---")

# --- 10. FEATURE HEALTH CHECK ---

message("--- FEATURE DENSITY REPORT ---")

diagnostic_peek <- lta_input_data %>%
  group_by(node_label) %>%
  slice_head(n = 2) %>%
  ungroup() %>%
  arrange(node_label, pat_id)

print(diagnostic_peek %>% 
        select(pat_id, node_label, starts_with(c("pre_", "feat_", "post_"))) %>%
        select(1:8))


# --- 10. REFINED FEATURE HEALTH DASHBOARD ---
#message("Full feature glimpse")
#message("\n" , paste(rep("=", 50), collapse = ""))
#message(paste(rep("=", 50), collapse = ""), "\n")

# Select a representative set of columns from each "Lens"
#diagnostic_summary <- lta_input_data %>%
#  group_by(node_label) %>%
#  slice_head(n = 2) %>%
#  ungroup() %>%
#  select(
#    pat_id, node_label, starts_with(c("pre_", "feat_", "intra", "post_"))
#  ) %>%
#  mutate(across(where(is.numeric), ~round(.x, 3)))
#print.data.frame(diagnostic_summary, row.names = FALSE)

message("\n" , paste(rep("-", 50), collapse = ""))
message("VARIABLE COMPLETION BY NODE TYPE (%)")
lta_input_data %>%
  group_by(node_label) %>%
  summarise(
    chronic_pct = mean(!is.na(pre_ph) & pre_ph != 0) * 100,
    anchor_pct  = mean(!is.na(feat_duration)) * 100,
    acute_pct   = mean(!is.na(post_ph_slope) & post_ph_slope != 0) * 100,
    .groups = "drop"
  ) %>%
  print.data.frame(row.names = FALSE)


# --- 1. THE GLOBALIZED METADATA LIST ---
# Using your global variables to ensure ID columns are never scaled
#metadata_to_skip <- c(
#  ENCOUNTER_ID_VAR, PATIENT_ID_VAR, 
#  "Class", "profile", "predicted_class", "node_label", "node_sequence", 
#  "patient_group", "reoperated_subset", "weight_sets_to_use",
#  "is_index_surgery", "encounter_num", ADMIT_TIME_VAR, 
#  "discharge_time", SURG_START_VAR, SURG_END_VAR, "event_start",
#  PCS_VAR, DX_VAR, ASA_VAR, AGE_VAR, GENDER_VAR, HEIGHT_VAR)

feature_prefixes <- c("pre_", "feat_", "intra_", "post_")
metadata_to_skip <- colnames(lta_input_data)[!grepl(paste0("^(", paste(feature_prefixes, collapse="|"), ")"), colnames(lta_input_data))]


# --- 1. IDENTIFY CLINICAL FEATURES ---
clinical_feats <- lta_input_data %>%
  select(where(is.numeric)) %>%
  names() %>%
  setdiff(metadata_to_skip) %>%
  .[!grepl("_confirm$", .)]

# --- 6. DYNAMIC FEATURE LISTS (EXPLICIT ALIGNMENT) ---
# Since weights are now "pre_ph_mean", we just look for "pre_ph" in the data
chronic_vars_final <- names(chronic_weight_set) %>%
  str_subset("^pre_.*_mean$") %>%   # Look for the explicit pre_ prefix
  str_remove("_mean$") %>%          # Get the variable name (e.g., "pre_ph")
  intersect(names(lta_input_data)) %>%
  setdiff(metadata_to_skip)

# Anchor and Acute remain specific to their own prefixes (intra_/post_)
anchor_vars_final <- names(anchor_weights) %>%
  str_subset("_mean$") %>%
  str_remove("_mean$") %>%
  intersect(names(lta_input_data)) %>%
  setdiff(metadata_to_skip)
acute_vars_final <- names(acute_weights) %>%
  str_subset("_mean$") %>%
  str_remove("_mean$") %>%
  intersect(names(lta_input_data))
# --- FINAL VERIFICATION ---
message("Features explicitly aligned for LTA:")
message("- Chronic (Pre-op): ", length(chronic_vars_final))
message("- Anchor (Intra-op): ", length(anchor_vars_final))
message("- Acute (Post-op):   ", length(acute_vars_final))


# --- LTA INPUT ALIGNMENT (THE TRIPLE YARDSTICK) ---

lta_input_scaled <- lta_input_data %>%
  mutate(across(-any_of(c(metadata_to_skip, "weight_sets_to_use")), ~as.numeric(.x))) %>% 
  mutate(
    # 1. SCALE CHRONIC LENS
    across(all_of(chronic_vars_final), ~{
      m <- chronic_baseline_stats[[paste0(cur_column(), "_m")]]
      s <- chronic_baseline_stats[[paste0(cur_column(), "_s")]]
      if(!is.null(s) && !is.na(s) && s > 0) (.x - m) / s else 0
    }),
    
    # 2. SCALE ANCHOR LENS
    across(all_of(anchor_vars_final), ~{
      m <- anchor_pop_stats[[paste0(cur_column(), "_m")]]
      s <- anchor_pop_stats[[paste0(cur_column(), "_s")]]
      if(!is.null(s) && !is.na(s) && s > 0) (.x - m) / s else 0
    }),
    
    # 3. SCALE ACUTE LENS
    across(all_of(acute_vars_final), ~{
      m <- acute_pop_stats[[paste0(cur_column(), "_m")]]
      s <- acute_pop_stats[[paste0(cur_column(), "_s")]]
      if(!is.null(s) && !is.na(s) && s > 0) (.x - m) / s else 0
    })
  )
# Overwrite for File 13 readiness
lta_input_data <- lta_input_scaled