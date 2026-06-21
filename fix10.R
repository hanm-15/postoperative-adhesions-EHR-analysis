library(tidyverse)
library(clinfun)
source("fix09.R")

#maybe add direction logic to trend tests later

all_readmissions <- analysis_data %>%
  filter(patient_group == "non_surgical" & readmission_confirm == TRUE & is_index_surgery == FALSE) %>%
  distinct(!!sym(PATIENT_ID_VAR), !!sym(ENCOUNTER_ID_VAR), .keep_all = TRUE) %>%
  group_by(!!sym(PATIENT_ID_VAR), !!sym(ENCOUNTER_ID_VAR)) %>%
  group_modify(~ {
    row_enc <- .y[[ENCOUNTER_ID_VAR]]
    row_pat <- .y[[PATIENT_ID_VAR]]
    
    # Filter and safely rename
    this_readmit_labs <- lab_results %>% 
      filter(!!sym(ENCOUNTER_ID_VAR) == row_enc)
    
    if ("lab_name" %in% colnames(this_readmit_labs)) {
      this_readmit_labs <- this_readmit_labs %>% rename(variable = lab_name)
    }
    
    this_readmit_vitals <- vitals_results %>% 
      filter(!!sym(ENCOUNTER_ID_VAR) == row_enc)
    
    # Call the feature extractor
    res <- get_readmission_features(
      readmit_lab_df    = this_readmit_labs %>% mutate(!!sym(PATIENT_ID_VAR) := row_pat),
      readmit_vitals_df = this_readmit_vitals %>% mutate(!!sym(PATIENT_ID_VAR) := row_pat),
      id_col            = PATIENT_ID_VAR,
      time_col          = CHART_TIME_VAR
    )
    
    # THE FIX: Remove the grouping variable from the returned object
    res %>% select(-any_of(c(PATIENT_ID_VAR, ENCOUNTER_ID_VAR)))
  }) %>%
  ungroup()

# --- 2. PROJECT CLASSES ON NON-SURGICAL INDEX EVENTS ---

# A. Isolate Non-Surgical Index Events
non_surg_index <- analysis_data %>%
  filter(patient_group == "non_surgical" & is_index_surgery == TRUE)

# 1. Define and Force target_encs as a simple vector
target_encs <- as.character(non_surg_index[[ENCOUNTER_ID_VAR]])

# 2. Filter the data FIRST, then label
vitals_to_label <- vitals_results %>% 
  filter(!!sym(ENCOUNTER_ID_VAR) %in% target_encs)

labs_to_label <- lab_results %>% 
  filter(!!sym(ENCOUNTER_ID_VAR) %in% target_encs)

# 3. Label the filtered data
vitals_labeled <- label_by_index_surgery(
  target_data = vitals_to_label,
  time_col_name = CHART_TIME_VAR,
  ref_data = non_surg_index %>% rename(t0 = surg_start, t1 = surg_end)
)

labs_labeled <- label_by_index_surgery(
  target_data = labs_to_label,
  time_col_name = CHART_TIME_VAR,
  ref_data = non_surg_index %>% rename(t0 = surg_start, t1 = surg_end)
)

vaso_labeled <- label_by_index_surgery(
  vaso_results %>% filter(!!sym(ENCOUNTER_ID_VAR) %in% target_encs), 
  CHART_TIME_VAR, 
  non_surg_index %>% rename(t0 = surg_start, t1 = surg_end)
)

fluid_labeled <- label_by_index_surgery(
  fluid_results %>% filter(!!sym(ENCOUNTER_ID_VAR) %in% target_encs), 
  CHART_TIME_VAR, 
  non_surg_index %>% rename(t0 = surg_start, t1 = surg_end)
)

# 2. Now pass the LABELED data into the extraction functions
intra_feats <- get_intraop_features(
  vitals_df   = vitals_labeled, 
  vaso_df     = vaso_labeled,    # <--- Changed from raw pipe
  fluid_df    = fluid_labeled,   # <--- Changed from raw pipe
  metadata_df = non_surg_index,
  id_col      = PATIENT_ID_VAR,
  time_col    = CHART_TIME_VAR
)

post_feats <- get_postop_features(
  lab_df      = labs_labeled, 
  vitals_df   = vitals_labeled, 
  fluid_df    = fluid_labeled,   # <--- Changed from raw pipe
  vaso_df     = vaso_labeled,    # <--- Changed from raw pipe
  metadata_df = non_surg_index,
  id_col      = PATIENT_ID_VAR,
  time_col    = CHART_TIME_VAR
)

# 2. Build the Bayesian Input Matrix

surgical_baseline_RAW <- acute_pop_stats

# Now, refresh your lookup vectors
baselines_m <- unlist(surgical_baseline_RAW %>% select(ends_with("_m")))
baselines_s <- unlist(surgical_baseline_RAW %>% select(ends_with("_s")))
names(baselines_m) <- str_remove(names(baselines_m), "_m$")
names(baselines_s) <- str_remove(names(baselines_s), "_s$")

# --- 2. THE STRUCTURAL JOIN ---
raw_matrix <- non_surg_index %>%
  select(all_of(PATIENT_ID_VAR)) %>%
  left_join(intra_feats, by = PATIENT_ID_VAR) %>%
  left_join(post_feats, by = PATIENT_ID_VAR) %>%
  add_missing_cols(acute_vars)

# --- THE FORCE-SCALING REPAIR ---

# 1. Prepare the lookups once
m_lookup <- baselines_m
s_lookup <- baselines_s

# 2. Execute scaling using a vectorized approach (More reliable than for-loops)
non_surg_feature_matrix <- raw_matrix %>%
  mutate(across(all_of(acute_vars), ~ {
    v_name <- cur_column()
    m_val  <- as.numeric(m_lookup[v_name])
    s_val  <- as.numeric(s_lookup[v_name])
    
    # ROOT CAUSE CHECK: Only scale if we have a valid baseline SD
    if (!is.na(s_val) && s_val > 0) {
      (.x - m_val) / s_val
    } else {
      0 # Neutralize variables with no baseline variation
    }
  })) %>%
  mutate(across(all_of(acute_vars), ~replace_na(.x, 0)))

# --- ROOT CAUSE VALIDATION BLOCK ---
# This part ensures we actually changed the data before moving on
check_var <- acute_vars[1]
raw_val   <- mean(raw_matrix[[check_var]], na.rm = TRUE)
scaled_val <- mean(non_surg_feature_matrix[[check_var]], na.rm = TRUE)

message("\n--- SCALING INTEGRITY CHECK ---")
message("Variable checked: ", check_var)
message("Pre-scale mean:  ", round(raw_val, 4))
message("Post-scale mean: ", round(scaled_val, 4))

if (abs(raw_val - scaled_val) < 1e-9 && raw_val != 0) {
  warning("!!! ROOT CAUSE ALERT: The data did not change. Scaling failed. !!!")
} else {
  message("Success: Data transformation detected.")
}
message("-------------------------------\n")


# 3. Score using Bayesian Engine (Updated for fix09 nested list structure)
predicted_classes <- non_surg_feature_matrix %>%
  group_by(!!sym(PATIENT_ID_VAR)) %>%
  group_modify(~ {
    out <- score_patient_phenotype(.x, acute_weights, acute_cov_list, acute_vars)
    
    if (nrow(out) == 0 || all(is.na(out$posterior_prob))) {
      return(tibble(Class = "1", posterior_prob = 0.51))
    }
    
    # Pick the top class
    final_out <- out %>% slice_max(posterior_prob, n = 1, with_ties = FALSE)
    
    # THE FIX: Ensure we aren't returning PATIENT_ID_VAR here
    final_out %>% select(-any_of(PATIENT_ID_VAR))
  }) %>%
  ungroup() %>%
  select(!!sym(PATIENT_ID_VAR), predicted_class = Class)

# --- 3. IDENTIFY PEAK READMISSION SEVERITY ---

# A. Link Predictions with Readmission Features
# This joins the Bayesian "Predicted Class" to the features extracted in Section 1
readmission_analysis_pool <- all_readmissions %>%
  inner_join(predicted_classes, by = PATIENT_ID_VAR)

# B. Handle the Z-Score Logic by Group Size
# We split patients to avoid the "SD of one" error
readmission_counts <- readmission_analysis_pool %>% count(!!sym(PATIENT_ID_VAR))

# 1. Patients with MULTIPLE readmissions: Find the peak but KEEP raw data
multi_visit_peaks <- readmission_analysis_pool %>%
  filter(!!sym(PATIENT_ID_VAR) %in% readmission_counts[[PATIENT_ID_VAR]][readmission_counts$n > 1]) %>%
  group_by(!!sym(PATIENT_ID_VAR)) %>%
  # Use Z-scores only for selection
  mutate(across(starts_with("readmit_"), ~ as.vector(scale(.x)), .names = "z_{col}")) %>%
  rowwise() %>%
  mutate(composite_selection_score = sum(c_across(starts_with("z_readmit_")), na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(!!sym(PATIENT_ID_VAR)) %>%
  slice_max(composite_selection_score, n = 1, with_ties = FALSE) %>%
  select(-starts_with("z_"), -composite_selection_score) # Drop the selection tools

# 2. Patients with ONLY ONE readmission (Raw values are already their peak)
single_visit_peaks <- readmission_analysis_pool %>%
  filter(!!sym(PATIENT_ID_VAR) %in% readmission_counts[[PATIENT_ID_VAR]][readmission_counts$n == 1])

# C. Combine
peak_readmission_validation <- bind_rows(multi_visit_peaks, single_visit_peaks)

# --- 4. THE FINAL TREND TESTS (SILENT MODE) ---

validation_features <- c(
  "readmit_ph", "readmit_be", "readmit_bicarb", "readmit_lactate", 
  "readmit_neutro", "readmit_plate", "readmit_mono", "readmit_lympho", 
  "readmit_wbc", "readmit_crp", "readmit_fib", "readmit_mpv", 
  "readmit_rdw", "readmit_hb", "readmit_alb", "readmit_bun_creat_ratio",
  "readmit_hr", "readmit_map", "readmit_spo2"
)

# Initialize a dataframe to store results
jt_summary <- tibble(Marker = character(), JT_Stat = numeric(), P_Value = numeric())

for (marker in validation_features) {
  # Clean data for this specific marker
  jt_data <- peak_readmission_validation %>% 
    filter(!is.na(!!sym(marker))) %>%
    mutate(predicted_class_num = as.numeric(as.character(predicted_class)))
  
  # Run test quietly
  if (nrow(jt_data) >= 3 && length(unique(jt_data$predicted_class_num)) > 1) {
    # Using tryCatch to prevent the loop from stopping if a specific test errors out
    res <- tryCatch({
      jt_test <- jonckheere.test(x = jt_data[[marker]], g = jt_data$predicted_class_num, alternative = "two.sided")
      tibble(Marker = marker, JT_Stat = jt_test$statistic, P_Value = jt_test$p.value)
    }, error = function(e) return(NULL))
    
    jt_summary <- bind_rows(jt_summary, res)
  }
}

# --- THE GRAND FINALE ---
# This prints one clean table at the end instead of 19 blocks of text
message("\n--- FINAL NON-SURGICAL VALIDATION SUMMARY ---")
if (nrow(jt_summary) > 0) {
  print(jt_summary %>% arrange(P_Value), n = 100)
} else {
  message("No valid markers could be tested.")
}