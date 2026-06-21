library(tidyverse)
library(stats)
library(corpcor)
library(ggalluvial)
library(reshape2)
source("12.R")

# Engine
score_patient_lta <- function(patient_data, weight_set, cov_list, feature_vars, lens_label) {
  classes <- names(cov_list)
  results <- map_df(classes, ~{
    current_class_name <- .x
    class_obj <- cov_list[[current_class_name]]
    if (is.null(class_obj) || !is.list(class_obj)) return(NULL)
    class_cov <- class_obj$mat
    class_center <- weight_set %>%
      filter(Class == current_class_name) %>%
      select(all_of(paste0(feature_vars, "_mean"))) %>% 
      as.matrix() %>% as.vector()
    patient_matrix <- patient_data %>%
      select(all_of(feature_vars)) %>%
      as.matrix()
    d2 <- tryCatch({
      stats::mahalanobis(patient_matrix, center = class_center, cov = class_cov)
    }, error = function(e) rep(NA_real_, nrow(patient_matrix)))
    tibble(
      !!sym(ENCOUNTER_ID_VAR) := patient_data[[ENCOUNTER_ID_VAR]],
      Class = current_class_name,
      likelihood = dchisq(d2, df = length(feature_vars))
    )
  })
  # --- Bayesian Normalization ---
  results %>%
    group_by(!!sym(ENCOUNTER_ID_VAR)) %>%
    mutate(
      total_l = sum(likelihood, na.rm = TRUE),
      prob = case_when(total_l == 0 ~ 0, TRUE ~ likelihood / total_l)
    ) %>%
    ungroup() %>%
    group_by(!!sym(ENCOUNTER_ID_VAR)) %>%
    slice_max(prob, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(
      !!sym(ENCOUNTER_ID_VAR), 
      !!paste0(lens_label, "_class") := Class, 
      !!paste0(lens_label, "_prob") := prob
    )
}
# 2. Sync Covariance Lists 
sync_cov <- function(cov_list, valid_vars) {
  map(cov_list, ~{
    clean_internal_vars <- str_remove(.x$vars, "_mean")
    idx <- which(clean_internal_vars %in% valid_vars)
    .x$mat <- .x$mat[idx, idx, drop = FALSE]
    .x$vars <- valid_vars[valid_vars %in% clean_internal_vars]
    .x
  })
}
chronic_cov_list <- sync_cov(chronic_cov_list, chronic_vars_final)
anchor_cov_list <- sync_cov(anchor_cov_list, anchor_vars_final)
acute_cov_list <- sync_cov(acute_cov_list, acute_vars_final)
# Generate the 3 score banks
chronic_scores <- score_patient_lta(lta_input_data, chronic_weight_set, chronic_cov_list, chronic_vars_final, "incoming")
anchor_scores <- score_patient_lta(lta_input_data, anchor_weights, anchor_cov_list, anchor_vars_final, "anchor")
acute_scores <- score_patient_lta(lta_input_data, acute_weights, acute_cov_list, acute_vars_final, "outgoing")
# --- 1. JOIN THE MAPPING FROM FILE 12 ---
lta_transition_records <- lta_skeleton %>%
  select(pat_id, enc_id, node_sequence, node_label) %>%
  left_join(lta_weight_mapping %>% select(enc_id, weight_sets_to_use), by = "enc_id") %>%
  left_join(chronic_scores, by = "enc_id") %>% 
  left_join(acute_scores, by = "enc_id") %>% 
  left_join(anchor_scores, by = "enc_id") %>% 
  mutate(
    anchor_class = if_else(map_lgl(weight_sets_to_use, ~"anchor_weights" %in% .x), anchor_class, NA_character_),
    outgoing_class = if_else(map_lgl(weight_sets_to_use, ~"acute_weights" %in% .x), outgoing_class, NA_character_),
    node_profile = paste0("Chr:", coalesce(incoming_class, "NA"), " Anc:", coalesce(anchor_class, "NA"), " Acu:", coalesce(outgoing_class, "NA"))
  )
# --- 4. THE JOURNEY ---
lta_patient_paths <- lta_transition_records %>%
  group_by(pat_id) %>%
  arrange(node_sequence) %>%
  summarise(journey = paste(node_profile, collapse = " -> "), .groups = "drop")
lta_transition_records <- lta_transition_records %>%
  mutate(directionality = case_when(
    !is.na(incoming_class) & !is.na(outgoing_class) ~ paste0(incoming_class, " -> ", outgoing_class),
    !is.na(incoming_class) & is.na(outgoing_class) ~ paste0(incoming_class, " -> End"),
    TRUE ~ NA_character_
  ))
calculate_entropy <- function(probs) { mean(probs, na.rm = TRUE) }
message("Precision Audit by Encounter Type:")
lta_transition_records %>%
  group_by(node_label) %>%
  summarise(
    N = n(),
    avg_incoming_prob = mean(incoming_prob, na.rm = TRUE),
    avg_anchor_prob = mean(anchor_prob, na.rm = TRUE),
    avg_outgoing_prob = mean(outgoing_prob, na.rm = TRUE)
  )
message("\n--- TRAJECTORY & TRANSITION AUDIT ---")
lta_transition_records %>%
  filter(node_label != "readmission_confirm") %>%
  count(incoming_class, outgoing_class) %>%
  group_by(incoming_class) %>% 
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  print(n = 5)
message("\n--- TOP 5 TRIPLE-WINDOW JOURNEYS ---")
lta_patient_paths %>% count(journey, sort = TRUE) %>% head(5) %>% print()

# --- REVISED CENTERS CHECK ---
message("\n--- CENTERS SAMPLE (Class 1 vs Class 2 Across Lenses) ---")
audit_centers <- function(df, label) {
  cat("\n---", label, "Lens ---\n")
  df %>%
    select(Class, matches("^(pre_|feat_|intra_|post_).*_mean$")) %>%
    group_by(Class) %>% slice(1) %>% ungroup() %>% print()
}
audit_centers(chronic_weight_set, "CHRONIC (Pre-op)")
audit_centers(anchor_weights, "ANCHOR (Intra-op)")
audit_centers(acute_weights, "ACUTE (Post-op)")
message("\n--- LIKELIHOOD SAMPLE ---")
lta_transition_records %>% 
  select(all_of(c(ENCOUNTER_ID_VAR, "incoming_class", "incoming_prob", "anchor_class", "anchor_prob", "outgoing_class", "outgoing_prob"))) %>% 
  head(3) %>% print()


# --- 5. BIOLOGICAL VALIDATION (MULTI-SURGERY CHAIN) ---
lta_validation <- lta_transition_records %>%
  group_by(pat_id) %>%
  arrange(node_sequence) %>%
  mutate(
    # 'predicted_class' is the Acute-Response from the PREVIOUS surgical exit
    predicted_class = lag(outgoing_class),
    predicted_prob  = lag(outgoing_prob)
  ) %>%
  # Handle back-to-back surgeries or intervening readmissions:
  # This ensures the prediction from the most recent surgery stays active 
  # until the next 'Truth' event (Lysis) occurs.
  fill(predicted_class, predicted_prob, .direction = "down") %>%
  
  # We ONLY evaluate the 'Truth' at Lysis events
  filter(node_label == "lysis_confirm") %>%
  filter(!is.na(anchor_class)) %>%
  ungroup() %>%
  mutate(
    # TARGET: The Ground Truth Anchor from the current row
    is_anchor_severe = if_else(anchor_class == "2", 1, 0),
    
    # CONCURRENT VALIDITY: (Current Chronic vs Current Anchor)
    # Does the pre-op 'look' match the intra-op 'reality'?
    chronic_error = if_else(incoming_class == anchor_class, 
                            (1 - incoming_prob)^2, 
                            (0 - incoming_prob)^2),
    
    # PREDICTIVE VALIDITY: (Previous Predicted vs Current Anchor)
    # Does the last surgery's recovery match this surgery's reality?
    predictive_error = if_else(predicted_class == anchor_class, 
                               (1 - predicted_prob)^2, 
                               (0 - predicted_prob)^2)
  )

# --- FINAL SUMMARY STATISTICS ---
message("\n--- MODEL CALIBRATION AUDIT ---")
summary_stats <- lta_validation %>%
  summarise(
    n_comparisons = n(),
    brier_chronic = mean(chronic_error, na.rm = TRUE),
    brier_predictive = mean(predictive_error, na.rm = TRUE)
  )

print(summary_stats)


# --- 5. BIOLOGICAL VALIDATION (MULTI-SURGERY CHAIN) ---
lta_validation <- lta_transition_records %>%
  group_by(pat_id) %>%
  arrange(node_sequence) %>%
  mutate(
    # 'predicted_class' is the Acute-Response from the PREVIOUS surgical exit
    predicted_class = lag(outgoing_class),
    predicted_prob  = lag(outgoing_prob)
  ) %>%
  # Handle back-to-back surgeries or intervening readmissions:
  # This ensures the prediction from the most recent surgery stays active 
  # until the next 'Truth' event (Lysis) occurs.
  fill(predicted_class, predicted_prob, .direction = "down") %>%
  
  # We ONLY evaluate the 'Truth' at Lysis events
  filter(node_label == "lysis_confirm") %>%
  filter(!is.na(anchor_class)) %>%
  ungroup() %>%
  mutate(
    # TARGET: The Ground Truth Anchor from the current row
    is_anchor_severe = if_else(anchor_class == "2", 1, 0),
    
    # CONCURRENT VALIDITY: (Current Chronic vs Current Anchor)
    # Does the pre-op 'look' match the intra-op 'reality'?
    chronic_error = if_else(incoming_class == anchor_class, 
                            (1 - incoming_prob)^2, 
                            (0 - incoming_prob)^2),
    
    # PREDICTIVE VALIDITY: (Previous Predicted vs Current Anchor)
    # Does the last surgery's recovery match this surgery's reality?
    predictive_error = if_else(predicted_class == anchor_class, 
                               (1 - predicted_prob)^2, 
                               (0 - predicted_prob)^2)
  )

# --- FINAL SUMMARY STATISTICS ---
message("\n--- MODEL CALIBRATION AUDIT ---")
summary_stats <- lta_validation %>%
  summarise(
    n_comparisons = n(),
    brier_chronic = mean(chronic_error, na.rm = TRUE),
    brier_predictive = mean(predictive_error, na.rm = TRUE)
  )

print(summary_stats)




# --- COVARIATE EXTRACTION ENGINE ---
extract_lta_covariates <- function(skeleton_df, raw_clinical_df, transition_records_df) {
  
  message("Extracting Patient-Level Statics...")
  # 1. COMORBIDITY & STATIC FEATURES (STRICTLY PATIENT-LEVEL)
  patient_statics <- raw_clinical_df %>%
    group_by(!!sym(PATIENT_ID_VAR)) %>%
    summarise(
      cov_diabetes = as.numeric(any(str_detect(!!sym(DX_VAR), "^E1[013]"))),
      cov_ckd      = as.numeric(any(str_detect(!!sym(DX_VAR), "^N18"))),
      cov_gender   = first(!!sym(GENDER_VAR)),
      .groups = "drop"
    )

  message("Calculating Encounter-Level Dynamics...")
  # 2. PROCEDURAL & DEMOGRAPHIC LOGIC (ENCOUNTER-LEVEL)
  lta_covariate_set <- skeleton_df %>%
    mutate(
      age_raw    = !!sym(AGE_VAR),
      weight_raw = !!sym(WEIGHT_VAR),
      height_raw = !!sym(HEIGHT_VAR),
      approach_char = substr(!!sym(PCS_VAR), 5, 5),
      approach_raw = case_when(
        approach_char == "0" ~ "Open",
        approach_char == "4" ~ "Laparoscopic",
        TRUE ~ "Other"
      ),
      urgency_raw = !!sym(URGENCY_VAR),
      asa_raw     = !!sym(ASA_VAR)
    ) %>%
    group_by(!!sym(PATIENT_ID_VAR)) %>%
    arrange(.data[[ADMIT_TIME_VAR]]) %>%
    mutate(
      cov_age      = age_raw,
      cov_bmi      = weight_raw / ((height_raw / 100)^2),
      cov_asa      = asa_raw,
      cov_last_approach = approach_raw,
      cov_last_urgency  = urgency_raw,
      cov_local_trauma    = cumsum(node_label %in% c("lysis_confirm", "intraperitoneal_confirm")),
      cov_systemic_trauma = cumsum(grepl("^.{4}[04].*", !!sym(PCS_VAR)))
    ) %>%
    fill(cov_age, cov_bmi, cov_asa, cov_last_approach, cov_last_urgency, .direction = "down") %>%
    ungroup() %>%
    mutate(days_between_nodes = replace_na(days_between_nodes, 0))

  message("Performing Master Join...")
  # 3. MASTER JOIN
  final_df <- transition_records_df %>%
    left_join(patient_statics, by = PATIENT_ID_VAR) %>%
    left_join(lta_covariate_set %>% 
                select(!!sym(ENCOUNTER_ID_VAR), starts_with("cov_"), days_between_nodes), 
              by = ENCOUNTER_ID_VAR) %>%
    mutate(
      cov_gender = factor(cov_gender),
      cov_asa    = factor(cov_asa, levels = c(1, 2, 3, 4, 5)),
      cov_last_approach = factor(cov_last_approach, levels = c("Laparoscopic", "Open", "Other")),
      cov_last_urgency  = factor(cov_last_urgency, levels = c("Elective", "Urgent", "Emergency"))
    )
  
  return(final_df)
}

# --- EXECUTION ---
final_lta_analysis <- extract_lta_covariates(lta_skeleton, raw_data, lta_transition_records)

message("\n--- GLOBAL COVARIATE AUDIT (ALL FEATURES) ---")

# --- FUNCTION 1: Console Audit (The "Table 1" Generator) ---
audit_lta_covariates <- function(analysis_df) {
  message("\n--- GLOBAL COVARIATE AUDIT (ALL FEATURES) ---")
  
  stats_table <- analysis_df %>%
    filter(!is.na(anchor_class)) %>%
    group_by(anchor_class) %>%
    summarise(
      N = n(),
      # Demographics
      Age = round(mean(cov_age, na.rm=TRUE), 1),
      BMI = round(mean(cov_bmi, na.rm=TRUE), 1),
      Pct_Male = round(mean(cov_gender == "Male", na.rm=TRUE)*100, 1),
      
      # Comorbidities
      Pct_Diabetes = round(mean(cov_diabetes, na.rm=TRUE)*100, 1),
      Pct_CKD = round(mean(cov_ckd, na.rm=TRUE)*100, 1),
      
      # Severity & Approach
      Avg_ASA = round(mean(as.numeric(as.character(cov_asa)), na.rm=TRUE), 2),
      Pct_Open = round(mean(cov_last_approach == "Open", na.rm=TRUE)*100, 1),
      Pct_Emergency = round(mean(cov_last_urgency == "Emergency", na.rm=TRUE)*100, 1),
      
      # Trauma & Timing
      Local_Trauma = round(mean(cov_local_trauma, na.rm=TRUE), 2),
      Systemic_Trauma = round(mean(cov_systemic_trauma, na.rm=TRUE), 2),
      Days_Since_Last = round(mean(days_between_nodes, na.rm=TRUE), 1),
      .groups = "drop"
    ) %>%
    t() # Transpose for readability
  
  print(stats_table)
}

# --- FUNCTION 2: Visualization Suite ---
visualize_lta_covariates <- function(analysis_df) {
  
  # 1. Continuous Comparison (Added ASA here)
  message("Generating Multi-Feature Comparison Plot...")
  plot_data <- analysis_df %>%
    filter(!is.na(anchor_class)) %>%
    mutate(cov_asa_num = as.numeric(as.character(cov_asa))) %>% # Convert factor to numeric for plot
    select(anchor_class, cov_age, cov_bmi, cov_asa_num, 
           cov_local_trauma, cov_systemic_trauma, days_between_nodes) %>%
    pivot_longer(cols = -anchor_class, names_to = "Feature", values_to = "Value")
  
  covariate_boxplots <- ggplot(plot_data, aes(x = anchor_class, y = Value, fill = anchor_class)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.1, size = 0.5) +
    facet_wrap(~Feature, scales = "free_y") + 
    theme_minimal() +
    scale_fill_brewer(palette = "Set1") +
    labs(title = "Clinical Feature Distribution by Anchor Class",
         subtitle = "Continuous variables including ASA and Trauma metrics",
         x = "LTA Anchor Class")
  
  print(covariate_boxplots)
  
  # 2. Categorical Heatmap (Added Gender/Male here)
  message("Generating Categorical Risk Heatmap...")
  cat_data <- analysis_df %>%
    filter(!is.na(anchor_class)) %>%
    group_by(anchor_class) %>%
    summarise(
      Diabetes = mean(cov_diabetes, na.rm=TRUE),
      CKD = mean(cov_ckd, na.rm=TRUE),
      Emergency = mean(cov_last_urgency == "Emergency", na.rm=TRUE),
      Open_Surg = mean(cov_last_approach == "Open", na.rm=TRUE),
      Male_Sex = mean(cov_gender == "Male", na.rm=TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(cols = -anchor_class)
  
  covariate_heatmap <- ggplot(cat_data, aes(x = anchor_class, y = name, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "white", high = "darkred", labels = scales::percent) +
    theme_minimal() +
    labs(title = "Risk Factor Prevalence by Phenotype",
         subtitle = "Includes Gender and Comorbidities",
         y = "Clinical Feature", x = "Anchor Class", fill = "Prevalence")
  
  print(covariate_heatmap)
}

# Visualize results
audit_lta_covariates(final_lta_analysis)
visualize_lta_covariates(final_lta_analysis)

message("Extracted Encounter-Level Covariates")