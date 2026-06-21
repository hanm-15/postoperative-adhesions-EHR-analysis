source("17.R")
library(tidyverse)
library(stats)
library(corpcor)
library(ggalluvial)
library(reshape2)

# --- 1. Scoring the EXTERNAL cohort ---
# Use 'lta_ext_scaled' instead of 'lta_input_data'

chronic_scores_ext <- score_patient_lta(lta_ext_scaled, chronic_weight_set, chronic_cov_list, chronic_vars_final, "incoming")
anchor_scores_ext  <- score_patient_lta(lta_ext_scaled, anchor_weights, anchor_cov_list, anchor_vars_final, "anchor")
acute_scores_ext   <- score_patient_lta(lta_ext_scaled, acute_weights, acute_cov_list, acute_vars_final, "outgoing")

# --- 2. Assemble the Transition Records ---
lta_transition_records_ext <- lta_skeleton %>%
  select(!!sym(PATIENT_ID_VAR), !!sym(ENCOUNTER_ID_VAR), node_sequence, node_label) %>%
  left_join(lta_weight_mapping %>% select(!!sym(ENCOUNTER_ID_VAR), weight_sets_to_use), by = ENCOUNTER_ID_VAR) %>%
  left_join(chronic_scores_ext, by = ENCOUNTER_ID_VAR) %>% 
  left_join(acute_scores_ext, by = ENCOUNTER_ID_VAR) %>% 
  left_join(anchor_scores_ext, by = ENCOUNTER_ID_VAR) %>% 
  mutate(
    anchor_class = if_else(map_lgl(weight_sets_to_use, ~"anchor_weights" %in% .x), anchor_class, NA_character_),
    outgoing_class = if_else(map_lgl(weight_sets_to_use, ~"acute_weights" %in% .x), outgoing_class, NA_character_),
    node_profile = paste0("Chr:", coalesce(incoming_class, "NA"), " Anc:", coalesce(anchor_class, "NA"), " Acu:", coalesce(outgoing_class, "NA"))
  )

# --- 4. THE JOURNEY ---
lta_patient_paths <- lta_transition_records_ext %>%
  group_by(pat_id) %>%
  arrange(node_sequence) %>%
  summarise(journey = paste(node_profile, collapse = " -> "), .groups = "drop")
lta_transition_records_ext <- lta_transition_records_ext %>%
  mutate(directionality = case_when(
    !is.na(incoming_class) & !is.na(outgoing_class) ~ paste0(incoming_class, " -> ", outgoing_class),
    !is.na(incoming_class) & is.na(outgoing_class) ~ paste0(incoming_class, " -> End"),
    TRUE ~ NA_character_
  ))
calculate_entropy <- function(probs) { mean(probs, na.rm = TRUE) }
message("Precision Audit by Encounter Type:")
lta_transition_records_ext %>%
  group_by(node_label) %>%
  summarise(
    N = n(),
    avg_incoming_prob = mean(incoming_prob, na.rm = TRUE),
    avg_anchor_prob = mean(anchor_prob, na.rm = TRUE),
    avg_outgoing_prob = mean(outgoing_prob, na.rm = TRUE)
  )
message("\n--- TRAJECTORY & TRANSITION AUDIT ---")
lta_transition_records_ext %>%
  filter(node_label != "readmission_confirm") %>%
  count(incoming_class, outgoing_class) %>%
  group_by(incoming_class) %>% 
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  print(n = 5)
message("\n--- FIRST 5 TRIPLE-WINDOW JOURNEYS ---")
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
lta_transition_records_ext %>% 
  select(all_of(c(ENCOUNTER_ID_VAR, "incoming_class", "incoming_prob", "anchor_class", "anchor_prob", "outgoing_class", "outgoing_prob"))) %>% 
  head(3) %>% print()


# --- 5. BIOLOGICAL VALIDATION (MULTI-SURGERY CHAIN) ---
lta_validation_ext <- lta_transition_records_ext %>%
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
summary_stats_ext <- lta_validation_ext %>%
  summarise(
    n_comparisons = n(),
    brier_chronic = mean(chronic_error, na.rm = TRUE),
    brier_predictive = mean(predictive_error, na.rm = TRUE)
  )

print(summary_stats_ext)



# --- 2. COVARIATE INTEGRATION (The New Function Call) ---
final_lta_analysis_ext <- extract_lta_covariates(
  skeleton_df = lta_skeleton, 
  raw_clinical_df = raw_surgery_ext,
  transition_records_df = lta_transition_records_ext
)

# --- 3. VALIDATION & VISUALIZATION (Portable Audit) ---

audit_lta_covariates(final_lta_analysis_ext)
visualize_lta_covariates(final_lta_analysis_ext)