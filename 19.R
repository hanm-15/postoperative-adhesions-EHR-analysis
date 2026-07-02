source("18.R")
library(tidyverse)
library(ggplot2)
library(lme4)
library(lmerTest)

message("\n--- SECTION 5.3: CROSS-COHORT VALIDATION ---")

# 1. Helper function to extract Raw Means and 95% CIs per incoming_class
extract_incoming_class_parameters <- function(data, scores, vars, cohort_label) {
  # Join the RAW data (data) with the incoming_class assignments (scores)
  data %>%
    inner_join(scores %>% select(!!sym(ENCOUNTER_ID_VAR), incoming_class), by = ENCOUNTER_ID_VAR) %>%
    select(incoming_class, all_of(vars)) %>%
    pivot_longer(cols = -incoming_class, names_to = "marker", values_to = "raw_value") %>%
    group_by(incoming_class, marker) %>%
    summarise(
      mean_val = mean(raw_value, na.rm = TRUE),
      sd_val   = sd(raw_value, na.rm = TRUE),
      n        = n(),
      se       = sd_val / sqrt(n),
      ci_lower = mean_val - (1.96 * se),
      ci_upper = mean_val + (1.96 * se),
      cohort   = cohort_label,
      .groups = "drop"
    )
}

# 2. Extract Parameters using the objects defined in 16.R and 18.R
# 'lta_input_copy' (Original Raw) and 'lta_ext_copy' (ext Raw) are from 16.R
# 'chronic_scores' (Original incoming_classes) and 'chronic_scores_ext' (ext incoming_classes) are from 12.R/18.R
params_orig <- extract_incoming_class_parameters(lta_input_copy, chronic_scores, chronic_vars_final, "Original")
params_ext  <- extract_incoming_class_parameters(lta_ext_copy, chronic_scores_ext, chronic_vars_final, "ext")

# --- 3. Combine and Audit Overlap (FIXED) ---
comparison_audit <- bind_rows(params_orig, params_ext) %>%
  # IMPORTANT: We drop n, se, and sd_val so pivot_wider can join by 
  # just incoming_class and marker
  select(incoming_class, marker, mean_val, ci_lower, ci_upper, cohort) %>%
  pivot_wider(names_from = cohort, values_from = c(ci_lower, ci_upper, mean_val)) %>%
  # Now that they are in the same row, we can compare them
  mutate(
    overlap_success = if_else(
      ci_lower_ext <= ci_upper_Original & ci_upper_ext >= ci_lower_Original,
      "SUCCESS", "FAIL"
    )
  )

# 4. Visualization
p_param_overlap <- ggplot(bind_rows(params_orig, params_ext), aes(x = cohort, y = mean_val, color = cohort)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.3, size = 1) +
  facet_wrap(marker ~ incoming_class, scales = "free_y") +
  theme_minimal() +
  scale_color_manual(values = c("Original" = "#2c3e50", "ext" = "#e74c3c")) +
  labs(title = "Chronic-State Parameter Overlap: Original vs. ext",
       subtitle = "Validation Success = Overlapping 95% Confidence Intervals",
       y = "Raw Physiological Value", x = "")

print(p_param_overlap) # <--- Fixed: Explicit print for the script execution

# --- 5. Summary Reporting (FIXED) ---
success_rate <- mean(comparison_audit$overlap_success == "SUCCESS", na.rm = TRUE) * 100
message(paste0("Parameter Overlap Success Rate: ", round(success_rate, 1), "%"))

if(!is.na(success_rate) && success_rate >= 80) {
  message("--- EXTERNAL VALIDATION STEP 3: PASSED ---")
} else {
  message("--- EXTERNAL VALIDATION STEP 3: MARGINAL/FAIL ---")
  # Filter for FAIL to see which marker/class combo is the outlier
  print(comparison_audit %>% 
          filter(overlap_success == "FAIL") %>% 
          select(incoming_class, marker, mean_val_Original, mean_val_ext))
}






message("\n--- STARTING SECTION 5.4: EXTERNAL FINGERPRINT & DRIFT VALIDATION ---")

# PART 1: Static physiological profiles

# 1. Bridge External Raw Data with External Acute Class assignments
static_analysis_ext <- lta_ext_copy %>%
  select(all_of(c(ENCOUNTER_ID_VAR, chronic_vars_final, acute_vars_final))) %>%
  inner_join(
    final_lta_analysis_ext %>% 
      select(!!sym(ENCOUNTER_ID_VAR), node_label, outgoing_class), 
    by = ENCOUNTER_ID_VAR
  ) %>%
  filter(node_label %in% c("lysis_confirm", "intraperitoneal_confirm")) %>%
  filter(!is.na(outgoing_class))

# 2. Calculate Medians/IQR
static_summary_ext <- static_analysis_ext %>%
  group_by(outgoing_class) %>%
  summarise(across(all_of(c(chronic_vars_final, acute_vars_final)), list(
    med = ~median(.x, na.rm = TRUE),
    q1  = ~quantile(.x, 0.25, na.rm = TRUE),
    q3  = ~quantile(.x, 0.75, na.rm = TRUE)
  )), .groups = "drop")

# 3. Visualization: The External "Environmental Fingerprint"
top_markers <- chronic_vars_final

plot_ext_fingerprint <- static_summary_ext %>%
  pivot_longer(-outgoing_class, names_to = "metric", values_to = "value") %>%
  separate(metric, into = c("variable", "stat"), sep = "_(?=[^_]+$)") %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  filter(variable %in% top_markers)

p_ext_fingerprint <- ggplot(plot_ext_fingerprint, aes(x = outgoing_class, y = med, color = outgoing_class)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = q1, ymax = q3), width = 0.2) +
  facet_wrap(~variable, scales = "free_y") + 
  theme_minimal() +
  scale_color_manual(values = c("1" = "#d73027", "2" = "#4575b4")) + 
  labs(title = "EXTERNAL VALIDATION: Clinical Fingerprint",
       subtitle = "INSPIRE Cohort: Medians and IQR")

print(p_ext_fingerprint)


# PART 2: External cohort drift analysis
drift_data_ext <- final_lta_analysis_ext %>%
  group_by(!!sym(PATIENT_ID_VAR)) %>%
  arrange(node_sequence) %>%
  mutate(months_since_index = cumsum(days_between_nodes) / 30.44) %>%
  # Join raw clinical values using the Global Encounter ID
  left_join(lta_ext_copy %>% 
              select(!!sym(ENCOUNTER_ID_VAR), all_of(chronic_vars_final)), 
            by = ENCOUNTER_ID_VAR) %>%
  filter(!is.na(directionality)) %>%
  ungroup()

# Defined covariates
ext_cov_list <- c(
  "cov_age", "cov_bmi", "cov_asa", "cov_systemic_trauma", "cov_diabetes",
  "cov_ckd", "cov_gender", "cov_last_approach", "cov_local_trauma", "cov_last_urgency")

# 1. Identify Valid Transitions (Matches Original Cohort's N >= 5 rule)
valid_transitions_ext <- drift_data_ext %>%
  count(directionality) %>%
  filter(!is.na(directionality), n >= 5) %>%
  pull(directionality)

# 2. Automated LME Engine
drift_results_ext <- map_df(chronic_vars_final, ~{
  message("Modeling Drift for: ", .x)
  
  # Filter for valid transitions and drop NAs for this specific marker
  analysis_df <- drift_data_ext %>%
    filter(directionality %in% valid_transitions_ext) %>%
    drop_na(all_of(c(.x, "months_since_index", ext_cov_list)))
  
  if(nrow(analysis_df) < 30) return(NULL)
  
  # Formula logic: 
  formula_str <- paste0(
    .x, " ~ months_since_index * directionality + ", 
    paste(ext_cov_list, collapse = " + "), 
    " + (1 + months_since_index | ", PATIENT_ID_VAR, ")"
  )
  
  mod <- tryCatch({
    lmer(as.formula(formula_str), data = analysis_df)
  }, error = function(e) {
    message("Model too complex for ", .x, "; falling back to Random Intercept.")
    simple_f <- gsub(paste0("\\(1 \\+ months_since_index \\| ", PATIENT_ID_VAR, "\\)"), 
                     paste0("(1 | ", PATIENT_ID_VAR, ")"), formula_str)
    lmer(as.formula(simple_f), data = analysis_df)
  })
  
  if(!is.null(mod)) {
    as.data.frame(summary(mod)$coefficients) %>%
      rownames_to_column("term") %>%
      filter(str_detect(term, "months_since_index:directionality")) %>%
      mutate(variable = .x)
  }
})


# PART 3: Report on directional consistency across cohorts
message("\n--- CROSS-COHORT DIRECTIONAL CONSISTENCY ---")

if(exists("drift_results_full")) {
  validation_compare <- drift_results_ext %>%
    select(variable, term, Estimate_Ext = Estimate, p_ext = `Pr(>|t|)`) %>%
    inner_join(drift_results_full %>% 
                 select(variable, term, Estimate_Orig = Estimate), 
               by = c("variable", "term")) %>%
    mutate(
      term = str_replace(term, "months_since_index:directionality", "Velocity: "),
      consistent = sign(Estimate_Ext) == sign(Estimate_Orig)
    )
  
  # Print as a data frame
  print(as.data.frame(validation_compare))
  
  success_count <- sum(validation_compare$consistent & validation_compare$p_ext < 0.05, na.rm = TRUE)
  message(paste0("\nDirectional Successes: ", success_count, " out of ", nrow(validation_compare)))
}


# PART 4: External cohort physiological drift visualization
message("\n--- GENERATING FOREST PLOT ---")

if(nrow(drift_results_ext) > 0) {
  drift_plot_ext <- drift_results_ext %>%
    mutate(
      # Extract the transition (e.g., 1 -> 2) from the term
      group = str_remove(term, "months_since_index:directionality"),
      # Clean up names for Y-axis
      variable_label = toupper(str_remove(variable, "pre_"))
    )
  
  p_velocity_ext <- ggplot(drift_plot_ext, aes(x = Estimate, y = variable_label, color = group)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_point(size = 3, position = position_dodge(width = 0.7)) +
    geom_errorbarh(aes(xmin = Estimate - `Std. Error`, xmax = Estimate + `Std. Error`), 
                   height = 0.3, position = position_dodge(width = 0.7)) +
    theme_minimal() +
    labs(title = "EXTERNAL VALIDATION: Physiological Drift Velocity",
         subtitle = "INSPIRE Cohort: Rate of change per month (Adjusted)",
         x = "Velocity (Units/Month)", 
         y = "Clinical Marker", 
         color = "Transition Path") +
    theme(axis.text.y = element_text(size = 10), legend.position = "bottom")
  
  print(p_velocity_ext)
}

# PART 5: Transition probability matrix
message("\n--- GENERATING EXTERNAL TRANSITION PROBABILITY MATRIX ---")

transition_matrix_ext <- drift_data_ext %>%
  group_by(!!sym(PATIENT_ID_VAR)) %>% 
  mutate(lead_class = lead(incoming_class)) %>%
  ungroup() %>%
  filter(!is.na(lead_class)) %>%
  count(incoming_class, lead_class) %>%
  group_by(incoming_class) %>%
  mutate(prob = round(n / sum(n), 3)) %>%
  select(-n) %>%
  pivot_wider(names_from = lead_class, values_from = prob, values_fill = 0, names_prefix = "To_")

print(transition_matrix_ext)

# PART 6: Validation summary
message("\n--- FINAL EXTERNAL VALIDATION SUMMARY ---")
message("1. Static Fingerprint: ", nrow(static_summary_ext), " classes profiled.")
message("2. Drift Analysis: ", length(unique(drift_results_ext$variable)), " markers modeled.")
if(exists("success_count")) {
  message("3. Cross-Cohort Consistency: ", success_count, " markers significant and directionally aligned.")
}
