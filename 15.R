source("14.R")
library(ggplot2)
library(lme4)
library(lmerTest)

# --- SECTION 4.A: UNIFIED STATIC ANALYSIS (ACUTE-RESPONSE LENS) ---
message("\n--- GENERATING STATIC PROFILES: LYSIS & INTRAPERITONEAL (ACUTE LENS) ---")

# 1. Identify the Target Encounters & Features
perioperative_vars <- c(chronic_vars_final, acute_vars_final)
# 2. Bridge Raw Data with the Acute Class assignments
static_analysis_unified <- lta_input_copy %>%
  select(all_of(c(ENCOUNTER_ID_VAR, perioperative_vars))) %>%
  inner_join(
    final_lta_analysis %>% 
      select(!!sym(ENCOUNTER_ID_VAR), node_label, outgoing_class, outgoing_prob), 
    by = ENCOUNTER_ID_VAR
  ) %>%
  filter(node_label %in% c("lysis_confirm", "intraperitoneal_confirm")) %>%
  filter(!is.na(outgoing_class))
# 3. Calculate Medians and IQR by the Acute Phenotype
static_summary_results <- static_analysis_unified %>%
  group_by(outgoing_class) %>%
  summarise(across(all_of(perioperative_vars), list(
    med = ~median(.x, na.rm = TRUE),
    q1  = ~quantile(.x, 0.25, na.rm = TRUE),
    q3  = ~quantile(.x, 0.75, na.rm = TRUE)
  )), .groups = "drop")
# 4. Format for Publication (Table 1 Style)
static_table_acute <- static_summary_results %>%
  pivot_longer(-outgoing_class, names_to = "metric", values_to = "value") %>%
  separate(metric, into = c("variable", "stat"), sep = "_(?=[^_]+$)") %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  mutate(summary_stat = paste0(round(med, 2), " [", round(q1, 2), "-", round(q3, 2), "]")) %>%
  select(outgoing_class, variable, summary_stat)
# 5. Summary Report
message("Unified Static Analysis complete for N = ", nrow(static_analysis_unified), " encounters.")
print(static_table_acute)
# --- SECTION 4.A.2: VISUALIZING THE ENVIRONMENTAL FINGERPRINT (WITH FIX) ---
# 1. NEW STEP: Build 'plot_final' from 'static_summary_results'
top_variables <- c("pre_ph", "pre_lactate", "pre_crp", "pre_wbc", 
                   "post_ph_slope", "post_lactate_slope", "feat_vis_score")
plot_final <- static_summary_results %>%
  pivot_longer(-outgoing_class, names_to = "metric", values_to = "value") %>%
  separate(metric, into = c("variable", "stat"), sep = "_(?=[^_]+$)") %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  filter(variable %in% top_variables) %>%
  mutate(
    variable_label = str_remove(variable, "pre_|feat_"),
    variable_label = str_replace_all(variable_label, "_", " "),
    variable_label = toupper(variable_label)
  )
# 2. Generate the plot object
p_fingerprint <- ggplot(plot_final, aes(x = outgoing_class, y = med, color = outgoing_class)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = q1, ymax = q3), width = 0.2, linewidth = 1) +
  facet_wrap(~variable_label, scales = "free_y") + 
  theme_minimal() +
  scale_color_manual(values = c("1" = "#d73027", "2" = "#4575b4")) + 
  labs(title = "Clinical Fingerprint of Surgical Phenotypes",
       subtitle = "Medians and Interquartile Ranges (Raw Units)",
       x = "Assigned Acute-Response Phenotype",
       y = "Clinical Value",
       color = "Phenotype") +
  theme(strip.text = element_text(face = "bold", size = 10),
        legend.position = "bottom",
        panel.spacing = unit(1.5, "lines"))
# 3. FORCE THE PRINT
print(p_fingerprint)
# 4. CONSOLE VERIFICATION
message("\n--- QUICK DELTA CHECK: SICK (1) vs STABLE (2) ---")
print(plot_final %>%
        select(outgoing_class, variable_label, med) %>%
        pivot_wider(names_from = outgoing_class, values_from = med, names_prefix = "Class_") %>%
        mutate(delta = Class_1 - Class_2,
               percent_diff = (delta / Class_2) * 100))


# --- 4.B. LONGITUDINAL INSIGHT: DRIFT & TRANSITIONS ---

# 1. PREPARE LONGITUDINAL DATA
# We ensure we use the global ID variables and the specific Chronic-State classes
longitudinal_drift_data <- final_lta_analysis %>%
  group_by(!!sym(PATIENT_ID_VAR)) %>%
  arrange(node_sequence) %>%
  mutate(
    # Time since the very first surgery in the record (Index)
    months_since_index = cumsum(days_between_nodes) / 30.44
  ) %>%
  # Join raw clinical values using the Global Encounter ID
  left_join(lta_input_copy %>% 
              select(!!sym(ENCOUNTER_ID_VAR), all_of(chronic_vars_final)), 
            by = ENCOUNTER_ID_VAR) %>%
  filter(!is.na(directionality)) %>%
  ungroup()

# 2. DEFINE EXPLICIT COVARIATES (Appendix 4.A)
# We list them manually to avoid 'grepl' confusion with 'covariance' objects
appendix_4a_covs <- c("cov_age", "cov_bmi", "cov_asa", "cov_systemic_trauma", "cov_diabetes", "cov_ckd",
                      "cov_gender", "cov_last_approach", "cov_local_trauma", "cov_last_urgency")

# 3. IDENTIFY VALID TRANSITIONS (N >= 5)
valid_transitions <- longitudinal_drift_data %>%
  count(directionality) %>%
  filter(!is.na(directionality), n >= 5) %>%
  pull(directionality)

# 4. AUTOMATED LME ENGINE (Standardized IDs & Explicit Covs)
run_refined_drift <- function(target_var, df, transitions, covs) {
  message("Modeling Drift for: ", target_var)
  
  # Ensure target, time, and ALL covariates are present for these rows
  analysis_df <- df %>%
    filter(directionality %in% transitions) %>%
    drop_na(all_of(c(target_var, "months_since_index", covs)))
  
  if(nrow(analysis_df) < 30) {
    message("Insufficient data for ", target_var)
    return(NULL)
  }
  
  # Formula: Time * Group + Appendix 4.A Fixed Effects + Random Slope/Intercept
  # We use PATIENT_ID_VAR for the random effect grouping
  formula_str <- paste0(
    target_var, " ~ months_since_index * directionality + ", 
    paste(covs, collapse = " + "), 
    " + (1 + months_since_index | ", PATIENT_ID_VAR, ")"
  )
  
  model <- tryCatch({
    lmer(as.formula(formula_str), data = analysis_df)
  }, error = function(e) {
    message("Model too complex for ", target_var, "; falling back to Random Intercept.")
    simple_f <- gsub(paste0("\\(1 \\+ months_since_index \\| ", PATIENT_ID_VAR, "\\)"), 
                     paste0("(1 | ", PATIENT_ID_VAR, ")"), formula_str)
    lmer(as.formula(simple_f), data = analysis_df)
  })
  return(model)
}

# --- 5. EXECUTE FULL PILLAR ANALYSIS (All 19 Chronic Vars) ---
# We use chronic_vars_final directly as defined in your earlier scripts
pillar_vars_full <- chronic_vars_final 

drift_results_full <- map_df(pillar_vars_full, ~{
  mod <- run_refined_drift(.x, longitudinal_drift_data, valid_transitions, appendix_4a_covs)
  if(!is.null(mod)) {
    as.data.frame(summary(mod)$coefficients) %>%
      rownames_to_column("term") %>%
      filter(str_detect(term, "months_since_index:directionality")) %>%
      mutate(variable = .x)
  }
})

# --- 6. VISUALIZE FULL VELOCITY (Expanded Forest Plot) ---
drift_plot_full <- drift_results_full %>%
  mutate(
    group = str_remove(term, "months_since_index:directionality"),
    # Clean up names: remove 'pre_' and make uppercase for the Y-axis
    variable_label = toupper(str_remove(variable, "pre_"))
  )

p_velocity_full <- ggplot(drift_plot_full, aes(x = Estimate, y = variable_label, color = group)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_point(size = 3, position = position_dodge(width = 0.7)) +
  geom_errorbarh(aes(xmin = Estimate - `Std. Error`, xmax = Estimate + `Std. Error`), 
                 height = 0.3, position = position_dodge(width = 0.7)) +
  theme_minimal() +
  labs(title = "Comprehensive Physiological Drift Velocity",
       subtitle = "Rate of change per month for all Chronic State variables (Adjusted)",
       x = "Velocity (Units/Month)", 
       y = "Clinical Marker", 
       color = "Transition Path") +
  theme(
    axis.text.y = element_text(size = 8), # Smaller text to fit all 19
    legend.position = "bottom"
  )

# PRO TIP: Use ggsave with a larger height to prevent squishing
print(p_velocity_full)

# 7. GENERATE TRANSITION PROBABILITY MATRIX
transition_matrix <- longitudinal_drift_data %>%
  group_by(!!sym(PATIENT_ID_VAR)) %>% # Ensures lead() stays within the same patient
  mutate(lead_class = lead(incoming_class)) %>%
  ungroup() %>%
  filter(!is.na(lead_class)) %>%
  count(incoming_class, lead_class) %>%
  group_by(incoming_class) %>%
  mutate(prob = round(n / sum(n), 3)) %>%
  select(-n) %>%
  pivot_wider(names_from = lead_class, values_from = prob, values_fill = 0, names_prefix = "To_")

message("\n--- TRANSITION PROBABILITY MATRIX ---")
print(transition_matrix)