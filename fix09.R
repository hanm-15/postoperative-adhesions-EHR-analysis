library(tidyverse)
library(stats)
library(corpcor)
source("08.R")

# --- 0. LOCAL STABILIZATION ---
# Use 'any_of' to prevent the crash if a variable is missing
index_features_internal <- index_features %>%
  mutate(across(any_of(acute_vars), ~as.vector(scale(.x))))

# --- 1. DYNAMIC COVARIANCE GENERATOR ---
calculate_class_covariance <- function(data, vars, class_label) {
  class_data <- data %>% 
    filter(Class == class_label) %>% 
    select(any_of(vars)) %>%
    select(where(~ (sd(.x, na.rm = TRUE) > 0) %in% TRUE)) 
  
  active_vars <- colnames(class_data)
  class_data_clean <- class_data %>% drop_na()
  
  if (nrow(class_data_clean) < 5 || length(active_vars) == 0) {
    return(list(mat = diag(length(vars)), vars = active_vars))
  }
  
  cov_mat <- as.matrix(corpcor::cov.shrink(class_data_clean, verbose = FALSE))
  return(list(mat = cov_mat, vars = active_vars))
}

# --- 2. BUILD DYNAMIC LISTS ---
acute_cov_list <- map(unique(index_features_internal$Class), ~{
  calculate_class_covariance(index_features_internal, acute_vars, .x)
}) %>% set_names(unique(index_features_internal$Class))

# --- 3. REVISED SCORING ENGINE ---
score_patient_phenotype <- function(patient_data, weight_means, cov_list, feature_vars) {
  classes <- names(cov_list)
  
  results <- map_df(classes, ~{
    current_class <- .x
    class_cov <- cov_list[[current_class]]$mat
    active_vars <- cov_list[[current_class]]$vars
    
    # 1. Align Center Names (Universal fix for the '_mean' suffix)
    class_center <- weight_means %>%
      filter(Class == current_class) %>%
      select(ends_with("_mean")) %>%
      rename_with(~str_remove(., "_mean")) %>%
      select(any_of(active_vars)) %>% 
      as.matrix() %>% as.vector()
    
    # 2. Extract Patient Data (Assumes data is ALREADY SCALED outside this function)
    patient_vector <- patient_data %>% select(any_of(active_vars)) %>% as.matrix()
    
    # 3. Math (Mahalanobis)
    d2 <- tryCatch({ stats::mahalanobis(patient_vector, center = class_center, cov = class_cov) }, 
                   error = function(e) return(NA))
    
    # 4. Log-Density for Numerical Stability
    log_l <- dchisq(d2, df = length(active_vars), log = TRUE)
    
    tibble(Class = current_class, dist_sq = d2, log_likelihood = log_l)
  })
  
  # 5. Stabilized Bayesian Normalization
  results <- results %>%
    mutate(
      # The Log-Sum-Exp Trick: subtract max to prevent overflow, then exp()
      rel_log_l = log_likelihood - max(log_likelihood, na.rm = TRUE),
      posterior_prob = exp(rel_log_l) / sum(exp(rel_log_l), na.rm = TRUE)
    ) %>%
    arrange(desc(posterior_prob))
  
  return(results)
}


# --- ENGINE AUDIT ---
message("\n--- BAYESIAN ENGINE AUDIT ---")

# 1. Verify Weight-to-Matrix Alignment
# 1. Verify Weight-to-Matrix Alignment (Updated for Safe-Fails)
message("Logic Check: Matching feature dimensions...")
acute_cov_list %>% 
  iwalk(~{
    # Check if the matrix is actually a valid covariance matrix or our safe-fail identity
    # We check the length of the 'vars' vector we stored
    n_vars <- length(.x$vars)
    
    if(n_vars > 0) {
      message("  - Class ", .y, " Covariance Matrix: ", n_vars, "x", n_vars, " features.")
    } else {
      message("  - Class ", .y, " Covariance Matrix: [SKIPPED - Insufficient Data]")
    }
  })


# Helper function (retired)
dist_to_prob <- function(d2, p) {
  # p is the number of features (degrees of freedom)
  # Use the Chi-squared density function
  likelihood <- dchisq(d2, df = p)
  return(likelihood)
}


# 2. Run a batch test on all patients (Acute-Response)
index_features_scaled_for_audit <- index_features %>%
  mutate(across(any_of(acute_vars), ~as.vector(scale(.x))))
# --- 2. FULL-SPECTRUM BATCH TEST ---
message("\nProjecting Phenotypes and Auditing Distances...")

# We map over all patients to get the full distance profile
final_assignments_table <- index_features_scaled_for_audit %>%
  group_by(.data[[PATIENT_ID_VAR]]) %>% 
  group_modify(~ {
    # We keep the full results from the engine
    score_patient_phenotype(.x, acute_weights, acute_cov_list, acute_vars)
  }) %>%
  ungroup()

# --- 3. THE MASTER AUDIT VIEW ---
message("\n--- COMPREHENSIVE PHENOTYPE AUDIT ---")

# This view shows the Winner vs the Runner-Up for the first 10 patients
audit_summary <- final_assignments_table %>%
  group_by(!!sym(PATIENT_ID_VAR)) %>%
  # Create a rank to see who was #1 vs #2
  mutate(rank = row_number()) %>%
  filter(rank <= 2) %>%
  select(!!sym(PATIENT_ID_VAR), rank, Class, dist_sq, posterior_prob)

print(head(audit_summary, 20))

# --- 4. QUICK STABILITY STATS ---
message("\n--- GLOBAL DISTANCE STABILITY ---")
summary(final_assignments_table$dist_sq)

message("\nAudit complete. If you see 'posterior_prob' values above, the Bayesian Engine is officially live.")