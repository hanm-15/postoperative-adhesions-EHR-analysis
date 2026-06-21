library(tidyverse)
library(tidyLPA)
library(ggrepel)
library(clinfun)
source("LPA_test.R")
source("fix07.R")


#change p threshold to 0.2 later, line 191
#add better visualization (e.g., all trajectories not just variance)

# 4. CONTINUE WITH VALIDATION
lca_final[[PATIENT_ID_VAR]] <- as.character(lca_final[[PATIENT_ID_VAR]])

#Phenotypic Validation goes here

# 1. Pull the best K
best_k <- stats %>%
  # --- THE FORCE: Only consider models with 2 or more classes ---
  filter(Classes >= 2) %>% 
  # Now find the best one among those survivors
  filter(BIC == min(BIC, na.rm = TRUE)) %>%
  pull(Classes)

# Safety check: if for some reason NO models exist >= 2, default to 2
if(length(best_k) == 0) { best_k <- 2 }

message("Proceeding with k=", best_k, " as the optimal phenotypic fit.")

# --- 2. Extract Patient-Class Assignments (INCLUDING CLASS 0) ---

# A. Bridge for the LPA results (Reoperated Patients)
id_bridge <- lca_final %>% 
  mutate(id = row_number()) %>% 
  select(id, !!sym(PATIENT_ID_VAR))

reop_assignments <- get_data(model_results) %>%
  filter(classes_number == best_k) %>%
  mutate(id = as.numeric(id)) %>% 
  left_join(id_bridge, by = "id") %>%
  select(Class, !!sym(PATIENT_ID_VAR)) %>%
  # Ensure Class is a character so it can merge with our new Class 0
  mutate(Class = as.character(Class))

# --- 2. Extract Patient-Class Assignments ---

# B. Extract Asymptomatic AND Non-Surgical
extra_assignments <- analysis_data %>%
  filter(patient_group %in% c("asymptomatic", "non_surgical")) %>% 
  select(!!sym(PATIENT_ID_VAR), patient_group) %>%
  distinct() %>%
  # We give Non-Surgical a high number (e.g., 99) so it doesn't interfere 
  # with the 0, 1, 2 sequence, but stays numeric.
  mutate(Class = if_else(patient_group == "asymptomatic", "0", "99")) %>%
  select(-patient_group) 

# C. Combine Assignments
final_assignments <- bind_rows(reop_assignments, extra_assignments) %>%
  distinct() %>%
  mutate(!!sym(PATIENT_ID_VAR) := as.character(!!sym(PATIENT_ID_VAR)))

# --- 3. Pull the Index Surgeries ---
index_validation_data <- analysis_data %>%
  mutate(across(all_of(c(PATIENT_ID_VAR, ENCOUNTER_ID_VAR)), as.character)) %>%
  filter(is_index_surgery == TRUE) %>% 
  left_join(final_assignments, by = PATIENT_ID_VAR) %>%
  mutate(
    t0 = !!sym(SURG_START_VAR),
    t1 = !!sym(SURG_END_VAR)
  ) %>%
  filter(!is.na(Class)) %>%
  # FILTER OUT the Non-Surgical (99) right here before the Trend Test 
  # so they don't skew the ordinal logic
  filter(Class != "99") %>% 
  mutate(Class = factor(Class, 
                        levels = sort(as.numeric(unique(Class))), 
                        ordered = TRUE))

# --- 4. DATA MAPPING (Updated for Stateless 07.R) ---
# We now pass 'index_validation_data' as the reference for t0/t1
lab_results_labeled     <- label_by_index_surgery(lab_results, CHART_TIME_VAR, index_validation_data)
vitals_results_labeled  <- label_by_index_surgery(vitals_results, CHART_TIME_VAR, index_validation_data)
intraop_vitals_labeled  <- vitals_results_labeled %>% filter(.data$period == "intra_op")
vaso_labeled            <- label_by_index_surgery(vaso_results, CHART_TIME_VAR, index_validation_data)
fluid_labeled           <- label_by_index_surgery(fluid_results, CHART_TIME_VAR, index_validation_data)




# --- 5. THE VALIDATION MATRIX GENERATOR ---
pre_features <- get_baseline_features(
  lab_df    = lab_results_labeled,
  vitals_df = vitals_results_labeled,
  id_col    = PATIENT_ID_VAR,
  time_col  = CHART_TIME_VAR
)

intra_features <- get_intraop_features(
  vitals_df   = intraop_vitals_labeled,
  vaso_df     = vaso_labeled,
  fluid_df    = fluid_labeled,
  metadata_df = index_validation_data,
  id_col      = PATIENT_ID_VAR,
  time_col    = CHART_TIME_VAR
)

post_features <- get_postop_features(
  lab_df      = lab_results_labeled,
  vitals_df   = vitals_results_labeled,
  fluid_df    = fluid_labeled,
  vaso_df     = vaso_labeled,
  metadata_df = index_validation_data,
  id_col      = PATIENT_ID_VAR,
  time_col    = CHART_TIME_VAR
)

index_features <- pre_features %>%
  left_join(intra_features, by = PATIENT_ID_VAR) %>%
  left_join(post_features, by = PATIENT_ID_VAR) %>%
  inner_join(final_assignments, by = PATIENT_ID_VAR)

# --- RESEARCH PLAN: TARGET FEATURES FOR TREND TESTING ---
features_to_test <- c(
  # 1. Preoperative Baseline
  "pre_ph", "pre_be", "pre_bicarb", "pre_lactate", "pre_neutro", 
  "pre_plate", "pre_mono", "pre_lympho", "pre_wbc", "pre_crp", 
  "pre_fib", "pre_mpv", "pre_rdw", "pre_hb", "pre_alb", 
  "pre_bun_creat_ratio", "pre_hr", "pre_map", "pre_spo2",
  
  # 2. Intraoperative
  "intra_map_cv", "intra_map_auc_65", "intra_hr_cv", "intra_max_vis", 
  "intra_spo2_auc_90", "intra_etco2_slope", "intra_fio2_auc_60", 
  "intra_temp_delta", "intra_wa_irrig_vol", "intra_wa_irrig_ret", 
  "intra_wa_net_fluid_bal", "intra_wa_cryst_vol", "intra_wa_transfusion", 
  "intra_proc_duration", "intra_wa_ebl", "intra_twa_svv",
  
  # 3. Postoperative
  "post_ph_slope", "post_ph_median", "post_be_slope", "post_lactate_slope", 
  "post_nlr_slope", "post_nlr_median", "post_plr_slope", "post_mlr_slope", 
  "post_wbc_slope", "post_crp_slope", "post_fib_delta", "post_alb_delta", 
  "post_bun_creat_delta", "post_hr_slope", "post_map_cv", 
  "post_vaso_wean_hrs", "post_wa_net_fluid_48", "post_wa_uop_hr", "post_sf_slope"
)

# --- 6. THE DYNAMIC TREND ENGINE (Fortified) ---

# 1. Check which features exist
available_features <- intersect(features_to_test, names(index_features))

message("Trend Test starting. Requested: ", length(features_to_test), 
        " | Found in data: ", length(available_features))

# 2. Run JT Test
trend_results <- map_dfr(available_features, function(feat) {
  
  # Remove NAs and ensure Class is numeric for the test
  clean_data <- index_features %>% 
    filter(!is.na(.data[[feat]])) %>%
    filter(!is.na(Class))
  
  # CRITICAL GUARD: Need at least 2 unique values in the feature AND the Class 
  # to calculate a trend. Dummy data often fails here.
  if(nrow(clean_data) < 5 || 
     length(unique(clean_data[[feat]])) < 2 || 
     length(unique(clean_data$Class)) < 2) return(NULL)
  
  tryCatch({
    jt_stat <- jonckheere.test(clean_data[[feat]], as.numeric(clean_data$Class))
    direction_rho <- cor(clean_data[[feat]], as.numeric(clean_data$Class), 
                         method = "spearman", use = "complete.obs")
    
    data.frame(
      Feature   = feat,
      Raw_P     = jt_stat$p.value,
      Statistic = jt_stat$statistic,
      Trend     = if_else(direction_rho > 0, "Increasing", "Decreasing"),
      n_obs     = nrow(clean_data)
    )
  }, error = function(e) return(NULL))
})

# 3. ONLY adjust and plot if we actually found trends
if(!is.null(trend_results) && nrow(trend_results) > 0) {
  trend_results <- trend_results %>%
    mutate(
      Adjusted_P = p.adjust(Raw_P, method = "fdr"),
      Significant = if_else(Adjusted_P < 0.05, "YES", "no")
    ) %>%
    arrange(Adjusted_P)
  
  print(trend_results)
  
  # Move the Plotting code INSIDE this 'if' block!
  plot_data <- trend_results %>% 
    mutate(log_p = -log10(Adjusted_P)) %>%
    filter(Adjusted_P < 1.0)
  
  if(nrow(plot_data) > 0) {
    # [Insert your ggplot code here]
    # --- DYNAMIC VISUALIZATION ---
    # We use the 'plot_data' created just above in the fortified block
    
    p <- ggplot(plot_data, aes(x = reorder(Feature, log_p), y = log_p)) +
      geom_segment(aes(xend = Feature, yend = 0, linetype = Trend), color = "grey70") +
      geom_point(aes(color = Trend, size = log_p)) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
      coord_flip() +
      labs(
        title = "Multi-Class Phenotypic Trend Analysis",
        subtitle = paste("Jonckheere-Terpstra Test across", length(unique(index_features$Class)), "Ordered Classes"),
        x = "Clinical Feature",
        y = "-log10(Adjusted P)"
      ) +
      theme_minimal()
    
    print(p) # This ensures the plot actually pops up
    print("Plotting results...")
  } else {
    message("Trends found, but none reached the p < 0.2 threshold for plotting.")
  }
} else {
  message("!!! NO VALID TRENDS FOUND !!!")
  message("This usually happens with dummy data because there isn't enough variation across Classes.")
}

# --- FINAL AUTOMATED SUMMARY ---
cat("\n", rep("=", 30), " 06.R SUMMARY ", rep("=", 30), "\n")

# 1. Check Data Presence
cat("Patients in index_features: ", nrow(index_features), "\n")
cat("Features found:            ", length(available_features), "\n")

# 2. Check Class Distribution (The JT Test needs variation here)
if("Class" %in% names(index_features)) {
  cat("Class Distribution:\n")
  print(table(index_features$Class))
} else {
  cat("ERROR: 'Class' column missing from index_features!\n")
}

# 3. Check for NAs in the first few features
cat("\nFirst 5 Features - NA Count:\n")
colSums(is.na(index_features[, head(available_features, 5)])) %>% print()

# 4. Final Status
if(exists("trend_results") && !is.null(trend_results) && nrow(trend_results) > 0) {
  cat("\nSUCCESS: Trend Test produced results. Check the Plots tab.\n")
} else {
  cat("\nWARNING: JT Test returned no results. Check if Class/Features have >1 unique value.\n")
}

cat(rep("=", 75), "\n")