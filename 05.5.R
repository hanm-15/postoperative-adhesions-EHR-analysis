library(tidyverse)
library(tidyLPA)
library(ggrepel)
library(clinfun)
source("LPA_test.R")

#consider dropping postop ileus as a feature

#check to make sure the dx_df logic aligns with that of the other files

# --- 1. Dynamic Best-K Selection ---
best_k <- stats %>%
  filter(BIC == min(BIC, na.rm = TRUE)) %>%
  slice(1) %>% 
  pull(Classes)

message("Lysis Validation: Proceeding with k=", best_k)

# --- 2. Extract Patient-Class Assignments (INDEX-FILE STYLE FIX) ---

# We recreate the bridge from 'raw_surgery_table' to match the 'lca_matrix' row order
id_bridge <- raw_surgery_table %>%
  filter(lysis_confirm == TRUE) %>%
  mutate(id = row_number()) %>%
  # Standardize type immediately to prevent "character vs numeric" join errors
  mutate(!!sym(PATIENT_ID_VAR) := as.character(!!sym(PATIENT_ID_VAR))) %>%
  select(id, !!sym(PATIENT_ID_VAR))

final_assignments <- get_data(model_results) %>%
  filter(classes_number == best_k) %>%
  mutate(id = as.numeric(id)) %>% 
  left_join(id_bridge, by = "id") %>%
  select(Class, !!sym(PATIENT_ID_VAR)) %>%
  distinct() %>%
  # Standardize for the final merge
  mutate(!!sym(PATIENT_ID_VAR) := as.character(!!sym(PATIENT_ID_VAR))) %>%
  mutate(Class = factor(Class, 
                        levels = sort(unique(as.numeric(as.character(Class)))), 
                        ordered = TRUE))

# --- 3. DEFINE TEMPORAL PILLARS (REPAIRED) ---
lysis_anchors <- raw_surgery_table %>%
  filter(lysis_confirm == TRUE) %>%
  select(
    !!sym(PATIENT_ID_VAR), 
    !!sym(ENCOUNTER_ID_VAR), 
    t1 = !!sym(SURG_END_VAR),
    discharge_time = !!sym(DISCHARGE_VAR),
    icu_admit_time = !!sym(ICU_IN_VAR),   # <--- ADD THIS
    icu_exit_time  = !!sym(ICU_OUT_VAR)    # <--- ADD THIS
  ) %>%
  mutate(t1_24 = t1 + hours(24), t1_48 = t1 + hours(48))

# --- 4. THE MATH ENGINES ---
calc_slope <- function(df, time_col, value_col) {
  if(!value_col %in% names(df) || all(is.na(df[[value_col]]))) return(0)
  df_clean <- df %>% filter(!is.na(!!sym(value_col)))
  if(nrow(df_clean) < 2) return(0)
  model <- lm(as.numeric(!!sym(value_col)) ~ as.numeric(!!sym(time_col)), data = df_clean)
  return(as.numeric(coef(model)[2] * 3600)) 
}

# --- 5. THE REFINED LYSIS VALIDATION FACTORY ---
get_lys_validation <- function(anchors, labs, vitals, vaso) {
  l_vars <- c(LACTATE_VAR, BUN_VAR, CREAT_VAR, BE_VAR, ALBUMIN_VAR, WBC_VAR, NEUTRO_VAR, LYMPHO_VAR, PLATE_VAR)
  v_vars <- c(MAP_VAR, HR_VAR)
  
  lab_s <- labs %>% add_missing_cols(l_vars)
  vit_s <- vitals %>% add_missing_cols(v_vars)
  # Ensure vaso table has the columns we check for weaning
  vas_s <- vaso %>% add_missing_cols(c(EPI_VAR, NORE_VAR, VASO_VAR, DOPE_VAR, DOBU_VAR, MIL_VAR))
  
  anchors %>%
    group_by(!!sym(PATIENT_ID_VAR)) %>%
    group_modify(~ {
      p_id <- .y[[1]]
      t1 <- .x$t1; t24 <- .x$t1_24; t48 <- .x$t1_48; tdis <- .x$discharge_time
      
      p_labs <- lab_s %>% filter(!!sym(PATIENT_ID_VAR) == p_id)
      p_vits <- vit_s %>% filter(!!sym(PATIENT_ID_VAR) == p_id)
      p_vaso <- vas_s %>% filter(!!sym(PATIENT_ID_VAR) == p_id)
      
      l_48 <- p_labs %>% filter(!!sym(CHART_TIME_VAR) >= t1, !!sym(CHART_TIME_VAR) <= t48)
      v_48 <- p_vits %>% filter(!!sym(CHART_TIME_VAR) >= t1, !!sym(CHART_TIME_VAR) <= t48)
      v_24 <- v_48 %>% filter(!!sym(CHART_TIME_VAR) <= t24)
      
      tibble(
        # 1. Recovery Slopes
        lys_lactate_slope = calc_slope(l_48, CHART_TIME_VAR, LACTATE_VAR),
        lys_be_slope      = calc_slope(l_48, CHART_TIME_VAR, BE_VAR),
        lys_alb_slope     = calc_slope(l_48, CHART_TIME_VAR, ALBUMIN_VAR),
        lys_wbc_slope     = calc_slope(l_48, CHART_TIME_VAR, WBC_VAR),
        lys_plt_slope     = calc_slope(l_48, CHART_TIME_VAR, PLATE_VAR),
        lys_hr_slope      = calc_slope(v_48, CHART_TIME_VAR, HR_VAR),
        lys_nlr_slope     = {
          tmp <- l_48 %>% mutate(nlr = !!sym(NEUTRO_VAR) / pmax(0.1, !!sym(LYMPHO_VAR)))
          calc_slope(tmp, CHART_TIME_VAR, "nlr")
        },
        
        # 2. Hemodynamic Volatility
        lys_map_cv = (sd(v_24[[MAP_VAR]], na.rm=T) / mean(v_24[[MAP_VAR]], na.rm=T)) * 100,
        
        # 3. Refined BUN/Creatinine Delta (First Pre-Op vs Post-Op Median)
        lys_bun_creat_delta = {
          pre_ratio <- p_labs %>% 
            filter(!!sym(CHART_TIME_VAR) < t1) %>% 
            arrange(!!sym(CHART_TIME_VAR)) %>% # Arranging ascending to get the FIRST value
            summarize(r = first(!!sym(BUN_VAR))/pmax(0.1, first(!!sym(CREAT_VAR)))) %>% pull(r)
          
          post_ratio <- l_48 %>% 
            summarize(r = median(!!sym(BUN_VAR)/pmax(0.1, !!sym(CREAT_VAR)), na.rm=T)) %>% pull(r)
          
          replace_na(post_ratio - pre_ratio, 0)
        },
        
        # 4. Refined Vaso Weaning
        lys_vaso_wean_hrs = {
          # Filter for active pressors
          pressor_times <- p_vaso %>%
            filter(!!sym(CHART_TIME_VAR) >= t1, !!sym(CHART_TIME_VAR) <= t48) %>%
            filter(!!sym(EPI_VAR) > 0 | !!sym(NORE_VAR) > 0 | !!sym(VASO_VAR) > 0 | 
                     !!sym(DOPE_VAR) > 0 | !!sym(DOBU_VAR) > 0 | !!sym(MIL_VAR) > 0) %>%
            pull(!!sym(CHART_TIME_VAR))
          
          if(length(pressor_times) > 0) {
            val <- as.numeric(difftime(max(pressor_times), t1, units = "hours"))
            pmax(0, val)
          } else {
            0
          }
        },
        
        # 5. Outcomes (Back to your preferred !!sym syntax)
        # Inside your tibble() in Section 5, update the ICU LOS line:
        lys_icu_los = replace_na(pmax(0, as.numeric(difftime(.x$icu_exit_time, .x$icu_admit_time, units = "days"))), 0),
        lys_total_los = as.numeric(difftime(tdis, t1, units = "days")),
        lys_ileus_bin = if_else(any(grepl("^K56\\.0|^K56\\.7|^K91\\.3", .x[[DX_VAR]])), 1, 0)
      )
    }) %>% ungroup()
}

# --- PASTE THE FINAL REPAIR BLOCK HERE ---
if(ncol(lab_results) > 0)    { lab_results    <- lab_results    %>% rename_with(~PATIENT_ID_VAR, contains("pat_id")) }
if(ncol(vitals_results) > 0) { vitals_results <- vitals_results %>% rename_with(~PATIENT_ID_VAR, contains("pat_id")) }

if (nrow(vaso_results) == 0 || ncol(vaso_results) == 0) {
  vaso_results <- tibble(
    !!sym(PATIENT_ID_VAR) := character(), 
    !!sym(CHART_TIME_VAR) := as.POSIXct(character())
  ) %>%
    add_missing_cols(c(EPI_VAR, NORE_VAR, VASO_VAR, DOPE_VAR, DOBU_VAR, MIL_VAR))
} else {
  vaso_results <- vaso_results %>% rename_with(~PATIENT_ID_VAR, contains("pat_id"))
}

lab_results[[PATIENT_ID_VAR]]    <- as.character(lab_results[[PATIENT_ID_VAR]])
vitals_results[[PATIENT_ID_VAR]] <- as.character(vitals_results[[PATIENT_ID_VAR]])
# ------------------------------------------

# --- 6. EXECUTION (The missing bridge) ---
validation_metrics <- get_lys_validation(lysis_anchors, lab_results, vitals_results, vaso_results)
final_validation_set <- final_assignments %>% inner_join(validation_metrics, by = PATIENT_ID_VAR)

# --- RESEARCH PLAN: TARGET LYSIS FEATURES ---
# This matches your specific research goals for the Reoperated subset
lysis_features_to_test <- c(
  "lys_lactate_slope", 
  "lys_be_slope", 
  "lys_alb_slope", 
  "lys_wbc_slope", 
  "lys_plt_slope", 
  "lys_hr_slope",
  "lys_nlr_slope", 
  "lys_map_cv", 
  "lys_bun_creat_delta",
  "lys_icu_los", 
  "lys_ileus_bin",
  "lys_vaso_wean_hrs",
  "lys_total_los"
)

# --- 7. THE DYNAMIC TREND ENGINE (Fortified & Manual) ---

# 1. Check which features actually exist in the final validation matrix
available_features <- intersect(lysis_features_to_test, names(final_validation_set))

message("Lysis Trend Test starting.")
message("Requested: ", length(lysis_features_to_test))
message("Found in data: ", length(available_features))

# 2. Run the Jonckheere-Terpstra Loop
trend_results <- map_dfr(available_features, function(feat) {
  
  # Remove NAs and ensure Class is available/ordered
  clean_data <- final_validation_set %>% 
    filter(!is.na(.data[[feat]]), !is.na(Class))
  
  # CRITICAL GUARD: Dummy data often lacks variation (all values the same)
  # We need at least 2 unique values to calculate a trend
  if(nrow(clean_data) < 5 || length(unique(clean_data[[feat]])) < 2) {
    return(NULL)
  }
  
  tryCatch({
    # Perform the JT Test
    jt_stat <- jonckheere.test(clean_data[[feat]], as.numeric(clean_data$Class))
    
    # Calculate direction via Spearman Correlation
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

# 3. ONLY adjust and plot if we actually found valid trends
if(!is.null(trend_results) && nrow(trend_results) > 0) {
  
  trend_results <- trend_results %>%
    mutate(
      Adjusted_P = p.adjust(Raw_P, method = "fdr"),
      Significant = if_else(Adjusted_P < 0.05, "YES", "no")
    ) %>%
    arrange(Adjusted_P)
  
  print(trend_results)
  
  # Prepare Plotting Data (Threshold p < 0.2 for "leaning" trends)
  plot_data <- trend_results %>% 
    mutate(log_p = -log10(Adjusted_P)) %>%
    filter(Adjusted_P < 0.2)
  
  if(nrow(plot_data) > 0) {
    # [Insert the ggplot lollipop code we discussed earlier]
    message("Generating Trend Plot...")
  } else {
    message("Trends found, but none reached p < 0.2 for plotting.")
  }
  
} else {
  message("!!! NO VALID TRENDS FOUND !!!")
  message("This is expected if dummy data is static or features are missing.")
}