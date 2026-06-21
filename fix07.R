library(tidyverse)
library(lubridate)

#change na -> 0 logic to mice later
#(and guardrails to convert missing features to 0 in calculations)
#use gsub on percentages so that as.numeric works properly
#check if the post_48h column is getting generated properly
#real data, must pivot non-op vitals to wide manually


add_missing_cols <- function(df, cols) {
  for (col in cols) {
    if (!col %in% names(df)) {
      df[[col]] <- NA_real_
    }
  }
  return(df)
}

# --- 1. THE TEMPORAL ENGINE (Fortified & Stateless) ---
label_by_index_surgery <- function(target_data, time_col_name, ref_data) {
  
  # 1. Standardize IDs using tidy eval
  target_data <- target_data %>% 
    mutate(!!sym(ENCOUNTER_ID_VAR) := as.character(!!sym(ENCOUNTER_ID_VAR)))
  
  if (PATIENT_ID_VAR %in% colnames(target_data)) {
    target_data <- target_data %>% 
      mutate(!!sym(PATIENT_ID_VAR) := as.character(!!sym(PATIENT_ID_VAR)))
  }
  
  # 2. Prepare the reference data passed into the function
  ref_clean <- ref_data %>%
    mutate(across(any_of(c(ENCOUNTER_ID_VAR, PATIENT_ID_VAR)), as.character)) %>%
    select(any_of(c(ENCOUNTER_ID_VAR, PATIENT_ID_VAR)), t0, t1) %>%
    distinct()
  
  join_vars <- if (PATIENT_ID_VAR %in% colnames(target_data)) c(ENCOUNTER_ID_VAR, PATIENT_ID_VAR) else ENCOUNTER_ID_VAR
  
  # 3. Join and Label
  target_data %>%
    left_join(ref_clean, by = join_vars) %>%
    mutate(
      # Standardize time for comparison
      tmp_time = as.POSIXct(!!sym(time_col_name)),
      t0 = as.POSIXct(t0),
      t1 = as.POSIXct(t1),
      hrs_post_op = as.numeric(difftime(tmp_time, t1, units = "hours")),
      
      period = case_when(
        is.na(tmp_time) | is.na(t0) ~ "UNKNOWN",
        tmp_time < t0 ~ "pre_op",
        tmp_time >= t0 & tmp_time <= t1 ~ "intra_op",
        hrs_post_op > 0  & hrs_post_op <= 24 ~ "post_24",
        hrs_post_op > 24 & hrs_post_op <= 48 ~ "post_48",
        TRUE ~ "post_late"
      )
    ) %>%
    select(-tmp_time)
}

# --- 2. THE MATH ENGINES ---
calc_slope <- function(df, time_col, value_col) {
  # 1. Shields: Return NA so we know the calculation failed
  if (is.null(df) || nrow(df) < 2) return(NA_real_)
  if(!value_col %in% names(df) || all(is.na(df[[value_col]]))) return(NA_real_)
  
  df_clean <- df[!is.na(df[[value_col]]), ]
  if(nrow(df_clean) < 2) return(NA_real_)
  
  df_clean <- df_clean[order(df_clean[[time_col]]), ]
  
  # 2. Catch model errors (e.g. if all Y values are identical)
  model <- try(lm(as.numeric(df_clean[[value_col]]) ~ as.numeric(df_clean[[time_col]])), silent = TRUE)
  if(inherits(model, "try-error")) return(NA_real_)
  
  return(as.numeric(coef(model)[2] * 3600)) 
}

calc_auc_threshold <- function(df, time_col, value_col, threshold, direction = "below") {
  # 1. Shields
  if (is.null(df) || nrow(df) < 2) return(0)
  if(!value_col %in% names(df) || all(is.na(df[[value_col]]))) return(0)
  
  # 2. Clean and Sort
  df_clean <- df[!is.na(df[[value_col]]), ]
  df_clean <- df_clean[order(df_clean[[time_col]]), ]
  if(nrow(df_clean) < 2) return(0)
  
  # 3. Area Calculation
  # We use simple diff/numeric vectors to avoid Tidyverse mask issues
  times <- as.numeric(df_clean[[time_col]])
  vals  <- as.numeric(df_clean[[value_col]])
  
  dt  <- diff(times) / 60 # convert seconds to minutes
  
  if(direction == "below") {
    gaps <- pmax(0, threshold - head(vals, -1))
  } else {
    gaps <- pmax(0, head(vals, -1) - threshold)
  }
  
  return(sum(gaps * dt, na.rm = TRUE))
}





calc_vis <- function(df) {
  # Helper to ensure we have a number and turn NAs to 0
  safe_val <- function(col_name) {
    if (col_name %in% names(df)) {
      val <- df[[col_name]]
      # If pivot_wider somehow created a list, grab the first number
      if (is.list(val)) val <- as.numeric(vapply(val, function(x) x[1], numeric(1)))
      val <- as.numeric(val)
      return(if_else(is.na(val), 0, val))
    } else {
      return(rep(0, nrow(df)))
    }
  }
  
  # Math using your Global Variables from saving02.R
  df %>% mutate(vis = 
                  safe_val(DOPE_VAR) + 
                  safe_val(DOBU_VAR) + 
                  (100 * safe_val(EPI_VAR)) + 
                  (100 * safe_val(NORE_VAR)) + 
                  (10000 * safe_val(VASO_VAR)) + 
                  (10 * safe_val(MIL_VAR)))
}




# --- 3. THE SUMMARY ENGINES ---
get_period_median <- function(labeled_df, id_col, target_period, value_col) {
  labeled_df %>%
    filter(period == target_period) %>%
    group_by(!!sym(id_col)) %>%
    summarize(median_val = median(!!sym(value_col), na.rm = TRUE), .groups = "drop")
}

calc_clinical_delta <- function(labeled_df, id_col, time_col, value_col) {
  id_sym <- sym(id_col)
  
  labeled_df %>%
    group_by(!!id_sym) %>%
    summarize(
      base_val = {
        # Subset the specific column for pre_op rows
        vals <- .data[[value_col]][.data$period == "pre_op"]
        vals <- vals[!is.na(vals)]
        if(length(vals) > 0) vals[1] else NA_real_
      },
      post_med = {
        # Subset for post_24 rows
        vals <- .data[[value_col]][.data$period == "post_24"]
        if(any(!is.na(vals))) median(vals, na.rm = TRUE) else NA_real_
      },
      .groups = "drop"
    ) %>%
    mutate(delta = replace_na(post_med - base_val, 0)) %>%
    select(!!id_sym, delta)
}


#INTRAOP INDEX FEATURES
ensure_wide <- function(df, id_col, time_col) {
  # If 'variable' exists, the data is in LONG format and needs pivoting
  if ("variable" %in% names(df)) {
    
    # FIX: Group by EVERYTHING except the actual 'value'
    # This ensures 'variable', 'period', 'pat_id', etc. are all preserved
    group_vars <- setdiff(names(df), "value")
    
    return(df %>% 
             dplyr::group_by(across(all_of(group_vars))) %>%
             dplyr::summarize(value = max(value, na.rm = TRUE), .groups = "drop") %>%
             # Now 'variable' still exists and can be pivoted!
             tidyr::pivot_wider(names_from = variable, values_from = value))
  }
  return(df) # If already wide, do nothing
}

get_intraop_features <- function(vitals_df, vaso_df, fluid_df, metadata_df, id_col, time_col, ...) {
  
  vaso_df  <- ensure_wide(vaso_df, id_col, time_col)
  fluid_df <- ensure_wide(fluid_df, id_col, time_col)
  
  # 1. Define all columns we expect
  vitals_vars <- c(MAP_VAR, HR_VAR, SPO2_VAR, FIO2_VAR, ETCO2_VAR, TEMP_VAR, SVV_VAR)
  fluid_vars  <- c(IRRIGATION_VAR, IRRIGATION_IN_VAR, CRYSTALLOID_VAR, 
                   COLLOID_VAR, BLOOD_PROD_VAR, URINE_VAR, EBL_VAR)
  
  # 2. HEMODYNAMICS & SVV (Wide-Compatible)
  vitals_calculated <- vitals_df %>%
    filter(.data$period == "intra_op") %>%
    # No pivot needed! add_missing_cols handles the wide columns
    add_missing_cols(vitals_vars) %>% 
    group_by(!!sym(id_col)) %>%
    summarize(
      # Directly reference the column names from your Global Vars
      intra_map_cv = if(n() > 1 && !all(is.na(!!sym(MAP_VAR)))) (sd(!!sym(MAP_VAR), na.rm = TRUE) / mean(!!sym(MAP_VAR), na.rm = TRUE)) * 100 else 0,
      intra_hr_cv  = if(n() > 1 && !all(is.na(!!sym(HR_VAR)))) (sd(!!sym(HR_VAR), na.rm = TRUE) / mean(!!sym(HR_VAR), na.rm = TRUE)) * 100 else 0,
      
      # Use the specific column names in your Math Engines
      intra_map_auc_65  = if(n() > 1) calc_auc_threshold(pick(everything()), time_col, MAP_VAR, 65) else 0,
      intra_spo2_auc_90 = if(n() > 1) calc_auc_threshold(pick(everything()), time_col, SPO2_VAR, 90) else 0,
      intra_fio2_auc_60 = if(n() > 1) calc_auc_threshold(pick(everything()), time_col, FIO2_VAR, 60, "above") else 0,
      intra_etco2_slope = if(n() > 1) calc_slope(pick(everything()), time_col, ETCO2_VAR) else 0,
      
      intra_temp_delta    = last(!!sym(TEMP_VAR), order_by = !!sym(time_col)) - first(!!sym(TEMP_VAR), order_by = !!sym(time_col)),
      intra_proc_duration = as.numeric(difftime(max(!!sym(time_col)), min(!!sym(time_col)), units = "mins")),
      intra_twa_svv = mean(!!sym(SVV_VAR), na.rm = TRUE),
      .groups = "drop"
    )
  
  # 3. VASOACTIVE LOAD (VIS)
  vaso_calculated <- vaso_df %>%
    filter(.data$period == "intra_op") %>%
    calc_vis() %>% 
    group_by(!!sym(id_col)) %>%
    summarize(intra_max_vis = max(vis, na.rm = TRUE), .groups = "drop") %>%
    mutate(intra_max_vis = ifelse(is.infinite(intra_max_vis), 0, intra_max_vis))
  
  # 4. FLUID DYNAMICS
  fluids_calculated <- fluid_df %>%
    filter(.data$period == "intra_op") %>%
    add_missing_cols(fluid_vars) %>% 
    left_join(metadata_df %>% select(!!sym(id_col), !!sym(WEIGHT_VAR)), by = id_col) %>%
    group_by(!!sym(id_col)) %>%
    summarize(
      tmp_weight     = first(!!sym(WEIGHT_VAR)),
      
      # We use your existing names for sums
      tmp_irrig_in   = sum(!!sym(IRRIGATION_IN_VAR), na.rm = TRUE),
      tmp_irrig_out  = sum(!!sym(IRRIGATION_VAR), na.rm = TRUE),
      tmp_total_in    = sum(!!sym(CRYSTALLOID_VAR), na.rm = TRUE) + 
        sum(!!sym(COLLOID_VAR), na.rm = TRUE) + 
        sum(!!sym(BLOOD_PROD_VAR), na.rm = TRUE),
      
      # Systemic Output
      tmp_total_out   = sum(!!sym(URINE_VAR), na.rm = TRUE) + 
        sum(!!sym(EBL_VAR), na.rm = TRUE),
      
      # --- THE FEATURES (Preserving your original names) ---
      intra_wa_ebl           = if_else(tmp_weight > 0, sum(!!sym(EBL_VAR), na.rm = TRUE) / tmp_weight, 0),
      intra_wa_cryst_vol     = if_else(tmp_weight > 0, sum(!!sym(CRYSTALLOID_VAR), na.rm = TRUE) / tmp_weight, 0),
      intra_wa_transfusion   = if_else(tmp_weight > 0, sum(!!sym(BLOOD_PROD_VAR), na.rm = TRUE) / tmp_weight, 0),
      intra_wa_irrig_vol     = if_else(tmp_weight > 0, tmp_irrig_in / tmp_weight, 0),
      intra_wa_irrig_ret     = if_else(tmp_weight > 0, (tmp_irrig_in - tmp_irrig_out) / tmp_weight, 0),
      intra_wa_net_fluid_bal = if_else(tmp_weight > 0, (tmp_total_in - tmp_total_out) / tmp_weight, 0),
      .groups = "drop"
    )
  
  # 5. FINAL JOIN
  vitals_calculated %>%
    left_join(vaso_calculated, by = id_col) %>%
    left_join(fluids_calculated, by = id_col) %>%
    mutate(across(where(is.numeric), ~replace_na(.x, 0)))
}



#PREOP FEATURES
get_baseline_features <- function(lab_df, vitals_df, id_col, time_col) {
  
  # 1. Process Labs: Handle LONG format
  pre_labs <- lab_df %>%
    filter(.data$period == "pre_op")
  
  # If it's LONG, pivot it WIDE so the rest of the logic works
  if("variable" %in% colnames(pre_labs)) {
    pre_labs <- pre_labs %>%
      pivot_wider(id_cols = c(!!sym(id_col), !!sym(time_col)), 
                  names_from = variable, values_from = value)
  }
  
  pre_labs <- pre_labs %>%
    add_missing_cols(c(PH_VAR, BE_VAR, BICARB_VAR, LACTATE_VAR, NEUTRO_VAR, 
                       PLATE_VAR, MONO_VAR, LYMPHO_VAR, WBC_VAR, CRP_VAR, 
                       FIBRINOGEN_VAR, MPV_VAR, RDW_VAR, HB_VAR, ALBUMIN_VAR, 
                       BUN_VAR, CREAT_VAR)) %>%
    
    arrange(!!sym(id_col), !!sym(time_col)) %>%
    group_by(!!sym(id_col)) %>%
    summarize(
      pre_ph      = if(PH_VAR %in% names(.)) first(na.omit(as.numeric(!!sym(PH_VAR)))) else NA_real_,
      pre_be      = if(BE_VAR %in% names(.)) first(na.omit(as.numeric(!!sym(BE_VAR)))) else NA_real_,
      pre_bicarb  = if(BICARB_VAR %in% names(.)) first(na.omit(as.numeric(!!sym(BICARB_VAR)))) else NA_real_,
      pre_lactate = if(LACTATE_VAR %in% names(.)) first(na.omit(as.numeric(!!sym(LACTATE_VAR)))) else NA_real_,
      pre_neutro  = if(NEUTRO_VAR %in% names(.)) first(na.omit(as.numeric(!!sym(NEUTRO_VAR)))) else NA_real_,
      pre_plate   = if(PLATE_VAR %in% names(.)) first(na.omit(as.numeric(!!sym(PLATE_VAR)))) else NA_real_,
      pre_mono    = if(MONO_VAR %in% names(.)) first(na.omit(as.numeric(!!sym(MONO_VAR)))) else NA_real_,
      pre_lympho  = if(LYMPHO_VAR %in% names(.)) first(na.omit(as.numeric(!!sym(LYMPHO_VAR)))) else NA_real_,
      pre_wbc     = if(WBC_VAR %in% names(.)) first(na.omit(as.numeric(!!sym(WBC_VAR)))) else NA_real_,
      pre_crp     = if(CRP_VAR %in% names(.)) first(na.omit(as.numeric(!!sym(CRP_VAR)))) else NA_real_,
      pre_fib     = if(FIBRINOGEN_VAR %in% names(.)) first(na.omit(as.numeric(!!sym(FIBRINOGEN_VAR)))) else NA_real_,
      pre_mpv     = if(MPV_VAR %in% names(.)) first(na.omit(as.numeric(!!sym(MPV_VAR)))) else NA_real_,
      pre_rdw     = if(RDW_VAR %in% names(.)) first(na.omit(as.numeric(!!sym(RDW_VAR)))) else NA_real_,
      pre_hb      = if(HB_VAR %in% names(.)) first(na.omit(as.numeric(!!sym(HB_VAR)))) else NA_real_,
      pre_alb     = if(ALBUMIN_VAR %in% names(.)) first(na.omit(as.numeric(!!sym(ALBUMIN_VAR)))) else NA_real_,
      tmp_bun     = if(BUN_VAR %in% names(.)) first(na.omit(as.numeric(!!sym(BUN_VAR)))) else NA_real_,
      tmp_creat   = if(CREAT_VAR %in% names(.)) first(na.omit(as.numeric(!!sym(CREAT_VAR)))) else NA_real_,
      .groups = "drop"
    ) %>%
    mutate(pre_bun_creat_ratio = as.numeric(tmp_bun) / pmax(0.1, as.numeric(tmp_creat))) %>%
    select(-starts_with("tmp_"))
  
  # 2. Process Vitals (Same logic)
  pre_vitals <- vitals_df %>%
    filter(period == "pre_op") %>%
    add_missing_cols(c(HR_VAR, MAP_VAR, SPO2_VAR)) %>%
    arrange(!!sym(id_col), !!sym(time_col)) %>%
    group_by(!!sym(id_col)) %>%
    summarize(
      pre_hr   = first(na.omit(!!sym(HR_VAR))),
      pre_map  = first(na.omit(!!sym(MAP_VAR))),
      pre_spo2 = first(na.omit(!!sym(SPO2_VAR))),
      .groups = "drop"
    )
  
  # 3. Final Join
  pre_labs %>%
    full_join(pre_vitals, by = id_col) %>%
    # Fill NAs with 0 only at the very end so the Trend Test can run
    mutate(across(where(is.numeric), ~replace_na(.x, 0)))
}



# --- 5. POSTOPERATIVE FEATURES (Fortified Version) ---
get_postop_features <- function(lab_df, vitals_df, fluid_df, vaso_df, metadata_df, id_col, time_col, ...) {
  
  # MANDATORY PIVOTING
  lab_df    <- ensure_wide(lab_df, id_col, time_col)
  vitals_df <- ensure_wide(vitals_df, id_col, time_col)
  vaso_df   <- ensure_wide(vaso_df, id_col, time_col)
  fluid_df  <- ensure_wide(fluid_df, id_col, time_col)
  
  # A. Define the variables we expect to find
  l_vars <- c(PH_VAR, BE_VAR, LACTATE_VAR, NEUTRO_VAR, LYMPHO_VAR, PLATE_VAR, 
              MONO_VAR, WBC_VAR, CRP_VAR, ALBUMIN_VAR, FIBRINOGEN_VAR, BUN_VAR, CREAT_VAR)
  v_vars <- c(SPO2_VAR, FIO2_VAR, MAP_VAR, HR_VAR)
  f_vars <- c(URINE_VAR, CRYSTALLOID_VAR, COLLOID_VAR, BLOOD_PROD_VAR, EBL_VAR)
  va_vars <- c(EPI_VAR, NORE_VAR, VASO_VAR, DOPE_VAR, DOBU_VAR, MIL_VAR)
  
  # B. Prepare data with Safety Guards
  lab_safe <- lab_df %>% add_missing_cols(l_vars)
  vit_safe <- vitals_df %>% add_missing_cols(v_vars)
  flu_safe <- fluid_df %>% add_missing_cols(f_vars)
  vas_safe <- vaso_df %>% add_missing_cols(va_vars)
  
  # 1. 48-HOUR TRAJECTORIES (Slopes)
  
  # A. Calculate S/F Ratio Slopes (Skip the Pivot!)
  sf_slope_table <- vit_safe %>%
    filter(period %in% c("post_24", "post_48")) %>%
    add_missing_cols(c(SPO2_VAR, FIO2_VAR)) %>%
    mutate(sf_ratio = !!sym(SPO2_VAR) / pmax(0.1, !!sym(FIO2_VAR) / 100)) %>%
    group_by(!!sym(id_col)) %>%
    summarize(
      post_sf_slope = if(n() >= 2) calc_slope(pick(everything()), all_of(time_col), "sf_ratio") else 0,
      .groups = "drop"
    )
  
  # B. Main Trajectory Factory (Swapping cur_data for pick)
  post_slopes <- lab_safe %>%
    filter(.data$period %in% c("post_24", "post_48")) %>%
    group_by(!!sym(id_col)) %>%
    summarize(
      post_ph_slope       = calc_slope(pick(everything()), time_col, PH_VAR),
      post_be_slope       = calc_slope(pick(everything()), time_col, BE_VAR),
      post_lactate_slope  = calc_slope(pick(everything()), time_col, LACTATE_VAR),
      post_nlr_slope      = calc_slope(pick(everything()) %>% mutate(tmp = !!sym(NEUTRO_VAR)/pmax(0.1, !!sym(LYMPHO_VAR))), time_col, "tmp"),
      post_plr_slope      = calc_slope(pick(everything()) %>% mutate(tmp = !!sym(PLATE_VAR)/pmax(0.1, !!sym(LYMPHO_VAR))), time_col, "tmp"),
      post_mlr_slope      = calc_slope(pick(everything()) %>% mutate(tmp = !!sym(MONO_VAR)/pmax(0.1, !!sym(LYMPHO_VAR))), time_col, "tmp"),
      post_wbc_slope      = calc_slope(pick(everything()), time_col, WBC_VAR),
      post_crp_slope      = calc_slope(pick(everything()), time_col, CRP_VAR),
      .groups = "drop"
    )
  
  # 2. 24-HOUR MEDIANS & AUTONOMIC
  post_24h_stats <- lab_safe %>%
    filter(.data$period == "post_24") %>%
    group_by(!!sym(id_col)) %>%
    summarize(
      post_ph_median  = median(!!sym(PH_VAR), na.rm = TRUE),
      post_nlr_median = median(!!sym(NEUTRO_VAR) / pmax(0.1, !!sym(LYMPHO_VAR)), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    full_join(
      vit_safe %>%
        filter(period == "post_24") %>%
        group_by(!!sym(id_col)) %>%
        summarize(
          post_map_cv   = (sd(!!sym(MAP_VAR), na.rm = TRUE) / mean(!!sym(MAP_VAR), na.rm = TRUE)) * 100,
          post_hr_slope = calc_slope(pick(everything()), time_col, HR_VAR),
          .groups = "drop"
        ), by = id_col
    )
  
  # 3. DELTAS (Reusing your clinical delta engine)
  post_deltas <- calc_clinical_delta(lab_safe, id_col, time_col, ALBUMIN_VAR) %>% 
    rename(post_alb_delta = delta) %>%
    full_join(
      calc_clinical_delta(lab_safe, id_col, time_col, FIBRINOGEN_VAR) %>% 
        rename(post_fib_delta = delta), by = id_col
    )
  
  # 4. BUN/CREAT RATIO DELTA (The "Pantry-Check" Version)
  post_bc_delta <- lab_safe %>%
    filter(period %in% c("pre_op", "post_48")) %>%
    group_by(!!sym(id_col), period) %>%
    summarize(
      tmp_bun   = median(!!sym(BUN_VAR), na.rm = TRUE),
      tmp_creat = median(!!sym(CREAT_VAR), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_wider(
      id_cols = !!sym(id_col),
      names_from = period,
      values_from = c(tmp_bun, tmp_creat)
    ) %>%
    # --- SAFETY INJECTOR ---
    # This ensures the columns exist even if the dummy data is missing that period
    add_missing_cols(c("tmp_bun_pre_op", "tmp_bun_post_48", 
                       "tmp_creat_pre_op", "tmp_creat_post_48")) %>%
    # -----------------------
  mutate(
    post_ratio = tmp_bun_post_48 / pmax(0.1, tmp_creat_post_48),
    pre_ratio  = tmp_bun_pre_op / pmax(0.1, tmp_creat_pre_op),
    post_bun_creat_delta = replace_na(post_ratio - pre_ratio, 0)
  ) %>%
    select(!!sym(id_col), post_bun_creat_delta)
  
  # 5. FLUID & RENAL DYNAMICS
  post_fluids <- flu_safe %>%
    left_join(metadata_df %>% select(!!sym(id_col), !!sym(WEIGHT_VAR)), by = id_col) %>%
    group_by(!!sym(id_col)) %>%
    summarize(
      tmp_weight = first(!!sym(WEIGHT_VAR)),
      uop_24_vol = sum(if_else(period == "post_24", !!sym(URINE_VAR), 0), na.rm = TRUE),
      post_wa_uop_hr = if_else(tmp_weight > 0, (uop_24_vol / tmp_weight) / 24, 0),
      total_in_48  = sum(if_else(period %in% c("post_24", "post_48"), 
                                 !!sym(CRYSTALLOID_VAR) + !!sym(COLLOID_VAR) + !!sym(BLOOD_PROD_VAR), 0), na.rm = TRUE),
      total_out_48 = sum(if_else(period %in% c("post_24", "post_48"), 
                                 !!sym(URINE_VAR) + !!sym(EBL_VAR), 0), na.rm = TRUE),
      post_wa_net_fluid_48 = if_else(tmp_weight > 0, (total_in_48 - total_out_48) / tmp_weight, 0),
      .groups = "drop"
    )
  
  # 6. VASOPRESSOR WEANING
  post_weaning <- vas_safe %>%
    filter(period %in% c("post_24", "post_48")) %>%
    filter(!!sym(EPI_VAR) > 0 | !!sym(NORE_VAR) > 0 | !!sym(VASO_VAR) > 0 | 
             !!sym(DOPE_VAR) > 0 | !!sym(DOBU_VAR) > 0 | !!sym(MIL_VAR) > 0) %>%
    group_by(!!sym(id_col)) %>%
    # The Fix: Use if/else to avoid -Inf from max() on empty sets
    summarize(last_pressor_time = if(n() > 0) max(!!sym(time_col), na.rm = TRUE) else NA, .groups = "drop") %>%
    filter(!is.na(last_pressor_time)) %>% # Only keep patients who actually had pressors
    left_join(metadata_df %>% select(!!sym(id_col), !!sym(SURG_END_VAR)), by = id_col) %>%
    mutate(
      post_vaso_wean_hrs = as.numeric(difftime(last_pressor_time, !!sym(SURG_END_VAR), units = "hours")),
      post_vaso_wean_hrs = pmax(0, replace_na(post_vaso_wean_hrs, 0))
    ) %>%
    select(!!sym(id_col), post_vaso_wean_hrs)
  
  # 7. FINAL MASTER MERGE (Now including the SF Slope table)
  post_slopes %>%
    left_join(post_24h_stats, by = id_col) %>%
    left_join(post_deltas, by = id_col) %>%
    left_join(post_bc_delta, by = id_col) %>%
    left_join(post_fluids %>% select(all_of(id_col), post_wa_uop_hr, post_wa_net_fluid_48), by = id_col) %>%
    left_join(post_weaning, by = id_col) %>%
    left_join(sf_slope_table, by = id_col) %>% # <--- The new safe join
    mutate(across(where(is.numeric), ~replace_na(.x, 0)))
}


# Readmission Features
get_readmission_features <- function(readmit_lab_df, readmit_vitals_df, id_col, time_col) {
# --- 1. Fix the Filter to be Case-Insensitive ---
readmit_labs <- readmit_lab_df %>%
  mutate(variable = as.character(variable)) %>%
  # Ensure we match the strings in your mapping variables
  filter(variable %in% c(PH_VAR, BE_VAR, BICARB_VAR, LACTATE_VAR, NEUTRO_VAR, 
                         PLATE_VAR, MONO_VAR, LYMPHO_VAR, WBC_VAR, CRP_VAR, 
                         FIBRINOGEN_VAR, MPV_VAR, RDW_VAR, HB_VAR, ALBUMIN_VAR, 
                         BUN_VAR, CREAT_VAR)) %>%
  pivot_wider(id_cols = any_of(c(id_col, time_col)), 
              names_from = variable, values_from = value)

# --- 2. Fix the Summarize to handle the 'lab_' prefix properly ---
readmit_labs <- readmit_labs %>%
  group_by(!!sym(id_col)) %>%
  arrange(!!sym(time_col), .by_group = TRUE) %>% 
  summarize(
    readmit_ph      = if(PH_VAR %in% names(.)) first(na.omit(as.numeric(!!sym(PH_VAR)))) else NA_real_,
    readmit_be      = if(BE_VAR %in% names(.)) first(na.omit(as.numeric(!!sym(BE_VAR)))) else NA_real_,
    readmit_bicarb  = if(BICARB_VAR %in% names(.)) first(na.omit(as.numeric(!!sym(BICARB_VAR)))) else NA_real_,
    readmit_lactate = if(LACTATE_VAR %in% names(.)) first(na.omit(as.numeric(!!sym(LACTATE_VAR)))) else NA_real_,
    readmit_neutro  = if(NEUTRO_VAR %in% names(.)) first(na.omit(as.numeric(!!sym(NEUTRO_VAR)))) else NA_real_,
    readmit_plate   = if(PLATE_VAR %in% names(.)) first(na.omit(as.numeric(!!sym(PLATE_VAR)))) else NA_real_,
    readmit_mono    = if(MONO_VAR %in% names(.)) first(na.omit(as.numeric(!!sym(MONO_VAR)))) else NA_real_,
    readmit_lympho  = if(LYMPHO_VAR %in% names(.)) first(na.omit(as.numeric(!!sym(LYMPHO_VAR)))) else NA_real_,
    readmit_wbc     = if(WBC_VAR %in% names(.)) first(na.omit(as.numeric(!!sym(WBC_VAR)))) else NA_real_,
    readmit_crp     = if(CRP_VAR %in% names(.)) first(na.omit(as.numeric(!!sym(CRP_VAR)))) else NA_real_,
    readmit_fib     = if(FIBRINOGEN_VAR %in% names(.)) first(na.omit(as.numeric(!!sym(FIBRINOGEN_VAR)))) else NA_real_,
    readmit_mpv     = if(MPV_VAR %in% names(.)) first(na.omit(as.numeric(!!sym(MPV_VAR)))) else NA_real_,
    readmit_rdw     = if(RDW_VAR %in% names(.)) first(na.omit(as.numeric(!!sym(RDW_VAR)))) else NA_real_,
    readmit_hb      = if(HB_VAR %in% names(.)) first(na.omit(as.numeric(!!sym(HB_VAR)))) else NA_real_,
    readmit_alb     = if(ALBUMIN_VAR %in% names(.)) first(na.omit(as.numeric(!!sym(ALBUMIN_VAR)))) else NA_real_,
    tmp_bun         = if(BUN_VAR %in% names(.)) first(na.omit(as.numeric(!!sym(BUN_VAR)))) else NA_real_,
    tmp_creat       = if(CREAT_VAR %in% names(.)) first(na.omit(as.numeric(!!sym(CREAT_VAR)))) else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(
    readmit_bun_creat_ratio = if_else(!is.na(tmp_bun) & !is.na(tmp_creat), 
                                      as.numeric(tmp_bun) / pmax(0.1, as.numeric(tmp_creat)), 
                                      NA_real_)
  ) %>%
  select(-starts_with("tmp_"))
  
  # 2. Vitals: UPDATED for Wide Format compatibility
  readmit_vitals <- readmit_vitals_df 
  
  # If it's still LONG (has a 'variable' column), pivot it. 
  # If it's already WIDE (from your new dummy), skip this.
  if("variable" %in% colnames(readmit_vitals)) {
    readmit_vitals <- readmit_vitals %>%
      pivot_wider(id_cols = c(!!sym(id_col), !!sym(time_col)), 
                  names_from = variable, values_from = value)
  }
  
  readmit_vitals <- readmit_vitals %>%
    add_missing_cols(unique(c(HR_VAR, MAP_VAR, SPO2_VAR))) %>%
    group_by(!!sym(id_col)) %>%
    arrange(!!sym(time_col), .by_group = TRUE) %>%
    summarize(
      readmit_hr   = if(HR_VAR %in% names(.)) first(na.omit(as.numeric(!!sym(HR_VAR)))) else NA_real_,
      readmit_map  = if(MAP_VAR %in% names(.)) first(na.omit(as.numeric(!!sym(MAP_VAR)))) else NA_real_,
      readmit_spo2 = if(SPO2_VAR %in% names(.)) first(na.omit(as.numeric(!!sym(SPO2_VAR)))) else NA_real_,
      .groups = "drop"
    )
  
  # 3. Join and Clean
  readmit_labs %>%
    full_join(readmit_vitals, by = id_col)
}