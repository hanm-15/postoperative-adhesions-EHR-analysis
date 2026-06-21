library(tidyverse)
library(lubridate)

# --- GLOBAL VARS RETAINED ---
# Global variables
DATE_VAR <- "surg_start"
ADMIT_TIME_VAR <- "admission_time"
ENCOUNTER_ID_VAR <- "enc_id"
PATIENT_ID_VAR   <- "pat_id"
PCS_VAR          <- "pcs_proc"
DX_VAR           <- "cm_diag"
DURATION_VAR     <- "weight_kg" 
WEIGHT_VAR       <- "weight_kg"
EBL_VAR          <- "intra_ebl"
AGE_VAR          <- "age_yrs"
CRYSTALLOID_VAR  <- "fluid_cryst_ml"
COLLOID_VAR      <- "intra_colloid"
BLOOD_PROD_VAR   <- "intra_blood"
IRRIGATION_VAR   <- "intra_irrig_out"
IRRIGATION_IN_VAR <- "intra_irrig_in"
SVV_VAR <- "stroke_vol"
DOPE_VAR         <- "med_dopamine"
DOBU_VAR         <- "med_dobutamine"
EPI_VAR          <- "med_epinephrine"
NORE_VAR         <- "med_norepinephrine"
MIL_VAR          <- "med_milrinone"
VASO_VAR         <- "med_vasopressin"
PATIENT_ID_VAR    <- "pat_id"
ENCOUNTER_ID_VAR  <- "enc_id"
CHART_TIME_VAR    <- "chart_time"
SURG_START_VAR    <- "surg_start"
SURG_END_VAR      <- "surg_end"
URINE_VAR <- "urine_ml"
DURATION_VAR <- "dur_min"
ICU_IN_VAR       <- "icu_admit_time"
ICU_OUT_VAR      <- "icu_exit_time"
DISCHARGE_VAR    <- "discharge_time"
URGENCY_VAR <- "urgency_status"
HEIGHT_VAR <- "hight_cm"
GENDER_VAR <- "gender"
ASA_VAR <- "asa_score"

# Clinical Mappings
MAP_VAR          <- "map_val" 
HR_VAR           <- "hr_val"
SPO2_VAR      <- "spo2_val"
FIO2_VAR      <- "fio2_val"
ETCO2_VAR     <- "etco2_val"
TEMP_VAR      <- "temp_val"

# Lab Mappings
PH_VAR        <- "lab_ph"
LACTATE_VAR   <- "lab_lactate"
WBC_VAR       <- "lab_wbc"
HB_VAR        <- "lab_hb"
ALBUMIN_VAR   <- "lab_alb"
CRP_VAR       <- "lab_crp"
BE_VAR        <- "lab_be"
BICARB_VAR    <- "lab_bicarb"
NEUTRO_VAR    <- "lab_neutro"
LYMPHO_VAR    <- "lab_lympho"
MONO_VAR      <- "lab_mono"
PLATE_VAR     <- "lab_plate"
FIBRINOGEN_VAR <- "lab_fib"
MPV_VAR       <- "lab_mpv"
RDW_VAR       <- "lab_rdw"
BUN_VAR       <- "lab_bun"
CREAT_VAR     <- "lab_creat"

# --- 1. THE ASYMMETRIC COHORTS (EXPANDED) ---
reop_ids     <- as.character(101:250) 
non_surg_ids <- as.character(251:300)


raw_data <- bind_rows(
  tibble(pat_id = reop_ids, patient_group = "reoperated", reoperated_gen = TRUE),
  tibble(pat_id = non_surg_ids, patient_group = "non_surgical", reoperated_gen = FALSE)
) %>%
  uncount(2) %>% 
  group_by(pat_id) %>%
  mutate(encounter_num = row_number()) %>%
  ungroup() %>%
  mutate(enc_id = as.character(row_number())) %>% 
  mutate(
    is_index_surgery = encounter_num == 1,
    
    # --- PROBABILITY-BASED PROFILES ---
    profile = case_when(
      !is_index_surgery & reoperated_gen ~ if_else(runif(n()) < 0.45, "Sick", "Healthy"),
      is_index_surgery & (runif(n()) < 0.25) ~ "Sick",
      !is_index_surgery & !reoperated_gen & (runif(n()) < 0.2) ~ "Sick",
      TRUE ~ "Healthy"
    ),
    
    # --- NOISY METADATA (Driven by the Profile above) ---
    age_yrs   = pmax(18, rnorm(n(), 55, 12)), 
    weight_kg = rnorm(n(), 80, 15), 
    dur_min   = runif(n(), 120, 360),
    !!sym(ASA_VAR) := case_when(
      profile == "Sick" ~ sample(c(2, 3, 4), n(), replace = TRUE, prob = c(0.2, 0.5, 0.3)),
      TRUE              ~ sample(c(1, 2, 3), n(), replace = TRUE, prob = c(0.4, 0.5, 0.1))
    ),
    !!sym(GENDER_VAR) := sample(c("M", "F"), n(), replace = TRUE),
    !!sym(HEIGHT_VAR) := if_else(
      !!sym(GENDER_VAR) == "M", 
      rnorm(n(), 175, 7), 
      rnorm(n(), 162, 7)
    ),
    !!sym(URGENCY_VAR) := case_when(
      profile == "Sick" ~ if_else(runif(n()) < 0.7, "Emergency", "Urgent"),
      TRUE ~ "Elective"
    ),
    # --- DATES & CODES ---
    event_start = if_else(is_index_surgery, ymd_hm("2020-01-01 08:00"), ymd_hm("2020-07-01 08:00")),
    
    # Only assign a surgery start if it's an Index or a Reop-subset second encounter
    # Non-surgical readmissions get NA for surgery times
    surg_start = case_when(
      is_index_surgery ~ event_start,
      !is_index_surgery & reoperated_gen ~ event_start,
      TRUE ~ as.POSIXct(NA) 
    ),
    surg_end = surg_start + minutes(as.integer(dur_min)),
    admission_time = event_start - hours(2),
    discharge_time = event_start + days(5),
    
    # ICD/PCS logic
    lysis_confirm = !is_index_surgery & reoperated_gen,
    pcs_proc = case_when(
      is_index_surgery ~ "0D060ZZ", 
      lysis_confirm    ~ "0DN60ZZ", 
      TRUE             ~ NA_character_ # Non-surgical readmissions have no PCS code
    ),
    cm_diag = if_else(is_index_surgery, "K255", "K660")
  )

# Extract the map for the measurements
severity_map <- raw_data %>% select(enc_id, profile, event_start)


# --- 2. THE LONG-FORMAT MEASUREMENTS (BALANCED SEVERITY) ---
set.seed(500)

master_measurements <- severity_map %>%
  uncount(48) %>% # Increased to 48 observations
  group_by(enc_id) %>%
  mutate(
    obs = row_number(),
    # Generate measurements every 1 hour, starting 12 hours before surgery
    chart_time = event_start + hours(obs - 1) - hours(12),
    
    # Generate distinct internal names so we can pivot them later
    map_internal = if_else(profile == "Sick", rnorm(n(), 55, 5), rnorm(n(), 85, 5)),
    hr_internal  = if_else(profile == "Sick", rnorm(n(), 115, 10), rnorm(n(), 75, 8)),
    
    # --- OTHER CLINICAL DATA ---
    !!sym(SPO2_VAR)  := rnorm(n(), 94, 2),
    !!sym(FIO2_VAR) := if_else(profile == "Sick", 
                               rnorm(n(), 75, 30), # Most sick patients now > 0.6
                               rnorm(n(), 50, 5)), # Healthy patients stay low
    !!sym(ETCO2_VAR) := rnorm(n(), 35, 3),
    !!sym(TEMP_VAR)  := rnorm(n(), 37, 0.5),
    !!sym(SVV_VAR)   := if_else(profile == "Sick", rnorm(n(), 18, 3), rnorm(n(), 10, 2)),
    
    # --- LABS ---
    !!sym(PH_VAR)      := if_else(profile == "Sick", rnorm(n(), 7.22, 0.05), rnorm(n(), 7.41, 0.02)),
    !!sym(LACTATE_VAR) := if_else(profile == "Sick", rnorm(n(), 5.5, 1.2), rnorm(n(), 1.1, 0.2)),
    !!sym(BE_VAR)     := if_else(profile == "Sick", rnorm(n(), -6, 2), rnorm(n(), 0, 1)),
    !!sym(BICARB_VAR) := if_else(profile == "Sick", rnorm(n(), 18, 2), rnorm(n(), 25, 1)),
    !!sym(WBC_VAR)     := if_else(profile == "Sick", rnorm(n(), 18, 3), rnorm(n(), 8, 1)),
    !!sym(HB_VAR)      := rnorm(n(), 12, 1.5), 
    !!sym(CRP_VAR)     := if_else(profile == "Sick", rnorm(n(), 50, 10), rnorm(n(), 5, 2)),
    !!sym(ALBUMIN_VAR) := if_else(profile == "Sick", rnorm(n(), 2.8, 0.3), rnorm(n(), 4.0, 0.3)),
    !!sym(NEUTRO_VAR)  := rnorm(n(), 70, 5),
    # --- DYNAMIC HEMATOLOGY & CHEMISTRY ---
    !!sym(PLATE_VAR)   := if_else(profile == "Sick", rnorm(n(), 140, 30), rnorm(n(), 250, 40)),
    !!sym(LYMPHO_VAR)  := rnorm(n(), 15, 5),
    !!sym(MONO_VAR)    := rnorm(n(), 8, 2),
    !!sym(FIBRINOGEN_VAR) := if_else(profile == "Sick", rnorm(n(), 450, 50), rnorm(n(), 300, 30)),
    !!sym(MPV_VAR)     := if_else(profile == "Sick", rnorm(n(), 12.5, 1), rnorm(n(), 9.5, 1)),
    !!sym(RDW_VAR)     := rnorm(n(), 14, 1.2),
    
    # --- RENAL FUNCTION (Highly Predictive) ---
    !!sym(BUN_VAR)     := if_else(profile == "Sick", rnorm(n(), 35, 8), rnorm(n(), 15, 3)),
    !!sym(CREAT_VAR)   := if_else(profile == "Sick", rnorm(n(), 1.8, 0.4), rnorm(n(), 0.9, 0.2)),
    
    # --- URINE & FLUIDS ---
    !!sym(URINE_VAR) := if_else(profile == "Sick", rnorm(n(), 100, 40), rnorm(n(), 800, 150)),
    !!sym(CRYSTALLOID_VAR) := if_else(profile == "Sick", rnorm(n(), 2500, 500), rnorm(n(), 500, 100)),
    !!sym(EBL_VAR)  := if_else(profile == "Sick", rnorm(n(), 800, 200), rnorm(n(), 100, 50)), # Ensures EBL is never 0
    !!sym(IRRIGATION_IN_VAR)  := if_else(profile == "Sick", rnorm(n(), 2000, 300), rnorm(n(), 500, 100)),
    !!sym(IRRIGATION_VAR)     := if_else(profile == "Sick", rnorm(n(), 1800, 200), rnorm(n(), 450, 50)),
    !!sym(BLOOD_PROD_VAR)     := if_else(profile == "Sick", rnorm(n(), 500, 100), 0),
    !!sym(COLLOID_VAR)        := if_else(profile == "Sick", rnorm(n(), 500, 50), 0),
    
    # --- VASOPRESSORS (The Fix) ---
    # We generate these as columns first to avoid the 'value' not found error
    med_norepinephrine = if_else(profile == "Sick" & obs > 12 & obs < 20, runif(n(), 0.05, 0.2), 0),
    med_vasopressin    = if_else(profile == "Sick" & obs > 12 & obs < 20, 0.04, 0),
    med_epinephrine    = 0, med_dopamine = 0, med_dobutamine = 0, med_milrinone = 0
  ) %>% ungroup()

# --- 3. THE PIPELINE SUBSETS (WIDE VITALS ALIGNMENT) ---

# 1. VITALS & HEMODYNAMICS (Keep this WIDE)
vitals_ext <- master_measurements %>%
  rename(!!sym(MAP_VAR) := map_internal, 
         !!sym(HR_VAR)  := hr_internal) %>%
  select(enc_id, chart_time, !!sym(MAP_VAR), !!sym(HR_VAR), 
         !!sym(SVV_VAR), # <--- Added this
         all_of(c(SPO2_VAR, FIO2_VAR, ETCO2_VAR, TEMP_VAR)))

# 2. LAB RESULTS (Keep this LONG)
lab_ext <- master_measurements %>%
  select(enc_id, chart_time, starts_with("lab_")) %>%
  pivot_longer(cols = starts_with("lab_"), 
               names_to = "variable", # Your functions look for 'variable'
               values_to = "value")

# 3. VASOPRESSORS (Fixed with Summarize)
vaso_ext <- master_measurements %>%
  select(enc_id, chart_time, starts_with("med_")) %>%
  pivot_longer(cols = starts_with("med_"), names_to = "variable", values_to = "value")

# Fluids (Long)
# Update this section in your generator:
fluid_ext <- master_measurements %>%
  select(enc_id, chart_time, 
         !!sym(CRYSTALLOID_VAR), !!sym(EBL_VAR), !!sym(URINE_VAR),
         !!sym(IRRIGATION_IN_VAR), !!sym(IRRIGATION_VAR), 
         !!sym(BLOOD_PROD_VAR), !!sym(COLLOID_VAR)) %>% # <--- Added these
  pivot_longer(cols = -c(enc_id, chart_time), names_to = "variable", values_to = "value") %>%
  group_by(enc_id, chart_time, variable) %>%
  summarize(value = max(value, na.rm = TRUE), .groups = "drop")

# --- 4. THE FINAL ALIGNMENT ---

timing_ext <- raw_data %>%
  select(!!sym(PATIENT_ID_VAR), !!sym(ENCOUNTER_ID_VAR), 
         !!sym(SURG_START_VAR), !!sym(SURG_END_VAR), admission_time) %>%
  mutate(
    # If there is no surgery, use Admission Time as the anchor (t0)
    t0 = if_else(is.na(!!sym(SURG_START_VAR)), admission_time, !!sym(SURG_START_VAR)),
    # For non-surgical, t1 can just be t0 + 4 hours (a dummy window)
    t1 = if_else(is.na(!!sym(SURG_END_VAR)), t0 + hours(4), !!sym(SURG_END_VAR))
  ) %>%
  select(!!sym(PATIENT_ID_VAR), !!sym(ENCOUNTER_ID_VAR), t0, t1)

raw_surgery_ext <- raw_data

message("External dataset created")