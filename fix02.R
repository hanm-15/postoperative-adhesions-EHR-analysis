library(tidyverse)
library(lubridate)
source("01.R")
source("saving02.R")

#Fix replace_na with mice later

# --- Encounter Definition ---
labeled_data <- raw_data %>%
  group_by(!!sym(ENCOUNTER_ID_VAR)) %>%
  mutate(
    # We use grepl here because it proved more stable with !!sym()
    intraperitoneal_confirm = grepl(intraperitoneal_surg, !!sym(PCS_VAR)),
    
    lysis_confirm = grepl(lysis_proc, !!sym(PCS_VAR)) & 
      grepl(adhesion_diag, !!sym(DX_VAR)),
    
    readmission_confirm = grepl(adhesion_diag, !!sym(DX_VAR)) & 
      (!grepl(surgical_intervention, !!sym(PCS_VAR)) | is.na(!!sym(PCS_VAR)))
  ) %>%
  ungroup() %>%
  # Fill NAs to prevent logic crashes
  mutate(across(ends_with("_confirm"), ~replace_na(.x, FALSE)))


# --- Exclusion Definitions ---
excl_cancer_rad <- "^(C7[789]|C80|Z923).*"
excl_insult <- "^([SVWXY]|T(?!(8[0-8]))).*"
excl_barrier <- paste0("^(0D.[678ABEFGHJLNUVW].[JK7]|",
                       "0F.[012456789].[JK7]|",
                       "07.[BP].[JK7]|",
                       "0U.[01245679CF].[JK7]).*")

# --- Patient Filtering
# 'raw data' is the placeholder
final_cohort <- labeled_data %>% 
  group_by(!!sym(PATIENT_ID_VAR)) %>% 
  filter(any(intraperitoneal_confirm == TRUE)) %>% 
  filter(all(!!sym(AGE_VAR) >= 18)) %>% 
  # Use %in% TRUE or handle NA explicitly
  filter(!any(str_detect(!!sym(DX_VAR), excl_cancer_rad) %in% TRUE)) %>% 
  filter(!any(str_detect(!!sym(DX_VAR), excl_insult) %in% TRUE)) %>% 
  filter(!any(str_detect(!!sym(PCS_VAR), excl_barrier) %in% TRUE)) %>% 
  ungroup()

# --- Index Surg Definition

if(nrow(final_cohort) > 0) {
  
  # Step A: Define Index Date (Using the stable temp_date method)
  patient_grouping <- final_cohort %>%
    group_by(!!sym(PATIENT_ID_VAR)) %>%
    arrange(!!sym(DATE_VAR)) %>%
    mutate(
      is_index_surgery = (intraperitoneal_confirm == TRUE & !duplicated(intraperitoneal_confirm)),
      temp_date = if_else(intraperitoneal_confirm == TRUE, !!sym(DATE_VAR), as.Date(NA))
    ) %>%
    mutate(
      index_date = min(temp_date, na.rm = TRUE)
    ) %>%
    select(-temp_date) %>% 
    ungroup()
  
  
  # --- Patient Grouping
  analysis_data <- patient_grouping %>%
    group_by(!!sym(PATIENT_ID_VAR)) %>%
    mutate(
      patient_group = case_when(
        any(lysis_confirm == TRUE) ~ "reoperated",
        any(readmission_confirm == TRUE) & !any(lysis_confirm == TRUE) ~ "non_surgical",
        any(intraperitoneal_confirm == TRUE) & !any(lysis_confirm == TRUE) & 
          !any(readmission_confirm == TRUE) ~ "asymptomatic",
        TRUE ~ "other_or_excluded"
      ),
      
      # Identify the subset logic
      post_index_surg = (intraperitoneal_confirm == TRUE & !!sym(DATE_VAR) > index_date),
      immediate_lysis = if_else(any(post_index_surg), 
                                lysis_confirm[which(post_index_surg)[1]] == TRUE, 
                                FALSE),
      reoperated_subset = (patient_group == "reoperated" & immediate_lysis == TRUE)
    ) %>%
    ungroup()
  
  message("!!! SUCCESS: analysis_data created !!!")
  
  # This provides a much clearer breakdown of the encounter-level logic
  analysis_data %>%
    dplyr::count(patient_group, is_index_surgery, lysis_confirm) %>% # Explicitly use dplyr
    mutate(label = case_when(
      is_index_surgery ~ "Index Surgeries",
      lysis_confirm    ~ "Lysis Operations",
      TRUE             ~ "Other/Excl"
    )) %>%
    select(patient_group, label, n) %>%
    print()
  
} else {
  stop("CRITICAL: final_cohort is empty. Check your filters.")
}