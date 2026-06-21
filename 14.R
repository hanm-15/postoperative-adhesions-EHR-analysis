source("LTA_test.R")
library(tidyLPA)
library(dplyr)

message("\n--- STARTING SECTION 3.D: LMR & ENTROPY STABILITY AUDIT ---")

# 1. Grab the RAW clinical values we just preserved
lpa_audit_data <- lta_input_copy %>% 
  filter(node_label == "lysis_confirm") %>%
#  select(all_of(c(chronic_vars_final, anchor_vars_final, acute_vars_final))) %>%
  select(all_of(c(anchor_vars_final))) %>%
  drop_na()

# 2. Local Scaling: Stretch the data to its own Max/Min
# This ensures Class 1 and Class 2 are perfectly visible to the LCA
lpa_audit_data <- lpa_audit_data %>%
  mutate(across(everything(), ~ as.numeric(scale(.))))

# --- NEW: Run the Actual Models first! ---
# This creates the 'all_models' list that was missing
message("Running expansion models (K=2 to K=5)...")
all_models <- lapply(2:5, function(k) {
  lpa_audit_data %>% 
    estimate_profiles(k, 
                      package = "mclust", 
                      control = mclust::emControl(itmax = c(500, 500))) # <--- THE "FORCE" BUTTON
})

# This creates the 'all_fits' list for Entropy checks
all_fits <- lapply(all_models, get_fit)


# --- SECTION 3.D: SEQUENTIAL STEP-UP AUDIT ---

# Function to calculate the p-value for the 'Step Up' (K vs K-1)
calc_lrt <- function(m_curr, m_prev) {
  f_curr <- get_fit(m_curr)
  f_prev <- get_fit(m_prev)
  lr_diff <- 2 * (f_curr$LogLik - f_prev$LogLik)
  df_diff <- f_curr$KIC - f_prev$KIC
  if (is.na(lr_diff) || lr_diff <= 0) return(1.0)
  return(1 - pchisq(lr_diff, df = df_diff))
}

# 1. Run the Step-Up Comparisons
# Note: all_models[[1]] is K=2, all_models[[2]] is K=3, etc.
p_2to3 <- calc_lrt(all_models[[2]], all_models[[1]])
p_3to4 <- calc_lrt(all_models[[3]], all_models[[2]])
p_4to5 <- calc_lrt(all_models[[4]], all_models[[3]])

# 2. Sequential "Success Gate" (BANDAGED)
# isTRUE forces any NA or NaN to return FALSE, preventing the crash
is_3_better_than_2 <- isTRUE(p_2to3 < 0.05 & all_fits[[2]]$Entropy >= 0.80)
is_4_better_than_3 <- isTRUE(p_3to4 < 0.05 & all_fits[[3]]$Entropy >= 0.80)
is_5_better_than_4 <- isTRUE(p_4to5 < 0.05 & all_fits[[4]]$Entropy >= 0.80)

message("--- SEQUENTIAL STABILITY REPORT ---")
# A small helper function to handle the math-gore
safe_p <- function(p) {
  if (is.nan(p)) return("NaN")
  return(format.pval(p, digits = 4, eps = 0.001))
}

message(paste("Step 2->3 (LMR p =", safe_p(p_2to3), "): Significant?", is_3_better_than_2))
message(paste("Step 3->4:(LMR p =", safe_p(p_3to4), "): Significant?", is_4_better_than_3))
message(paste("Step 4->5:(LMR p =", safe_p(p_4to5), "): Significant?", is_5_better_than_4))

# 3. Final Verdict (BANDAGED)
if (!is_3_better_than_2) {
  message("--- AUDIT SUCCESS: No evidence to move beyond 2nd-class Legacy. ---")
  use_legacy <- TRUE
} else {
  message("--- AUDIT NOTICE: Expansion possible. Proceeding with Legacy for Invariance. ---")
  use_legacy <- TRUE
}


# --- SECTION 3.D: THE TRANSPARENCY UPDATE ---
options(scipen = 999) # Stop scientific notation (show 0.00001 instead of 1e-05)

message("\n--- RAW STATISTICAL BREAKDOWN ---")

# Look at the raw Log-Likelihoods to see the 'Information Gain'
for(i in 1:4) {
  k_val <- i + 1
  fit <- all_fits[[i]]
  message(paste0("K=", k_val, " | LogLik: ", round(fit$LogLik, 2), 
                 " | BIC: ", round(fit$BIC, 2), 
                 " | Entropy: ", round(fit$Entropy, 4)))
}

message("\n--- DETAILED STEP-UP COMPARISONS ---")
message(paste("Step 2->3: raw p-value =", format.pval(p_2to3, digits = 10, eps = 0.0000000001)))
message(paste("Step 3->4: raw p-value =", format.pval(p_3to4, digits = 10, eps = 0.0000000001)))
message(paste("Step 4->5: raw p-value =", format.pval(p_4to5, digits = 10, eps = 0.0000000001)))

# --- 4. CLASS DISTRIBUTION BREAKDOWN (FORCED PRINT) ---
message("\n--- CLASS DISTRIBUTION BREAKDOWN ---")

for(i in 1:length(all_models)) {
  k_val <- i + 1
  
  # Use tidyLPA's built-in assignment tool to get the N counts
  # This avoids the "Object not found" or "Skipped" errors
  assignments <- tryCatch({
    get_data(all_models[[i]])
  }, error = function(e) return(NULL))
  
  if(!is.null(assignments)) {
    # tidyLPA usually names the column 'Class'
    dist_table <- assignments %>%
      count(Class) %>%
      mutate(Percentage = (n / sum(n)) * 100)
    
    message(paste0("\n--- Results for K = ", k_val, " ---"))
    for(j in 1:nrow(dist_table)) {
      message(paste0("  Class ", dist_table$Class[j], ": n = ", dist_table$n[j], 
                     " (", round(dist_table$Percentage[j], 1), "%)"))
    }
  } else {
    message(paste0("\n--- K = ", k_val, ": Could not retrieve assignments ---"))
  }
}


if("package:mclust" %in% search()) detach("package:mclust", unload=TRUE)