library(tidyverse)
library(tidyLPA)
library(clinfun)
library(ggrepel)
source("fix06.R")

# --- 1. DERIVE ANCHOR WEIGHTS (LYSIS DATA) ---
message("Deriving Robust Anchor Weights...")

# A. Capture Global Baseline for Lysis Population
anchor_pop_stats <- lca_features %>%
  summarise(across(all_of(viable_cols), 
                   list(m = ~mean(.x, na.rm = TRUE), s = ~sd(.x, na.rm = TRUE))))

# B. Derive and Standardize
anchor_weights <- lca_features %>%
  inner_join(final_assignments, by = PATIENT_ID_VAR) %>%
  filter(as.numeric(as.character(Class)) > 0) %>% 
  group_by(Class) %>%
  summarise(across(all_of(viable_cols), 
                   list(mean = ~mean(.x, na.rm = TRUE), 
                        sd = ~sd(.x, na.rm = TRUE)), 
                   .names = "{col}_{fn}")) %>%
  pivot_longer(cols = -Class, names_to = "metric", values_to = "raw_val") %>%
  mutate(
    base_var = str_remove(metric, "_(mean|sd)$"),
    type = str_extract(metric, "(mean|sd)$")
  ) %>%
  rowwise() %>%
  mutate(
    m_ref = anchor_pop_stats[[paste0(base_var, "_m")]],
    s_ref = anchor_pop_stats[[paste0(base_var, "_s")]]
  ) %>%
  mutate(z_val = case_when(
    type == "mean" ~ (raw_val - m_ref) / s_ref,
    type == "sd"   ~ raw_val / s_ref,
    TRUE ~ raw_val
  )) %>%
  ungroup() %>%
  select(Class, metric, z_val) %>%
  pivot_wider(names_from = metric, values_from = z_val)


# --- 2. DERIVE ACUTE WEIGHTS (INDEX SURGERY) ---
message("Deriving Robust Acute Weights...")

# This looks for any string starting with 'intra_' OR 'post_'
acute_vars <- str_subset(available_features, "^(intra_|post_)")

# A. Capture Global Baseline for Index Population
acute_pop_stats <- index_features %>%
  summarise(across(all_of(acute_vars), 
                   list(m = ~mean(.x, na.rm = TRUE), s = ~sd(.x, na.rm = TRUE))))

# B. Derive and Standardize
acute_weights <- index_features %>%
  group_by(Class) %>%
  summarise(across(all_of(acute_vars), 
                   list(mean = ~mean(.x, na.rm = TRUE), 
                        sd = ~sd(.x, na.rm = TRUE)), 
                   .names = "{col}_{fn}")) %>%
  pivot_longer(cols = -Class, names_to = "metric", values_to = "raw_val") %>%
  mutate(
    base_var = str_remove(metric, "_(mean|sd)$"),
    type = str_extract(metric, "(mean|sd)$")
  ) %>%
  rowwise() %>%
  mutate(
    m_ref = acute_pop_stats[[paste0(base_var, "_m")]],
    s_ref = acute_pop_stats[[paste0(base_var, "_s")]]
  ) %>%
  mutate(z_val = case_when(
    type == "mean" ~ (raw_val - m_ref) / s_ref,
    type == "sd"   ~ raw_val / s_ref,
    TRUE ~ raw_val
  )) %>%
  ungroup() %>%
  select(Class, metric, z_val) %>%
  pivot_wider(names_from = metric, values_from = z_val)

# --- FINAL INTERPRETATION & VISUALIZATION ---

message("\n--- PHENOTYPIC WEIGHT SUMMARY ---")
message("Acute-Response Features found: ", length(acute_vars))
message("Anchor (Lysis) Features found: ", length(viable_cols))

# 1. Console Summary: Compare a key feature across classes
# Let's pick a representative intra-op and post-op feature to peek at
example_vars <- head(acute_vars, 3) 
acute_summary_table <- acute_weights %>%
  select(Class, contains(example_vars[1]), contains(example_vars[2]))

print(acute_summary_table)

# 2. Visual Profile Plot (The "Lollipop" of Means)
# This helps you see if Class 0 is truly different from the others
plot_weights <- acute_weights %>%
  pivot_longer(cols = -Class, names_to = "metric", values_to = "value") %>%
  filter(str_detect(metric, "_mean$")) %>%
  mutate(feature = str_remove(metric, "_mean$"))

# Create the visualization
weight_plot <- ggplot(plot_weights, aes(x = feature, y = value, color = Class, group = Class)) +
  geom_line(size = 1, alpha = 0.5) +
  geom_point(size = 3) +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Acute-Response Phenotypic Profiles",
    subtitle = "Mean values across Class 0 (Asymptomatic) vs. Latent Classes",
    x = "Clinical Feature",
    y = "Standardized Mean Value",
    color = "Phenotype"
  ) +
  theme(legend.position = "bottom")

print(weight_plot)

message("\n!!! SUCCESS: Weight derivation complete. All systems stable for the night. !!!")