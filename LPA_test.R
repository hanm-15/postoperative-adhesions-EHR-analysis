library(tidyverse)
library(tidyLPA)
library(ggrepel)
source("fix04.R")

# 1. Prepare the Data
# tidyLPA works best with a clean data frame of your Z-scores
lca_data <- lca_matrix 

# 2. Run the LPA Models
# In the real study, you might test 1:6.
message("Estimating Latent Profiles...")
model_results <- lca_data %>%
  estimate_profiles(1:2, 
                    package = "mclust", 
                    control = mclust::emControl(itmax = c(500, 500)))

# 3. Model Fit Statistics
# This table xhows BIC, SABIC, and Entropy.
# REMEMBER: You want the lowest BIC and Entropy closest to 1.0.
stats <- get_fit(model_results)
print(stats)

# 4. Visualize the BIC 'Elbow' to choose the best K
ggplot(stats, aes(x = Classes, y = BIC)) +
  geom_line(group = 1, color = "steelblue") +
  geom_point(size = 3, color = "darkblue") +
  labs(title = "Model Selection: BIC by Number of Classes",
       x = "Number of Phenotypes (K)",
       y = "BIC (Lower is Better)") +
  theme_minimal()

# --- NEW: COMPOSITE SEVERITY RANKING (FIXED) ---
# 1. Calculate the 'Composite Deviation' for each class
severity_map <- get_estimates(model_results) %>%
  filter(Category == "Means") %>%
  group_by(Class) %>%
  summarize(composite_deviation = mean(abs(Estimate), na.rm = TRUE)) %>%
  arrange(composite_deviation) %>%
  mutate(New_Ordered_Class = row_number()) %>%
  # FIX: Convert Class to character so it matches the other table
  mutate(Class = as.character(Class)) %>% 
  select(Class, New_Ordered_Class)

# 2. Update the assignments bridge with the new labels
final_assignments <- get_data(model_results) %>%
  # Ensure Class is character here too
  mutate(Class = as.character(Class)) %>%
  left_join(severity_map, by = "Class") %>%
  mutate(Class = factor(New_Ordered_Class, 
                        levels = sort(unique(New_Ordered_Class)), 
                        ordered = TRUE)) %>%
  select(-New_Ordered_Class)

# --- THE SIMPLEST PRINT ---
dist_summary <- final_assignments %>%
  count(Class) %>%
  mutate(pct = round(n / sum(n) * 100, 1))
# This prints just the lines you want
pwalk(dist_summary, ~message("Phenotype ", ..1, ": (", ..3, "%)"))

# 5. Visualize the Phenotypes
plot_profiles(model_results)