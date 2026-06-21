library(tidyverse)
library(lubridate)
# --- Intraperitoneal Surg Definition ---
gi_intra <- "0D[018BJNQSTW][678ABEFGHJLNUVW][04].*"
hepato_intra <- "0F[018BJNQTW][012456789][04].*"
lymph_intra <- "07[BNQTW][BP][04].*"
female_repro_intra <- "0U[1BJLNQSTW][01245679CF][04].*"
intraperitoneal_surg <- paste0("^(", gi_intra, "|", hepato_intra, "|",
                               lymph_intra, "|", female_repro_intra, ")")
# --- Lysis Surg Definition ---
gi_lys <- "(0D[N][678ABEFGHJLNUVW]|0D8W).*"
hepato_lys <- "0FN[012456789].*"
lymph_lys <- "07N[BP].*"
female_repro_lys <- "0UN[01245679CF].*"
lysis_proc <- paste0(
  "^(", gi_lys, "|",  hepato_lys, "|",
  lymph_lys, "|", female_repro_lys,")")
# --- Adhesion Related Diagnosis ---
adhesion_diag <- "^(K660|K565\\.?[012]|N994).*"
# --- Intervention Definition ---
surgical_intervention <- "^.{4}[04].*"