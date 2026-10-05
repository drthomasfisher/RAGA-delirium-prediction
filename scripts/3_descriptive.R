# Data Cleaning and Cohort Preparation

# 1. LIBRARIES & SETUP
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(naniar)
  library(flextable)
  library(tableone)
  library(xgboost)
  library(visdat)
})

# 2. LOAD & DEEP JOIN
message("Step 1: Loading Data...")
master_list <- readRDS(here::here("data_interim", "master_list.rds"))

extract_dfs <- function(x) {
  if (is.data.frame(x)) return(list(x))
  if (is.list(x)) return(purrr::map(x, extract_dfs) |> unlist(recursive = FALSE))
  list()
}

master <- purrr::reduce(extract_dfs(master_list), full_join, by = "id")

# --- Baseline Demographics & Comorbidities (T0) ---
t0_vars <- c(
  "age_years", "sex", "bmi", "asa_status", "fracture_type",
  "cardiovascular", "respiratory", "cns", "urinary_tract",
  "cog_pre_dementia_flag", "education_level", "ph_missing_flag"
)

# --- Pre-op Vitals, Labs & Cognitive Scores (T0) ---
lab_vars <- c(
  "temp_t0", "sysbp_t0", "diabp_t0", "hr_t0", "spo2_t0",
  "hb_t0", "albumin_t0", "creatinine_t0", "urea_t0",
  "na_t0", "k_t0", "cl_t0", "alt_t0", "ast_t0",
  "ph_t0", "po2_kpa_t0", "pco2_kpa_t0", "hco3_t0",
  "cog_mmse_pre", "urea_creat_ratio_t0"
)

# 3. SWIVEL AND CLEAN INTRA-OPERATIVE ANESTHETIC DRUGS & ADJUNCTS
message("Processing intra-operative drug logs...")

drug_cols <- names(master) %>% str_subset("^anaesthetic_drug_tidy_")

drugs_processed_wide <- master %>%
  select(id, all_of(drug_cols)) %>%
  pivot_longer(
    cols = all_of(drug_cols),
    names_to = "source_col",
    values_to = "raw_val"
  ) %>%
  filter(!is.na(raw_val)) %>%
  mutate(
    clean_drug = case_when(
      str_detect(tolower(raw_val), "propofol") ~ "intraop_propofol",
      str_detect(tolower(raw_val), "etomidate") ~ "intraop_etomidate",
      str_detect(tolower(raw_val), "sufentanil") ~ "intraop_sufentanil",
      str_detect(tolower(raw_val), "remifentanil") ~ "intraop_remifentanil",
      str_detect(tolower(raw_val), "fentanyl") & !str_detect(tolower(raw_val), "remi|sufen") ~ "intraop_fentanyl",
      str_detect(tolower(raw_val), "morphine") ~ "intraop_morphine",
      str_detect(tolower(raw_val), "butorphanol") ~ "intraop_butorphanol",
      str_detect(tolower(raw_val), "tramadol") ~ "intraop_tramadol",
      str_detect(tolower(raw_val), "hydromorphone") ~ "intraop_hydromorphone",
      str_detect(tolower(raw_val), "oxycodone") ~ "intraop_oxycodone",
      str_detect(tolower(raw_val), "ketorolac") ~ "intraop_ketorolac",
      str_detect(tolower(raw_val), "noyan|nuoyang") ~ "intraop_nalbuphine", 
      str_detect(tolower(raw_val), "lidocaine") ~ "intraop_lidocaine",
      str_detect(tolower(raw_val), "ropivacaine") ~ "intraop_ropivacaine",
      str_detect(tolower(raw_val), "levobupivacaine|levuvivivivaine|l-bub") ~ "intraop_levobupivacaine",
      str_detect(tolower(raw_val), "bupivacaine") & !str_detect(tolower(raw_val), "levo") ~ "intraop_bupivacaine",
      str_detect(tolower(raw_val), "chloroprocaine|clopcaine") ~ "intraop_chloroprocaine",
      str_detect(tolower(raw_val), "atracurium|atrakuammonium|tracurium") & !str_detect(tolower(raw_val), "homeopathic") ~ "intraop_atracurium",
      str_detect(tolower(raw_val), "rocuronium") ~ "intraop_rocuronium",
      str_detect(tolower(raw_val), "vecuronium|vicuronium") ~ "intraop_vecuronium",
      str_detect(tolower(raw_val), "atropine") ~ "intraop_atropine", 
      str_detect(tolower(raw_val), "adrenaline|epinephrine") ~ "intraop_adrenaline",
      str_detect(tolower(raw_val), "digoxin") ~ "intraop_digoxin",
      str_detect(tolower(raw_val), "nixon|nyxin") ~ "intraop_anisodamine_nyxin", 
      str_detect(tolower(raw_val), "dexamethasone") ~ "intraop_dexamethasone",
      str_detect(tolower(raw_val), "dexmedetomidine") ~ "intraop_dexmedetomidine",
      str_detect(tolower(raw_val), "midazolam|diazepine") ~ "intraop_benzodiazepine",
      str_detect(tolower(raw_val), "wanwen|voluven|starch") ~ "intraop_hydroxyethyl_starch", 
      TRUE ~ NA_character_ 
    )
  ) %>%
  filter(!is.na(clean_drug)) %>%
  distinct(id, clean_drug) %>% 
  mutate(present = 1L) %>%
  pivot_wider(
    names_from = clean_drug,
    values_from = present,
    values_fill = list(present = 0L)
  )

new_drug_vars <- setdiff(names(drugs_processed_wide), "id")

intraop_vars <- c(
  "anaes_tiva", "anaes_volatile", "anaes_regional",
  "intraop_hypotension", "hypotension_medicated",
  "estimated_blood_loss", "duration_hrs", "blood_transfusion",
  new_drug_vars 
)

# 4. COHORT FILTERING (REWRITTEN FOR COHORT-SPECIFIC TRACKING)
early_deaths <- master_list$outcomes$survival %>%
  filter(survival_time_days < 7) %>%
  pull(id)

cohort_filtered <- master %>%
  mutate(
    # Set cohort assignment immediately to allow split exclusion tracking
    cohort_assignment = if_else(
      grepl("^1", as.character(id)), 
      "Derivation", 
      "Validation"
    ),
    delirium_7d = case_when(
      cog_delirium_7d_binary %in% c(TRUE, 1, "1") ~ 1,  
      cog_delirium_7d_binary %in% c(FALSE, 0, "0") ~ 0, 
      TRUE ~ NA_real_
    ),
    delirium_preop_occurs = coalesce(delirium_t0 == 1,
                                     cog_pre_delirium_exists == TRUE,
                                     FALSE)
  )

# Process the filtering steps sequentially
n0 <- nrow(cohort_filtered)
c1 <- cohort_filtered %>% filter(delirium_preop_occurs == FALSE)
n1 <- nrow(c1)
c2 <- c1 %>% filter(!id %in% early_deaths)
n2 <- nrow(c2)
c3 <- c2 %>% filter(!is.na(delirium_7d))
n3 <- nrow(c3)

# Derivation cohort tracking
deriv_n0 <- cohort_filtered %>% filter(cohort_assignment == "Derivation") %>% nrow()
deriv_n1 <- c1 %>% filter(cohort_assignment == "Derivation") %>% nrow()
deriv_n2 <- c2 %>% filter(cohort_assignment == "Derivation") %>% nrow()
deriv_n3 <- c3 %>% filter(cohort_assignment == "Derivation") %>% nrow()

# Validation cohort tracking
valid_n0 <- cohort_filtered %>% filter(cohort_assignment == "Validation") %>% nrow()
valid_n1 <- c1 %>% filter(cohort_assignment == "Validation") %>% nrow()
valid_n2 <- c2 %>% filter(cohort_assignment == "Validation") %>% nrow()
valid_n3 <- c3 %>% filter(cohort_assignment == "Validation") %>% nrow()

# Build comprehensive split exclusion summary
exclusion_summary <- tibble(
  Step = c("Initial", "Baseline Delirium Excl", "Early Death Excl", "Missing Outcome Excl"),
  
  Total_Remaining = c(n0, n1, n2, n3),
  Total_Excluded  = c(0, n0 - n1, n1 - n2, n2 - n3),
  
  Derivation_Remaining = c(deriv_n0, deriv_n1, deriv_n2, deriv_n3),
  Derivation_Excluded  = c(0, deriv_n0 - deriv_n1, deriv_n1 - deriv_n2, deriv_n2 - deriv_n3),
  
  Validation_Remaining = c(valid_n0, valid_n1, valid_n2, valid_n3),
  Validation_Excluded  = c(0, valid_n0 - valid_n1, valid_n1 - valid_n2, valid_n2 - valid_n3)
)

print(exclusion_summary)

c3 <- c3 %>%
  mutate(
    age_years = ifelse(age_years < 50, NA_real_, age_years),
    sex = case_when(
      sex %in% c(0, "0", "M", "m", "Male", "male")      ~ "Male",
      sex %in% c(1, "1", "F", "f", "Female", "female")  ~ "Female",
      TRUE ~ NA_character_
    ),
    sex = factor(sex, levels = c("Male", "Female"))
  )

cohort_raw_descriptive <- c3

winsorise_vec <- function(x, ref = TRUE, lower = 0.01, upper = 0.99) {
  q <- quantile(x[ref], probs = c(lower, upper), na.rm = TRUE)
  x[x < q[1]] <- q[1]
  x[x > q[2]] <- q[2]
  x
}

# 5. CLEANING & PHYSIOLOGICAL STRUCTURES 
clean_binary_flag <- function(x) {
  case_when(
    x %in% c(TRUE, 1, "1", "Yes", "yes", "Y", "y") ~ 1L,
    x %in% c(FALSE, 0, "0", "No", "no", "N", "n")   ~ 0L,
    TRUE ~ NA_integer_
  )
}

cohort_cleaned <- c3 %>%
  mutate(across(any_of(c(
    "cardiovascular", "respiratory", "cns", "urinary_tract", 
    "cog_pre_dementia_flag", "anaes_tiva", "anaes_volatile", "anaes_regional"
  )), clean_binary_flag)) %>%
  
  mutate(cohort_assignment = if_else(
    grepl("^1", as.character(id)), 
    "Derivation", 
    "Validation"
  )) %>%
  mutate(fracture_type = case_when(
    fracture_type == "Femoral neck fracture" ~ "Intracapsular",
    fracture_type == "Fracture of femoral head" ~ "Intracapsular",
    fracture_type %in% c("Intertrochanteric", "Subtrochanteric") ~ "Extracapsular",
    TRUE ~ NA_character_
  )) %>%
  mutate(
    blood_transfusion = case_when(
      blood_transfusion == TRUE  ~ 1L,
      blood_transfusion == FALSE ~ 0L,
      TRUE ~ NA_integer_
    )
  ) %>%
  mutate(across(any_of(c(lab_vars, "age_years", "bmi", "asa_status")), 
                ~suppressWarnings(as.numeric(as.character(.))))) %>%
  mutate(
    age_years = if_else(age_years < 40 | age_years > 110, NA_real_, age_years)
  ) %>%
  mutate(across(any_of(c(lab_vars, "duration_hrs", "estimated_blood_loss", "bmi")),
                ~winsorise_vec(.x, cohort_assignment == "Derivation"))) %>%
  mutate(
    urea_creat_ratio_t0 = (urea_t0 / creatinine_t0) * 1000
  ) %>%
  mutate(across(any_of("urea_creat_ratio_t0"), ~winsorise_vec(.x, cohort_assignment == "Derivation"))) %>%
  mutate(ph_missing_flag = factor(if_else(is.na(ph_t0), "Missing", "Measured")))

cohort_cleaned_preimp <- cohort_cleaned

cohort_cleaned <- cohort_cleaned %>%
  left_join(drugs_processed_wide, by = "id") %>%
  mutate(across(any_of(new_drug_vars), ~ replace_na(.x, 0L)))
cohort_raw_unimputed <- cohort_cleaned

# 6. MISSINGNESS HEATMAP 
master_flat_raw <- readRDS(here::here("data_interim", "master_flat_raw.rds"))
heat_vars <- names(master_flat_raw)[names(master_flat_raw) %in% names(cohort_cleaned)]
if (length(heat_vars) < 5) {
  heat_vars <- names(master_flat_raw)[!names(master_flat_raw) %in% c("id", "centre")]
}

heat_df <- master_flat_raw %>%
  filter(id %in% cohort_raw_descriptive$id) %>%
  select(any_of(heat_vars))

p_heat <- visdat::vis_miss(heat_df, warn_large_data = FALSE) +
  theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1)) +
  labs(title = "Missingness Pattern (Pre-Imputation, Analytic Cohort)")

# 7. LAB THRESHOLDS & ABNORMAL OUTLIER FLAGGING
lab_thresholds <- tibble::tribble(
  ~variable,    ~lower, ~upper,
  "temp",       28,     48,
  "sysbp",      50,     300,
  "diabp",      30,     150,
  "hr",         20,     200,
  "spo2",       50,     100,
  "hb",         30,     210,
  "albumin",    10,     70,
  "creatinine", 10,     1500,
  "ph",         6.3,    7.7,
  "po2_kpa",    3,      80,
  "pco2_kpa",   2,      15,
  "hco3",       5,      50,
  "cl",         70,     130,
  "na",         80,     180,
  "k",          2,      10,
  "alt",        0,      100000,
  "ast",        0,      100000,
  "urea",       0.5,    100
)

target_labs <- lab_thresholds$variable
lab_pattern <- paste0("^(", paste(target_labs, collapse = "|"), ")_t0$")
lab_cols_to_flag <- grep(lab_pattern, names(c3), value = TRUE)

bio_flagged_readable <- c3 %>%
  select(id, all_of(lab_cols_to_flag)) %>%
  pivot_longer(-id, names_to = "var_raw", values_to = "value") %>%
  mutate(variable = str_replace(var_raw, "_t0$", "")) %>%
  left_join(lab_thresholds, by = "variable") %>%
  filter(!is.na(value) & (value < lower | value > upper)) %>%
  group_by(variable) %>%
  summarise(abnormal_values = paste0(id, " (Val: ", value, ")", collapse = "; "), .groups = "drop")

# 8. LAB-DELIRIUM RELATIONSHIP PLOTS
plot_labs <- c3 %>%
  select(delirium_7d, all_of(lab_cols_to_flag)) %>%
  pivot_longer(-delirium_7d, names_to = "lab", values_to = "value") %>%
  filter(!is.na(value)) %>%
  ggplot(aes(x = value, y = delirium_7d)) +
  geom_jitter(height = 0.05, alpha = 0.2, size = 0.5) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), color = "red", linetype = "dashed") +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), color = "blue") +
  coord_cartesian(ylim = c(0, 1)) +
  facet_wrap(~lab, scales = "free_x") +
  theme_minimal() +
  labs(title = "Lab Values vs. Delirium Risk", subtitle = "Red = Linear | Blue = Non-linear (GAM)", y = "Delirium Probability")

# 9. DESCRIPTIVE TABLES
t1_summary <- tableone::CreateTableOne(
  vars = t0_vars, 
  strata = "delirium_7d", 
  data = cohort_raw_unimputed
)

table1_flex <- as.data.frame(print(
  t1_summary, 
  printToggle = FALSE, 
  missing = TRUE,        
  showAllLevels = TRUE
)) %>%
  rownames_to_column("Variable") %>%
  flextable() %>%
  theme_booktabs() %>%
  autofit()

by_centre_table <- cohort_raw_unimputed %>%
  group_by(centre) %>%
  summarise(
    N = n(), 
    Incidence = mean(delirium_7d == 1, na.rm = TRUE), 
    Missing_Outcome_N = sum(is.na(delirium_7d)),
    .groups = "drop"
  ) %>%
  flextable() %>%
  autofit()

# 10. DERIVATION & VALIDATION COHORT SPLITTING
message("Splitting cohorts into Derivation and Validation sets...")

cohort_preop_df  <- cohort_raw_unimputed %>% select(-any_of(intraop_vars))
cohort_postop_df <- cohort_raw_unimputed

pre_data  <- cohort_preop_df %>% filter(cohort_assignment == "Derivation")
pre_val   <- cohort_preop_df %>% filter(cohort_assignment == "Validation")

post_data <- cohort_postop_df %>% filter(cohort_assignment == "Derivation")
post_val  <- cohort_postop_df %>% filter(cohort_assignment == "Validation")

core_data_subset <- pre_data %>% 
  select(id, centre, delirium_7d, age_years, urea_t0, cog_pre_dementia_flag)

core_val <- pre_val %>% 
  select(id, centre, delirium_7d, age_years, urea_t0, cog_pre_dementia_flag)


save(
  # Split Datasets 
  pre_data, 
  pre_val,
  post_data, 
  post_val,
  core_data_subset, 
  core_val,
  
  # Primary Modeling Track
  cohort_raw_unimputed,
  cohort_preop_df,   
  cohort_postop_df,  
  
  # Structural Frames additional
  cohort_raw_descriptive,          
  cohort_cleaned_preimp,           
  
  # Reference Indexes
  t0_vars, 
  lab_vars,
  intraop_vars,
  exclusion_summary, 
  table1_flex, 
  by_centre_table,
  bio_flagged_readable, 
  plot_labs, 
  p_heat,
  file = here::here("data_interim", "model_cohorts.RData")
)

message("Complete")

save(
  exclusion_summary, 
  table1_flex, 
  by_centre_table,
  bio_flagged_readable,
  file = here::here("data_interim", "results_descriptive.RData")
)

message("✅ Success: Descriptive validation structures preserved.")